#-*- coding: utf-8 -*-
'''
Name: MaxReport
Author: Tao Zhou, zhoutao@njmu.edu.cn
License:
This script could be used, modified and redistributed freely for academic research purpose.
Description:
MaxReport is an enhenced proteomic result reporting tool for MaxQuant.
Highlighted features:
Results are well organized and indexed.
Optimized and minimal protein reporting.
Generate site table of protein N-terminal modifications.
Provided general descriptive statistical analyses and figures.
Support isobaric labeling quantification at protein, peptide and site levels.
#change log (20160226):
>Supports "I=L" setting and other unknown amino acid in the sequence database.
>Fixes xlsxwriter bug for handling NAN and INF values.
>Supports MaxQuant 1.2.0.18 with no MS/MS IDs in the protein and site files.
>Detects MaxQuant version for bug report automatically.
>Improves multiple files selection.
'''

__version__='MaxReport 2.1 (2016-02-26)'

import os
import re
import argparse
import datetime
import numpy as np
import pygal
from pygal.style import RedBlueStyle
import webbrowser

#to solve coding problem with xlsxwriter
import xlsxwriter
import sys
reload(sys)
sys.setdefaultencoding('utf-8')

#reset svg background
RedBlueStyle.background='#FFFFFF'
RedBlueStyle.plot_background='#FFFFFF'

#functions for tsv files
def file2txt(fina):
    'read file to txt'
    tipf=open(fina,'rU')
    rst=tipf.read()
    tipf.close()
    return rst

def file2lis(fina):
    'read file lines to a list, will not retain empty lines'
    #use rU mode to solve line string of different system \n \r \r\n
    tipf=open(fina,'rU')
    rst=re.findall('.+',tipf.read())
    tipf.close()
    return rst

def file2tab(fina,etits=[],ktits=[],symb='\t'):
    'convert tab (or other symbols) separated text file to two-dimensional list data, based on file2lis()'
    rst=[]
    if etits==[] and ktits==[]:
        templis=file2lis(fina)
        for et in templis:
            rst.append(re.split(symb,et))
    else:
        tipf=open(fina,'rU')
        eline=tipf.readline()
        tmptits=re.split(symb,re.sub('\n','',eline))
        tmpidxlis=getlispos(etits+getkeystit(ktits,etits,tmptits),tmptits)
        while eline:
            rst.append(selectlis(re.split(symb,re.sub('\n','',eline)),tmpidxlis))
            eline=tipf.readline()
        tipf.close()
    return rst

def file2dic(fina,sn=0,symb='\t'):
    'convert tab (or other symbols) separated text to two-dimensional list data, based on file2lis(), sn: the key index'
    rst={}
    templis=file2lis(fina)
    keylis=[]
    for et in templis:
        tmprlis=re.split(symb,et)
        tmpkey=tmprlis[sn]
        if tmpkey not in keylis:
            rst[tmpkey]=tmprlis
            keylis.append(tmpkey)
    return rst

def unilis(alis):
    'remove redundancies in a list'
    rst=[]
    for ea in alis:
        if ea not in rst:
            rst.append(ea)
    return rst

def joinlis(alis,midstr='\t',endstr='\n'):
    'join a list with middle string (midstr) and end string (endstr)'
    tn=len(alis);jl='';i=0
    while i<tn:
        if i<>tn-1:
            jl+=str(alis[i])+midstr
        else:
            jl+=str(alis[i])+endstr
        i+=1
    return jl

def addlis(tlis,olis):
    'add elements of olis into tlis, non-redudant mode'
    for eo in olis:
        if eo not in tlis:
            tlis.append(eo)
    return tlis


def checklis(tlis,olis):
    'check the existence of tlis in olis'
    for et in tlis:
        if et not in olis:
            raise Exception('Configuration file is not complete, dictionary key: "%s" is lost!' %et)

def lis2file(alis,fina):
    'write alis to fina, one element per line'
    topf=open(fina,'w')
    topf.writelines(joinlis(alis,'\n','\n'))
    topf.close()

def tab2file(atab,fina,transp=0):
    'write atab (2D list) to tsd (tab seperated data) file (fina), transp: transpose'
    if transp==1:
        atab=map(list,zip(*atab))
    topf=open(fina,'w')
    for et in atab:
        topf.writelines(joinlis(et))
    topf.close()

def selectlis(tlis,ilis):
    'select elements of tlis according to index numbers in ilis, will not check redundancy'
    rst=[]
    tnu=len(tlis)
    for ei in ilis:
        if ei<=tnu-1 and ei>=0-tnu:
            rst.append(tlis[ei])
        else:
            raise IndexError('String index out of range')
    return rst

def dic2tab(adic):
    'convert a dict (adic) to a 2D list'
    atab=[]
    for ea in adic.keys():
        tmpval=adic[ea]
        if type(tmpval)==list:
            atab.append([ea]+tmpval)
        else:
            atab.append([ea,tmpval])
    return atab

def cpxmapping(acpx,adic,symb=';',csymb=';'):
    'mapping acpx based on adic'
    rst=[]
    for ec in cpx2lis(acpx,symb):
        rst.append(adic[ec])
    return joinlis(rst,csymb,'')

def unicpx(acpx,symb=';',csymb=';'):
    'remove redundancy in a cpx based on unilis()'
    tmplis=unilis(cpx2lis(acpx,symb))
    return joinlis(tmplis,csymb,'')

def selectcpx(acpx,slis=[0],symb=';',csymb=';'):
    'get the elements in acpx according to indexes in slis, based on selectlis()'
    rst=''
    rst=joinlis(selectlis(re.split(symb,acpx),slis),csymb,'')
    return rst

def getcpxpos(tcpx,ocpx,symb=';'):
    'get the positions of ocpx in tcpx'
    rst=[]
    tlis=cpx2lis(tcpx,symb)
    for ec in cpx2lis(ocpx,symb):
        rst.append(tlis.index(ec))
    return rst

def arr2cpxlis(atab,csymb=';'):
    'join a tab/array to a list of cpx'
    rst=[]
    tmparr=map(list,zip(*atab))
    for et in tmparr:
        rst.append(joinlis(et,csymb,''))
    return rst

def filtercpx(acpx,alis,fityp=1,symb=';',csymb=';'):
    'filtering acpx according to alis, fityp: 1(keep); 0(remove)'
    rstlis=[]
    tmplis=re.split(symb,acpx)
    for et in tmplis:
        if fityp==1:
            if et in alis:
                rstlis.append(et)
        else:
            if et not in alis:
                rstlis.append(et)
    return joinlis(rstlis,csymb,'')

def cpx2lis(acpx,symb='[;|,]'):
    'convert complex string to a list'
    return re.split(symb,acpx)

def getpos(adat,alis):
    'get the first position of adat in alis'
    return alis.index(adat)

def getlispos(olis,tlis):
    'get the first element positions of olis in tlis, based on getpos()'
    poslis=[]
    for eo in olis:
        poslis.append(getpos(eo,tlis))
    return poslis

def getkeytit(keystr,extlis,titlis):
    'find elements with key string (keystr) in titlis, extlis: excluded elements list'
    rst=[]
    for et in titlis:
        if keystr in et and (not et in extlis):
            rst.append(et)
    return rst

def gettitkey(keystr,extlis,titlis):
    'find the key strings in addition to keystr'
    rst=[]
    for et in titlis:
        if keystr in et and (not et in extlis):
            rst.append(re.sub(keystr,'',et))
    return rst

def getkeystit(keylis,extlis,titlis):
    'find elements with key strings in a list (substring), based on getkeytit()'
    rst=[]
    for et in titlis:
        for ek in keylis:
            if ek in et and (not et in extlis):
                rst.append(et)
    return rst

#absolute file and folder checking
def checkdir(adir):
    'check the existence and full path of a directory'
    if not os.path.isabs(adir):
        sys.stderr.write('Directory: "%s" is not a full path!\n' %adir)
        sys.exit(1)
    if not os.path.isdir(adir):
        sys.stderr.write('Directory: "%s" is not exist!\n' %adir)
        sys.exit(1)
    return adir

def checkfile(afile):
    'check the existence and full path of a file'
    if not os.path.isabs(afile):
        sys.stderr.write('File: "%s" is not a full path!\n' %afile)
        sys.exit(1)
    if not os.path.isfile(afile):
        sys.stderr.write('File: "%s" is not exist!\n' %afile)
        sys.exit(1)
    return afile

#logging function
def fmtime(dtobj='',dtpat='%m/%d/%Y, %H:%M:%S'):
    'dtpat: format date and time of datetime.datetime.now()'
    if not dtobj:
        dtobj=datetime.datetime.now()
    return dtobj.strftime(dtpat)

def loginfo(logcont='',isexit=0,msec=0):
    'report a log note'
    if msec==0:
        logtm=fmtime(dtpat='%x, %X')
    else:
        logtm=fmtime(dtpat='%x, %X, %f')
    sys.stdout.write('> '+logtm+'\n'+logcont+'\n')
    if isexit==1:
        sys.stdout.write('Exit.')
        sys.exit(0)
    return

#functions for MaxQuant report
def cpxmaxnum(acpx,symb=';'):
    'the count of element with max number'
    alis=[]
    for ec in re.split(symb,acpx):
        alis.append(int(ec))
    tmpmax=max(alis)
    maxnum=0
    for ea in alis:
        if ea==tmpmax:
            maxnum+=1
        else:
            break
    return maxnum

def getbestpgid(acpx,adic,bestsc='l'):
    'select best score id, bestsc: l(larger); s(lower)'
    tmplis=cpx2lis(acpx)
    rst=tmplis[0]
    for et in tmplis:
        if (bestsc=='l' and adic[et]>adic[rst]) or (bestsc=='s' and adic[et]<adic[rst]):
            rst=et
    return rst

def getobjkey(aval,adic):
    'get the dict (list dict) key of a val, !first key'
    rst=''
    for ea in adic.keys():
        if aval in adic[ea]:
            rst=ea
            break
    return rst

def countexp(aobjlis,aobjdic,alis,sn=1):
    'assign and count the appearence in different experiments, aobjlis(pg/pep/site, ms2 list), aobjdic(ms2: experiment), alis(experiment list)'
    rsttab=[['Experiments']]
    ctdic={}
    for ea in alis:
        ctdic[ea]=0
        rsttab[0].append(ea)
    for eo in aobjlis:
        tmplis=[eo[0]]
        for ea in alis:
            for ee in cpx2lis(eo[sn]):
                tmpstr=''
                if aobjdic[ee]==ea:
                    tmpstr='+'
                    break
            tmplis.append(tmpstr)
            if tmpstr=='+':
                ctdic[ea]+=1
        rsttab.append(tmplis)
    return rsttab,ctdic

#functions for configuration file
def readcfg(fina):
    'read configuration file to a dict {key, val}'
    rstdic={}
    tmplis=file2lis(fina)
    for et in tmplis:
        if et[0]!='#':
            etlis=re.split('\t',et)
            if len(etlis)>1:
                rstdic[etlis[0]]=etlis[1]
            else:
                rstdic[etlis[0]]=''
    return rstdic

def readcfgs(fina,clen=3):
    'read configuration file to a dict {key, [*val]}'
    rstdic={}
    tmplis=file2lis(fina)
    for et in tmplis:
        if et[0]!='#':
            etlis=re.split('\t',et)
            tmplen=len(etlis)
            if tmplen==1:
                rstdic[etlis[0]]=['']*clen
            else:
                rstdic[etlis[0]]=etlis[1:clen+1]+['']*(clen-len(etlis[1:]))
    return rstdic

#functions for fasta files
def file2seqs(fina,pats):
    'read a fasta file to a dict, {ID: sequence, header, [gene name, description]}'
    idpat=pats[0]
    anpats=pats[1:]
    seqdic={}
    idlis=[]
    tipf=open(fina,'r')
    while True:
        eline=tipf.readline()
        if not eline:
            break
        if eline[0]=='>':
            tmpid=re.findall(idpat,eline)[0]
            tmphd=re.findall('(.+)',re.sub('\t',' ',eline))[0]
            seqdic[tmpid]=['',tmphd,hdparse(tmphd,anpats)]
            if tmpid in idlis:
                raise Exception('Repeated identity')
            idlis.append(tmpid)
        else:
            tmpseq=re.sub('\s','',eline)
            seqdic[tmpid][0]+=tmpseq.upper()
    tipf.close()
    return seqdic

def formatseq(seqlis,linenu=100):
    'convert a sequence (list format: sequence, header) to string format for writing to fasta file'
    tmptxt=seqlis[1]+'\n'
    cn=0
    tmpseq=seqlis[0]
    seqlen=len(tmpseq)
    while cn<seqlen:
        tmptxt+=tmpseq[cn:cn+linenu]+'\n'
        cn+=linenu
    return tmptxt

def getpeploc(apep,aseq):
    'get the position of a peptide, !support ambiguous AA (X), also return previous and next AA'
    try:
        rst=aseq.index(apep)
    except Exception,em:
        rst=0
        isloca=0
        peplen=len(apep)
        prolen=len(aseq)
        while rst<=prolen-peplen:
            tmpseq=re.sub('[L|U|O|B|Z|J|X|*]','.',aseq[rst:rst+peplen])
            if re.findall(tmpseq,apep):
                isloca=1
                break
            rst+=1
        if isloca==0:
            raise Exception('Peptide (%s) not found in protein (%s)!' %(apep,aseq))
    if rst==0:
        paas='_'
    else:
        paas=aseq[rst-1]
    if rst+len(apep)==len(aseq):
        naas='_'
    else:
        naas=aseq[rst+len(apep)]
    return [rst+1,paas,naas]

def getseqwin(aseq,pos,flen=6):
    'get sequence window using a site position'
    siteidx=pos-1
    if siteidx-flen<0:
        rst='_'*(flen-len(aseq[0:siteidx]))+aseq[0:siteidx]
    else:
        rst=aseq[(siteidx-flen):siteidx]
    rst+=aseq[siteidx]
    if siteidx+flen+1>len(aseq):
        rst+=aseq[siteidx+1:len(aseq)]+'_'*(flen-len(aseq[siteidx+1:len(aseq)]))
    else:
        rst+=aseq[siteidx+1:siteidx+flen+1]
    return rst

def hdparse(hdstr,patlis):
    'get fasta header features based on pattern string list'
    rst=[]
    for ep in patlis:
        if not ep:
            rst.append('')
        else:
            tmplis=re.findall(ep,hdstr)
            if not tmplis:
                rst.append('')
            else:
                rst.append(tmplis[0])
    return rst

#functions for quantification
def floatlis(alis):
    'convert elements in alis to float type'
    rst=[]
    for ea in alis:
        if ea=='':
            rst.append(0)
        else:
            rst.append(float(ea))
    return rst

def intlis(alis):
    'convert elements in alis to int type'
    rst=[]
    for ea in alis:
        if ea=='':
            rst.append(0)
        else:
            rst.append(int(ea))
    return rst

def file2ary(fina):
    'read matrix file (2-d list),data type: np.array'
    rst=[]
    rawlis=file2lis(fina)
    for et in rawlis:
        tmplis=[]
        if et[0]!='#':
            etlis=re.split('\t',et)
            for ee in etlis:
                tmplis.append(float(ee))
            rst.append(tmplis)
    return np.array(rst)

def datcorrct(adic,corary):
    'correct adic list (1 D) in adat based on correction matrix, data type: np.array'
    rst={}
    for ea in adic.keys():
        rst[ea]=np.dot(np.linalg.matrix_power(corary*0.01,-1),adic[ea])
    return rst

def arysigma(alis,sigma=1.96):
    'descript summary of a dat array, data type: np.array'
    rst=[]
    tmpary=np.transpose(alis)
    for et in tmpary:
        tmpavg=np.mean(et)
        tmpstd=np.std(et)
        rst.append([tmpavg-tmpstd*sigma,tmpavg+tmpstd*sigma])
    return rst

def sigmarg(alis,valrgs):
    'check sigma range of alis'
    rst=[]
    outrg=0
    cnu=0
    for ea in alis:
        if ea<valrgs[cnu][0] or ea>valrgs[cnu][1]:
            outrg=1
            break
        cnu+=1
    return outrg

def statmtx(amtx):
    'calc median, mean, std values of a mtx'
    tmpary=np.transpose(amtx)
    tmpmeds=[]
    tmpavgs=[]
    tmpstds=[]
    for et in tmpary:
        tmpmeds.append(np.median(et))
        tmpavgs.append(np.mean(et))
        tmpstds.append(np.std(et))
    return [tmpmeds,tmpavgs,tmpstds]

def datcalibra(adic):
    'overall calibration: assures that the overall labeling becomes equal, data type: np.array'
    rst={}
    datlen=len(adic.values()[0])
    sumlis=np.array([0]*datlen)
    for ea in adic.keys():
        sumlis+=adic[ea]
    datavg=np.average(sumlis)
    clblis=datavg/sumlis
    for ea in adic.keys():
        rst[ea]=adic[ea]*clblis
    return (rst,clblis)

def thfilter(alis,trhold):
    'filter alis based on threshold value, will also remove empty values'
    outrg=0
    for ea in alis:
        if ea<trhold:
            outrg=1
            break
    if np.sum(alis)==0:
        outrg=1
    return outrg

def datlibra(adat,trhold=20,sigma=1.96):
    'for proteomics quantification based on modified libra algorithm, data type: np.array'
    #s1: normalization
    rst={}
    ctdic={}
    tmp1rst={}
    tmp2rst={}
    cidic={}
    for ea in adat.keys():
        tmplis=[]
        ctdic[ea]=[0,0]
        for ee in adat[ea]:
            if not thfilter(ee,trhold):
                tmplis.append(ee/np.sum(ee))
        tmp1rst[ea]=tmplis
        if len(tmplis)>0:
            cidic[ea]=arysigma(tmplis,sigma)
            ctdic[ea][0]=len(tmplis)
    #s2: remove outliers
    for ea in tmp1rst.keys():
        tmplis=[]
        for ee in tmp1rst[ea]:
            if not sigmarg(ee,cidic[ea]):
                tmplis.append(ee)
        tmp2rst[ea]=tmplis
        ctdic[ea][1]=len(tmplis)
    #s3: report quantification
    for ea in tmp2rst.keys():
        rst[ea]=[]
        if ctdic[ea][1]>0:
            rst[ea]=statmtx(tmp2rst[ea])
    return (rst,ctdic)

#functions for html parsing
def txt2link(ahref,atxt='',atar='_blank'):
    'creat a link for atxt'
    if atxt=='':
        atxt=ahref
    return '<a href="%s" target="%s">%s</a>'%(ahref,atar,atxt)

def addsvg(fina,atxt):
    'add a image in html'
    return '<embed type="image/svg+xml" src="%s" title="%s" pluginspage="http://www.adobe.com/svg/viewer/install/" />'%(fina,atxt)

#functions for svg ploting
def ctpepdis(atab):
    'Distribution of peptide counts'
    rst={}
    rglis=['1','2','3','4','5-10','11-20','>20']
    for ea in atab[0][1:]:
        rst[ea]={}
        for er in rglis:
            rst[ea][er]=0
    for ea in atab[1:]:
        cn=1
        for et in atab[0][1:]:
            if int(ea[cn])==1:
                rst[et]['1']+=1
            elif int(ea[cn])==2:
                rst[et]['2']+=1
            elif int(ea[cn])==3:
                rst[et]['3']+=1
            elif int(ea[cn])==4:
                rst[et]['4']+=1
            elif (10>=int(ea[cn])>=5):
                rst[et]['5-10']+=1
            elif (20>=int(ea[cn])>=11):
                rst[et]['11-20']+=1
            elif int(ea[cn])>20:
                rst[et]['>20']+=1
            cn+=1
    for ea in atab[0][1:]:
        for er in rglis:
            if rst[ea][er]==0:
                del rst[ea][er]
    return rst

def ctsitedis(atab):
    'distribution of site count'
    rst={}
    rglis=['1','2','3','4','5-10','>10']
    for er in rglis:
        rst[er]=0
    for ea in atab[1:]:
        if int(ea[1])==1:
            rst['1']+=1
        elif int(ea[1])==2:
            rst['2']+=1
        elif int(ea[1])==3:
            rst['3']+=1
        elif int(ea[1])==4:
            rst['4']+=1
        elif (10>=int(ea[1])>=5):
            rst['5-10']+=1
        elif int(ea[1])>10:
            rst['>10']+=1
    for er in rglis:
        if rst[er]==0:
            del rst[er]
    return rst

def tab2barsvg(atab,fina,xyaxes=['',''],transp=0,valtyp='i'):
    'plot svg figure using atab (must be classified); valtyp: i(int); f(float)'
    if transp==1:
        atab=map(list,zip(*atab))
    barsvg=pygal.Bar(pretty_print=True,x_title=xyaxes[0],y_title=xyaxes[1],width=640,height=360,explicit_size=True,style=RedBlueStyle)
    barsvg.x_labels=atab[0][1:]
    for ec in atab[1:]:
        if valtyp=='i':
            tmpvals=intlis(ec[1:])
        else:
            tmpvals=floatlis(ec[1:])
        barsvg.add(ec[0],tmpvals)
    barsvg.render_to_file(fina)

def dic2piesvg(adic,fina):
    'plot svg figure using adic'
    piesvg=pygal.Pie(pretty_print=True,width=450,height=360,explicit_size=True,style=RedBlueStyle)
    tmptab=sorted(dic2tab(adic),key=lambda es: es[1],reverse=True)
    for et in tmptab:
        piesvg.add(et[0],et[1])
    #for ea in adic.keys():
    #    piesvg.add(ea,adic[ea])
    piesvg.render_to_file(fina)

#functions for writing xlsx file
def xlsadd(wkbook,atab,sname='Sheet'):
    'write a table to a sheet'
    tmpst=wkbook.add_worksheet(sname)
    rnu=0
    for ea in atab:
        cnu=0
        for ee in ea:
            #fix xlsxwriter bug for NAN and INF values
            try:
                tmpst.write(rnu,cnu,ee)
            except Exception,em:
                tmpst.write(rnu,cnu,str(ee))
                #loginfo('xlsxwriter error: %s'%str(em))
            cnu+=1
        rnu+=1
    return

def tabs2xls(fina,stlis=[],tablis=[]):
    'write a list (*tablis) of tab to xls file, ?limitation of sheet name length'
    xlsfile=xlsxwriter.Workbook(fina)
    cn=0
    for et in tablis:
        tmpstnm='Sheet%s'%cn
        if cn<len(stlis):
            tmpstnm=stlis[cn]
        cn+=1
        xlsadd(xlsfile,rgfmt(et),tmpstnm)
    xlsfile.close()
    return

def rgfmt(atab):
    'format element type in atab for xls in a rude and guessing way?'
    rst=[]
    for ea in atab:
        tmplis=[]
        for eb in ea:
            tmpdat=eb
            try:
                tmpdat=int(eb)
            except Exception,em:
                try:
                    tmpdat=float(eb)
                except Exception,em:
                    pass
            tmplis.append(tmpdat)
        rst.append(tmplis)
    return rst

if __name__=='__main__':
    parser=argparse.ArgumentParser(description='MaxReport: enhanced proteomic result reporting tool for MaxQuant')
    #path issue: full path, ?chinese/specific coding
    #sdir: source directory
    parser.add_argument('sdir',help='Source folder for a MaxQuant project')
    #rcfg: parameter file for filtering
    parser.add_argument('rcfg',help='Configuration file for result reporting')
    #fasta rules configuration file
    parser.add_argument('-r','--frules',required=True,nargs='+',help='Rule files for parsing fasta sequences')
    #main protein sequence, further support multiple fasta files (change file2seqs results: ID header sequence, gene name, and description), merged 
    parser.add_argument('-s','--fseqs',required=True,nargs='+',help='Protein sequence files')
    #experimental design template file
    parser.add_argument('-e','--exptmp',help='File for experimental design')
    #peptide type for quantification
    parser.add_argument('-a','--allpep',action='store_true',help='Use all peptides for quantification')
    #correction file for isobaric quantification
    parser.add_argument('-c','--cormtx',help='Correction matrix file for isobaric quantification')
    #Minimum threshold value for isobaric quantification
    parser.add_argument('-m','--minith',default=5,help='Minimum threshold value for isobaric quantification')
    #Write combined identification and quantification results to xlsx files
    parser.add_argument('-x','--xlsxop',action='store_true',help='Write combined identification and quantification results to xlsx files')
    #print version
    parser.add_argument('-v','--version', action='version', version=__version__)
    args=parser.parse_args()
    #check input arguments
    sdir=checkdir(args.sdir)
    #obtain project name:
    prona=os.path.basename(sdir)
    rcfg=checkfile(args.rcfg)
    frules=[]
    for er in args.frules:
        frules.append(checkfile(er))
    fseqs=[]
    for es in args.fseqs:
        fseqs.append(checkfile(es))
    if args.exptmp:
        exptmp=checkfile(args.exptmp)
    allpep=args.allpep
    cormtx=np.array([])
    if args.cormtx:
        cormtx=file2ary(checkfile(args.cormtx))
    minith=float(args.minith)
    xlsxop=args.xlsxop
    #creat report folder
    rdir=os.path.join(sdir,'maxreport')
    if os.path.isdir(rdir):
        sys.stderr.write('Directory: "%s" is already exist!\n' %rdir)
        sys.exit(1)
    else:
        os.makedirs(rdir)
        os.makedirs(os.path.join(rdir,'report'))
        os.makedirs(os.path.join(rdir,'summary'))
        os.makedirs(os.path.join(rdir,'quant'))
        os.makedirs(os.path.join(rdir,'misc'))
        os.makedirs(os.path.join(rdir,'figures'))
    #Redirect output and error to a log file
    lgopf=open(os.path.join(sdir,'maxreport','log.txt'),'w')
    sys.stdout=lgopf
    sys.stderr=lgopf
    loginfo('Initial parameters checked. Working directory generated. Log file created.')
    #space issue
    mrcmds=[]
    for esa in sys.argv:
        if ' ' in esa:
            mrcmds.append('"%s"' %esa)
        else:
            mrcmds.append(esa)
    loginfo('Commands: %s'%joinlis(mrcmds,' ',''))
    #check keys
    cfglis=['sortid','pgids','pepids','pepid','msmsids','sitekeyids','rawfile','expdes','pgfixtits','pgkeys','proids','pepcts','unipepcts','revkey','conkey','proscore','bestsc','ms2fixtits','ms2keys','ms2rawfile','ms2scannu','pepfixtits','pepkeys','modfixtits','modkeys','modpros','modpos','modsite','modflkseq','modrplp','modrplpmin','nterms','ntermtar','pepseqtit','ntpeptits','qtionkey']
    frulelis=['prules']
    #read and check configuration file
    cfgdic=readcfg(rcfg)
    checklis(cfglis,cfgdic.keys())
    #reading fasta file
    fseqdic={}
    cn=0
    for efr in frules:
        fruledic=readcfgs(efr)
        checklis(frulelis,fruledic.keys())
        loginfo('Reading protein sequence file: [%s] ...'%fseqs[cn])
        fseqdic.update(file2seqs(fseqs[cn],fruledic['prules']))
        cn+=1
    loginfo('Protein sequences loaded.')
    #For id mapping and filtering among tables (msms, peptides, proteins and sites)
    sortid=cfgdic['sortid']
    pgids=cfgdic['pgids']
    pepids=cfgdic['pepids']
    pepid=cfgdic['pepid']
    msmsids=cfgdic['msmsids']
    sitekeyids=cfgdic['sitekeyids']
    #For extraction of experimental designs based on summary.txt
    rawfile=cfgdic['rawfile']
    expdes=cfgdic['expdes']
    #For proteinGroups reporting
    pgfixtits=cpx2lis(cfgdic['pgfixtits'])
    pgkeys=cpx2lis(cfgdic['pgkeys'])
    #For filtering and optimizing protein reporting
    proids=cfgdic['proids']
    pepcts=cfgdic['pepcts']
    unipepcts=cfgdic['unipepcts']
    revkey=cfgdic['revkey']
    conkey=cfgdic['conkey']
    proscore=cfgdic['proscore']
    bestsc=cfgdic['bestsc']
    #For msms reporting
    ms2fixtits=cpx2lis(cfgdic['ms2fixtits'])
    ms2keys=cpx2lis(cfgdic['ms2keys'])
    ms2rawfile=cfgdic['ms2rawfile']
    ms2scannu=cfgdic['ms2scannu']
    #For peptide reporting
    pepfixtits=cpx2lis(cfgdic['pepfixtits'])
    pepkeys=cpx2lis(cfgdic['pepkeys'])
    #For modification reporting
    modfixtits=cpx2lis(cfgdic['modfixtits'])
    modkeys=cpx2lis(cfgdic['modkeys'])
    #For modification optimizing
    modpros=cfgdic['modpros']
    modpos=cfgdic['modpos']
    modsite=cfgdic['modsite']
    modflkseq=int(cfgdic['modflkseq'])
    modrplp=cfgdic['modrplp']
    modrplpmin=float(cfgdic['modrplpmin'])
    #For protein N-term modification
    nterms=cpx2lis(cfgdic['nterms'])
    ntermtar=cfgdic['ntermtar']
    pepseqtit=cfgdic['pepseqtit']
    ntpeptits=cpx2lis(cfgdic['ntpeptits'])
    #For quantification based on reporter ions
    qtionkey=cfgdic['qtionkey']
    loginfo('Reporting configuration file and sequence parsing file checked.')
    #extraction of experimental designs from summary.txt or specified experimentalDesignTemplate.txt
    loginfo('Detecting experimental designs ...')
    if args.exptmp:
        osumtab=file2tab(exptmp)
        loginfo('Using experimental template file: %s' %exptmp)
        rawfile='Name'
        expdes='Experiment'
    else:
        osumtab=file2tab(os.path.join(sdir,'combined','txt','summary.txt'))
        loginfo('Using default summary.txt file to extract experimental design.')
    #experiment: raw files dict
    expdic={}
    if expdes in osumtab[0]:
        for eo in osumtab[1:]:
            try:
                tmpexp=eo[getpos(expdes,osumtab[0])]
                if tmpexp<>'':
                    if tmpexp not in expdic.keys():
                        expdic[tmpexp]=[eo[getpos(rawfile,osumtab[0])]]
                    else:
                        expdic[tmpexp].append(eo[getpos(rawfile,osumtab[0])])
            except Exception,em:
                loginfo(em)
    loginfo('%s experimental replications detected.' %len(expdic.keys()))
    del osumtab
    #check MaxQuant version for bug report
    mqver=re.findall('Version\t(.+)',file2txt(os.path.join(sdir,'combined','txt','parameters.txt')))[0]
    loginfo('MaxQuant version: %s'%mqver)
    #extract MS/MS ID mapping relationships based on evidence.txt
    loginfo('Extracting MS/MS ID mapping relationships based on evidence.txt ...')
    evpgms2dic={}
    evpepms2dic={}
    evsitesdic={}
    oevitab=file2tab(os.path.join(sdir,'combined','txt','evidence.txt'))
    sitecls=gettitkey(sitekeyids,[],oevitab[0])
    for esc in sitecls:
        evsitesdic[esc]={}
    for eo in oevitab[1:]:
        tmpmslis=re.split(';',eo[getpos(msmsids,oevitab[0])])
        for epg in re.split(';',eo[getpos(pgids,oevitab[0])]):
            if epg not in evpgms2dic.keys():
                evpgms2dic[epg]=[]
            for ems in tmpmslis:
                if ems not in evpgms2dic[epg]:
                    evpgms2dic[epg].append(ems)
        if eo[getpos(pepid,oevitab[0])] not in evpepms2dic.keys():
            evpepms2dic[eo[getpos(pepid,oevitab[0])]]=[]
        for ems in tmpmslis:
            if ems not in evpepms2dic[eo[getpos(pepid,oevitab[0])]]:
                evpepms2dic[eo[getpos(pepid,oevitab[0])]].append(ems)
        for esc in sitecls:
            for est in re.split(';',eo[getpos(esc+sitekeyids,oevitab[0])]):
                if est not in evsitesdic[esc].keys():
                    evsitesdic[esc][est]=[]
                for ems in tmpmslis:
                    if ems not in evsitesdic[esc][est]:
                        evsitesdic[esc][est].append(ems)
    del oevitab
    for epg in evpgms2dic.keys():
        tmplis=sorted(evpgms2dic[epg])
        evpgms2dic[epg]=joinlis(tmplis,';','')
    for epe in evpepms2dic.keys():
        tmplis=sorted(evpepms2dic[epe])
        evpepms2dic[epe]=joinlis(tmplis,';','')
    for esc in sitecls:
        for est in evsitesdic[esc].keys():
            tmplis=sorted(evsitesdic[esc][est])
            evsitesdic[esc][est]=joinlis(tmplis,';','')
    #initial filtering for proteinGroups: protein reporting optimize
    loginfo('Initial protein groups optimizing: \n(1). remove contaminants and reverse hits; \n(2). reporting minimum proteins for each protein group.')
    oprotab=file2tab(os.path.join(sdir,'combined','txt','proteinGroups.txt'))
    if len(oprotab)==0:
        sys.stderr.write('Empty proteinGroup file.')
        sys.exit(1)
    proidsidx=getpos(proids,oprotab[0])
    pepctsidx=getpos(pepcts,oprotab[0])
    unipepctsidx=getpos(unipepcts,oprotab[0])
    #[protein group id, msms ids] tab
    pgms2lis=[]
    #id: ms2 list, first protein ID, gene name, description
    pgantdic={}
    #id: score
    pgscdic={}
    #excluded id list
    pgrmlis=[]
    for eo in oprotab[1:]:
        tmpmnu=min([cpxmaxnum(eo[pepctsidx]),cpxmaxnum(eo[unipepctsidx])])
        tmpproids=selectcpx(eo[proidsidx],range(tmpmnu))
        #!and  ...: also to exclude empty MS/MS IDs for site mapping (MaxQuant errors: empty MS/MS IDs?)
        #if (revkey not in tmpproids) and (conkey not in tmpproids) and eo[getpos(msmsids,oprotab[0])]<>'':
        if (revkey not in tmpproids) and (conkey not in tmpproids) and evpgms2dic[eo[getpos(sortid,oprotab[0])]]<>'':
            pgms2lis.append([eo[getpos(sortid,oprotab[0])],evpgms2dic[eo[getpos(sortid,oprotab[0])]]])
            try:
                pgantdic[eo[getpos(sortid,oprotab[0])]]=[cpx2lis(evpgms2dic[eo[getpos(sortid,oprotab[0])]]),selectcpx(eo[proidsidx])]+fseqdic[selectcpx(eo[proidsidx])][-1]
            except Exception,em:
                loginfo('Initial filtering protein groups error ...\n%s\nPlease check the proteingroups.txt.'%str(em),1)
            pgscdic[eo[getpos(sortid,oprotab[0])]]=float(eo[getpos(proscore,oprotab[0])])
        else:
            pgrmlis.append(eo[getpos(sortid,oprotab[0])])
    del oprotab
    #msms filtering and generate optimized id list for msms, protein, peptide and protein IDs (caution: modification sites IDs may be missing due to MaxQuant errors?)
    loginfo('Filtering MS/MS scans ...')
    #oms2tab=file2tab(os.path.join(sdir,'combined','txt','msms.txt'))
    oms2tab=file2tab(os.path.join(sdir,'combined','txt','msms.txt'),[sortid,pepid,pgids]+ms2fixtits,ms2keys+[qtionkey,sitekeyids])
    ms2opf=open(os.path.join(sdir,'maxreport','report','msms.txt'),'w')
    #check quantification
    isquant=0
    if getkeytit(qtionkey,[],oms2tab[0]):
        isquant=len(getkeytit(qtionkey,[],oms2tab[0]))
    if isquant==0:
        loginfo('Found no reporter ions for quantification.')
    else:
        loginfo('Found %s reporter ions for quantification.'%isquant)
    ms2tits=getkeytit(sitekeyids,[],oms2tab[0])+ms2fixtits+getkeystit(ms2keys,ms2fixtits,oms2tab[0])+getkeytit(qtionkey,[],oms2tab[0])
    ms2titsidx=getlispos(ms2tits,oms2tab[0])
    ms2addtits=['id','Peptide ID','Protein group IDs','First proteins']+nterms
    ms2opf.writelines(joinlis(ms2addtits+ms2tits))
    #sitecls=gettitkey(sitekeyids,[],oms2tab[0])
    #nterm mod tit: id list (ms2, pep, pg)
    ntermms2dic={}
    ntermpepdic={}
    ntermpgdic={}
    for en in nterms:
        ntermms2dic[en]=[]
        ntermpepdic[en]=[]
        ntermpgdic[en]=[]
    #msms id: [reporter intensity array], using np.array()
    ms2rpions={}
    #msms id: unique proteingroup
    ms2unipg={}
    #msms id: experiment
    ms2expdic={}
    #msms id: [raw file, scan number]
    ms2idxlis={}
    for eo in oms2tab[1:]:
        pgscpx=filtercpx(eo[getpos(pgids,oms2tab[0])],pgrmlis,0)
        if pgscpx<>'':
            if len(expdic.keys())>0:
                ms2expdic[eo[getpos(sortid,oms2tab[0])]]=getobjkey(eo[getpos(ms2rawfile,oms2tab[0])],expdic)
            tmpunifet='YES'
            if len(cpx2lis(pgscpx))>1:
                tmpunifet='NO'
            ms2rpions[eo[getpos(sortid,oms2tab[0])]]=np.array(floatlis(selectlis(eo,getlispos(getkeytit(qtionkey,[],oms2tab[0]),oms2tab[0]))))
            ms2unipg[eo[getpos(sortid,oms2tab[0])]]=tmpunifet
            ms2idxlis[eo[getpos(sortid,oms2tab[0])]]=[eo[getpos(ms2rawfile,oms2tab[0])],eo[getpos(ms2scannu,oms2tab[0])]]
            tmpfps=[]
            for ep in cpx2lis(pgscpx):
                tmpfps.append(pgantdic[ep][1])
            tmpmodlis=cpx2lis(eo[getpos(ntermtar,oms2tab[0])])
            tmpntlis=[]
            for en in nterms:
                if en in tmpmodlis:
                    tmpntlis.append('+')
                    ntermms2dic[en].append(eo[getpos(sortid,oms2tab[0])])
                    ntermpepdic[en]=addlis(ntermpepdic[en],[eo[getpos(pepid,oms2tab[0])]])
                    ntermpgdic[en]=addlis(ntermpgdic[en],cpx2lis(pgscpx))
                else:
                    tmpntlis.append('')
            tmplis=[eo[getpos(sortid,oms2tab[0])],eo[getpos(pepid,oms2tab[0])],pgscpx,joinlis(tmpfps,';','')]+tmpntlis+selectlis(eo,ms2titsidx)
            ms2opf.writelines(joinlis(tmplis))
    tmpnterms=nterms[:]
    for en in tmpnterms:
        if len(ntermms2dic[en])==0:
            loginfo('Empty modification file: %s.'%(en+'Sites.txt'))
            nterms.remove(en)   
    del oms2tab
    ms2opf.close()
    loginfo('MS/MS scans reporting file created: %s' %('report/msms.txt'))
    #proteinGroups filtering
    loginfo('Filtering protein groups ...')
    #protein IDs for generate separated fasta sequences
    profilis=[]
    oprotab=file2tab(os.path.join(sdir,'combined','txt','proteinGroups.txt'))
    pgopf=open(os.path.join(sdir,'maxreport','report','proteinGroups.txt'),'w')
    pgtits=getkeytit(sitekeyids,[],oprotab[0])+pgfixtits+getkeystit(pgkeys,pgfixtits,oprotab[0])
    proaddtits=['id','Peptide IDs','MS/MS IDs','Protein accessions','All peptide counts','Unique peptide counts','First protein']+nterms
    pgopf.writelines(joinlis(proaddtits+pgtits+['Gene name','Description']))
    propeptab=[['Protein ID','Peptides','Unique peptides']]
    for eo in oprotab[1:]:
        if eo[getpos(sortid,oprotab[0])] not in pgrmlis:
            tmpmnu=min([cpxmaxnum(eo[pepctsidx]),cpxmaxnum(eo[unipepctsidx])])
            tmpproids=selectcpx(eo[proidsidx],range(tmpmnu))
            profilis=addlis(profilis,cpx2lis(tmpproids))
            tmplis=[eo[getpos(sortid,oprotab[0])],eo[getpos(pepids,oprotab[0])],evpgms2dic[eo[getpos(sortid,oprotab[0])]],tmpproids,selectcpx(eo[pepctsidx],range(tmpmnu)),selectcpx(eo[unipepctsidx],range(tmpmnu)),selectcpx(eo[proidsidx])]
            for en in nterms:
                if eo[getpos(sortid,oprotab[0])] in ntermpgdic[en]:
                    tmplis.append('+')
                else:
                    tmplis.append('')
            tmplis+=selectlis(eo,getlispos(pgtits,oprotab[0]))
            tmpproet=selectcpx(eo[proidsidx])
            tmplis+=pgantdic[eo[getpos(sortid,oprotab[0])]][-2:]
            pgopf.writelines(joinlis(tmplis))
            propeptab.append([tmpproet,selectcpx(eo[pepctsidx]),selectcpx(eo[unipepctsidx])])
    del oprotab
    pgopf.close()
    loginfo('Protein groups reporting file created: %s' %('report/proteinGroups.txt'))
    #peptides filtering
    loginfo('Filtering peptides ...')
    #nterm mod: peptide id: nterm peptide information
    ntpepsdic={}
    #id: ms2 list, best protein ID, gene namem description
    pepantdic={}
    for en in nterms:
        ntpepsdic[en]={}
    pepms2lis=[]
    unipepidlis=[]
    opeptab=file2tab(os.path.join(sdir,'combined','txt','peptides.txt'))
    pepopf=open(os.path.join(sdir,'maxreport','report','peptides.txt'),'w')
    peptits=getkeytit(sitekeyids,[],opeptab[0])+pepfixtits+getkeystit(pepkeys,pepfixtits,opeptab[0])
    pepaddtits=['id','Protein group IDs','MS/MS IDs','First proteins','Best protein group ID','Best protein']+nterms
    pepopf.writelines(joinlis(pepaddtits+peptits+['Position','Previous amino acid','Next amino acid','Gene name','Description']))
    for eo in opeptab[1:]:
        pgscpx=filtercpx(eo[getpos(pgids,opeptab[0])],pgrmlis,0)
        #!and  ...: also to exclude empty MS/MS IDs for site mapping (MaxQuant errors: empty MS/MS IDs?)
        if pgscpx<>'' and eo[getpos(msmsids,opeptab[0])]<>'':
            pepms2lis.append([eo[getpos(sortid,opeptab[0])],eo[getpos(msmsids,opeptab[0])]])
            tmpfps=[]
            if len(cpx2lis(pgscpx))==1:
                unipepidlis.append(eo[getpos(sortid,opeptab[0])])
            for ep in cpx2lis(pgscpx):
                tmpfps.append(pgantdic[ep][1])
            tmpbestpgid=getbestpgid(pgscpx,pgscdic,bestsc)
            tmpbestpro=pgantdic[tmpbestpgid][1]
            tmplis=[eo[getpos(sortid,opeptab[0])],pgscpx,eo[getpos(msmsids,opeptab[0])],joinlis(tmpfps,';',''),tmpbestpgid,tmpbestpro]
            for en in nterms:
                if eo[getpos(sortid,opeptab[0])] in ntermpepdic[en]:
                    tmplis.append('+')
                    ntpepseq=eo[getpos(pepseqtit,opeptab[0])]
                    ntpepsdic[en][eo[getpos(sortid,opeptab[0])]]=[ntpepseq[0],ntpepseq,eo[getpos(sortid,opeptab[0])]]+selectlis(eo,getlispos(ntpeptits,opeptab[0]))+[eo[getpos(msmsids,opeptab[0])],pgscpx,tmpbestpgid,tmpbestpro]
                else:
                    tmplis.append('')
            tmplis+=selectlis(eo,getlispos(peptits,opeptab[0]))
            tmpprocnt=fseqdic[tmpbestpro]
            try:
                tmplis+=getpeploc(eo[getpos(pepseqtit,opeptab[0])],tmpprocnt[0])
            except Exception,em:
                loginfo('Get peptide location error ...\n%s\nPlease check the peptides.txt file and sequence database.'%str(em),1)
            tmplis+=pgantdic[tmpbestpgid][-2:]
            pepopf.writelines(joinlis(tmplis))
            pepantdic[eo[getpos(sortid,opeptab[0])]]=[cpx2lis(eo[getpos(msmsids,opeptab[0])]),tmpbestpro]+pgantdic[tmpbestpgid][-2:]+[eo[getpos(pepseqtit,opeptab[0])]]
    del opeptab
    pepopf.close()
    loginfo('Peptides reporting file created: %s' %('report/peptides.txt'))
    #regenerating N-terminal modifications, based on proteingroups and peptides
    #nterm mod: protein IDs
    ntermprodic={}
    #nterm mod: pro-pos: site info
    ntsitedic={}
    #nterm mod: protein ID: site count
    prontsdic={}
    #nterm mod: [pro-pos, seqwin] list
    ntseqwindic={}
    #nterm mod: ms2 ids
    ntsitems2dic={}
    #nterm mod: id (pro-site): ms2 list, best protein id, gene name, description, position, site, sequence window
    ntsantdic={}
    for en in nterms:
        loginfo('Generating N-terminal modification: %s ...' %en)
        ntermprodic[en]=[]
        ntsitedic[en]={}
        ntsantdic[en]={}
        prontsdic[en]={}
        ntsitems2dic[en]=[]
        ntseqwindic[en]=[['Protein:Site position','Sequence window']]
        for entp in ntpepsdic[en].keys():
            entlis=ntpepsdic[en][entp]
            tmpproet=entlis[-1]
            if tmpproet not in ntermprodic[en]:
                ntermprodic[en].append(tmpproet)
            tmpprocnt=fseqdic[tmpproet]
            tmpanno=pgantdic[entlis[-2]][-2:]
            tmppeploc=getpeploc(entlis[1],tmpprocnt[0])[0]
            tmpseqwin=getseqwin(tmpprocnt[0],tmppeploc,modflkseq)
            #?rare case: ambiguous amino acid in sequence window
            if entlis[0]<>tmpseqwin[modflkseq]:
                tmpseqwin=tmpseqwin[0:modflkseq]+entlis[0]+tmpseqwin[modflkseq+1:]
            tmpsiteid=entlis[-1]+':'+str(tmppeploc)
            if tmpsiteid not in ntsitedic[en].keys():
                ntsitedic[en][tmpsiteid]=[entlis[-2],entlis[-1]]+tmpanno+[entlis[0],tmppeploc,tmpseqwin,[]]
            ntsitedic[en][tmpsiteid][-1].append(entlis[1:-2])
        tmpopf=open(os.path.join(sdir,'maxreport','report',en+'Sites.txt'),'w')
        ntermstits=['id','Best protein group ID','Best protein ID','Gene name','Description','Site','Position','Sequence window','Peptides','Peptide IDs']+ntpeptits+['MS/MS IDs','Protein group IDs','Unique (Groups)']
        tmpopf.writelines(joinlis(ntermstits))
        for enst in ntsitedic[en].keys():
            #tmpntslis: id, best pg id, best protein id, gene name, description, site, position, sequence window
            tmpntslis=ntsitedic[en][enst][0:-1]
            if tmpntslis[1] not in prontsdic[en].keys():
                prontsdic[en][tmpntslis[1]]=1
            else:
                prontsdic[en][tmpntslis[1]]+=1
            ntseqwindic[en].append([enst,tmpntslis[6]])
            #tmppeplis: pepseq, pep id, ... tits ..., ms2 ids, pg ids
            tmppeplis=arr2cpxlis(ntsitedic[en][enst][-1])
            tmppeplis[-1]=unicpx(tmppeplis[-1])
            tmppeplis[-2]=unicpx(tmppeplis[-2])
            tmpunipg='No'
            if len(cpx2lis(tmppeplis[-1]))==1:
                tmpunipg='Yes'
            tmpopf.writelines(joinlis([enst]+tmpntslis+tmppeplis+[tmpunipg]))
            ntsitems2dic[en].append([enst,tmppeplis[-2]])
            ntsantdic[en][enst]=[cpx2lis(tmppeplis[-2])]+tmpntslis[1:]
        tmpopf.close()
        loginfo('N-terminal modification: %s reporting file created: %s' %(en,'report/'+en+'Sites.txt'))
    #filtering sites (caution: MS/MS IDs may be missing due to low localization probability)
    #site tit: ms2 list
    stms2dic={}
    #site tit: peptide list
    stpepdic={}
    #site tit: protein ID list
    stprodic={}
    #site tit: proteingroup list
    stpgdic={}
    #site tit: site id
    stsitedic={}
    #site tit: protein: site count
    prositedic={}
    #site tit: prosite, seqwin
    siteseqwdic={}
    #site tit: site id, ms2 ids
    modsitems2dic={}
    #site tit: id : ms2 list, best protein id, gene name, description, position, site, sequence window
    stsantdic={}
    #to avoid remove element trap for a list in a loop
    tempcls=sitecls[:]
    for es in tempcls:
        loginfo('Filtering modification: %s ...' %es)
        stms2dic[es]=[]
        stpepdic[es]=[]
        stprodic[es]=[]
        stpgdic[es]=[]
        stsitedic[es]=[]
        modsitems2dic[es]=[]
        prositedic[es]={}
        stsantdic[es]={}
        siteseqwdic[es]=[['Protein:Site position','Sequence window']]
        ositetab=file2tab(os.path.join(sdir,'combined','txt',es+'Sites.txt'))
        if len(ositetab)==0:
            loginfo('Empty modification file: %s.'%(es+'Sites.txt'))
            sitecls.remove(es)
            continue
        tmpopf=open(os.path.join(sdir,'maxreport','report',es+'Sites.txt'),'w')
        modtits=modfixtits+getkeystit(modkeys,modfixtits,ositetab[0])
        modaddtits=['id','Peptide IDs','Protein group IDs','MS/MS IDs','First proteins','Positions','Best protein group ID','Best protein','Best protein position','Localization prob','Class']
        tmpopf.writelines(joinlis(modaddtits+modtits+['Sequence window','Gene name','Description','Unique (Groups)']))
        for eo in ositetab[1:]:
            pgscpx=filtercpx(eo[getpos(pgids,ositetab[0])],pgrmlis,0)
            if pgscpx<>'':
                #tmplis=[eo[getpos(sortid,ositetab[0])],eo[getpos(pepids,ositetab[0])],pgscpx,eo[getpos(msmsids,ositetab[0])]]
                tmplis=[eo[getpos(sortid,ositetab[0])],eo[getpos(pepids,ositetab[0])],pgscpx,evsitesdic[es][eo[getpos(sortid,ositetab[0])]]]
                pgscpxpos=getcpxpos(eo[getpos(pgids,ositetab[0])],pgscpx)
                tmpfps=selectcpx(eo[getpos(modpros,ositetab[0])],pgscpxpos)
                tmplis.append(tmpfps)
                tmpfpspos=selectcpx(eo[getpos(modpos,ositetab[0])],pgscpxpos)
                tmplis.append(tmpfpspos)
                tmpbestpgid=getbestpgid(pgscpx,pgscdic,bestsc)
                tmpbestpos=getcpxpos(pgscpx,tmpbestpgid)
                tmplis.append(tmpbestpgid)
                tmplis.append(selectcpx(tmpfps,tmpbestpos))
                tmplis.append(selectcpx(tmpfpspos,tmpbestpos))
                tmpmodpb=float(eo[getpos(modrplp,ositetab[0])])
                tmplis.append(eo[getpos(modrplp,ositetab[0])])
                if tmpmodpb>modrplpmin:
                    tmplis.append('Specific')
                else:
                    tmplis.append('Ambiguous')
                tmplis+=selectlis(eo,getlispos(modtits,ositetab[0]))
                tmpproet=selectcpx(tmpfps,tmpbestpos)
                tmpprocnt=fseqdic[tmpproet]
                tmpseqwin=getseqwin(tmpprocnt[0],int(selectcpx(tmpfpspos,tmpbestpos)),modflkseq)
                if eo[getpos(modsite,ositetab[0])]<>tmpseqwin[modflkseq]:
                    tmpseqwin=tmpseqwin[0:modflkseq]+eo[getpos(modsite,ositetab[0])]+tmpseqwin[modflkseq+1:]
                tmplis.append(tmpseqwin)
                tmplis+=pgantdic[tmpbestpgid][-2:]
                if tmpproet not in prositedic[es].keys():
                    prositedic[es][tmpproet]=1
                else:
                    prositedic[es][tmpproet]+=1
                stsitedic[es].append(eo[getpos(sortid,ositetab[0])])
                stprodic[es]=addlis(stprodic[es],[tmpproet])
                stpgdic[es]=addlis(stpgdic[es],cpx2lis(pgscpx))
                stpepdic[es]=addlis(stpepdic[es],cpx2lis(eo[getpos(pepids,ositetab[0])]))
                stms2dic[es]=addlis(stms2dic[es],cpx2lis(evsitesdic[es][eo[getpos(sortid,ositetab[0])]]))
                siteseqwdic[es].append([tmpproet+':'+selectcpx(tmpfpspos,tmpbestpos),tmpseqwin])
                tmpunipg='No'
                if len(cpx2lis(pgscpx))==1:
                    tmpunipg='Yes'
                tmpopf.writelines(joinlis(tmplis+[tmpunipg]))
                #!to exclude empty MS/MS IDs for site mapping (MaxQuant errors: empty MS/MS IDs?)
                if evsitesdic[es][eo[getpos(sortid,ositetab[0])]]<>'':
                    modsitems2dic[es].append([eo[getpos(sortid,ositetab[0])],evsitesdic[es][eo[getpos(sortid,ositetab[0])]]])
                    if tmpmodpb>modrplpmin:
                        stsantdic[es][eo[getpos(sortid,ositetab[0])]]=[cpx2lis(evsitesdic[es][eo[getpos(sortid,ositetab[0])]]),tmpproet]+pgantdic[tmpbestpgid][-2:]+[selectcpx(tmpfpspos,tmpbestpos),eo[getpos(modsite,ositetab[0])],tmpseqwin]
        del ositetab
        tmpopf.close()
        loginfo('Modification: %s reporting file created: %s' %(es,'report/'+es+'Sites.txt'))
    #summary and counting protein, peptides and sites
    sumfistab=[]
    tmpopf=open(os.path.join(sdir,'maxreport','summary','proteome-summary.txt'),'w')
    tmpopf.writelines(joinlis(['Class','Protein groups','Peptides','MS/MS spectra']))
    #tmpopf2=open(os.path.join(sdir,'maxreport','summary','proteome-experiments-summary.txt'),'w')
    #tmpopf2.writelines(joinlis(['Experiment','Protein groups','Peptides','MS/MS spectra']))
    tmpopf.writelines(joinlis(['Count',len(pgms2lis),len(pepms2lis),len(ms2idxlis.keys())]))
    tmpopf.close()
    loginfo('Proteome summary file created: %s' %('summary/proteome-summary.txt'))
    sumfistab.append(['Proteome summary',txt2link('summary/proteome-summary.txt')])
    tab2file(propeptab,os.path.join(sdir,'maxreport','summary','protein-peptides-count.txt'))
    loginfo('Protein-peptides count file created: %s' %('summary/protein-peptides-count.txt'))
    sumfistab.append(['Protein-peptides count',txt2link('summary/protein-peptides-count.txt')])
    expsumtab=[['Experiments']]
    if len(expdic.keys())>0:
        loginfo('Generating experiment and modification specific files ...')
        for ets in expdic.keys():
            expsumtab[0].append(ets)
        tmppgexpt,tmppgexps=countexp(pgms2lis,ms2expdic,expdic.keys())
        tab2file(tmppgexpt,os.path.join(sdir,'maxreport','summary','pg-exp-assignment.txt'))
        loginfo('Protein group-experiment assignment file created: %s' %('summary/pg-exp-assignment.txt'))
        sumfistab.append(['Protein group-experiment assignment',txt2link('summary/pg-exp-assignment.txt')])
        tmppepexpt,tmppepexps=countexp(pepms2lis,ms2expdic,expdic.keys())
        tab2file(tmppepexpt,os.path.join(sdir,'maxreport','summary','pep-exp-assignment.txt'))
        loginfo('Peptide-experiment assignment file created: %s' %('summary/pep-exp-assignment.txt'))
        sumfistab.append(['Peptide-experiment assignment',txt2link('summary/pep-exp-assignment.txt')])
        tmpms2expt,tmpms2exps=countexp(dic2tab(ms2expdic),ms2expdic,expdic.keys(),0)
        tab2file(tmpms2expt,os.path.join(sdir,'maxreport','summary','MS2-exp-assignment.txt'))
        loginfo('MS/MS scans-experiment assignment file created: %s' %('summary/MS2-exp-assignment.txt'))
        sumfistab.append(['MS/MS scans-experiment assignment',txt2link('summary/MS2-exp-assignment.txt')])
        #for ee in expdic.keys():
            #tmpopf2.writelines(joinlis([ee,tmppgexps[ee],tmppepexps[ee],tmpms2exps[ee]]))
        tmpexplis=['Protein groups']
        for ee in expdic.keys():
            tmpexplis.append(tmppgexps[ee])
        expsumtab.append(tmpexplis)
    #tmpopf2.close()
    tmpopf=open(os.path.join(sdir,'maxreport','summary','sites-summary.txt'),'w')
    tmpopf.writelines(joinlis(['Class','Sites','Peptides','Proteins','MS/MS spectra']))
    for en in nterms:
        if len(expdic.keys())>0:
            tmpntexpt,tmpntexps=countexp(ntsitems2dic[en],ms2expdic,expdic.keys())
            tab2file(tmpntexpt,os.path.join(sdir,'maxreport','summary',en+'-exp-assignment.txt'))
            loginfo('%s-experiment assignment file created: %s' %(en,'summary/'+en+'-exp-assignment.txt'))
            sumfistab.append(['%s site-experiment assignment'%en,txt2link('summary/%s-exp-assignment.txt'%en)])
            tmpsitelis=[en+' sites']
            for ets in expdic.keys():
                tmpsitelis.append(tmpntexps[ets])
            expsumtab.append(tmpsitelis)
        tmpopf.writelines(joinlis([en,len(ntsitedic[en].keys()),len(ntermpepdic[en]),len(ntermprodic[en]),len(ntermms2dic[en])]))
        tab2file([['Protein','Count']]+dic2tab(prontsdic[en]),os.path.join(sdir,'maxreport','summary','protein-'+en+'-sites-count.txt'))
        loginfo('Protein-%s sites count file created: %s' %(en,'summary/protein-'+en+'-sites-count.txt'))
        sumfistab.append(['Protein-%s sites count'%en,txt2link('summary/protein-%s-sites-count.txt'%en)])
        tab2file(ntseqwindic[en],os.path.join(sdir,'maxreport','summary',en+'-sequence-window.txt'))
        loginfo('%s sequence window file created: %s' %(en,'summary/'+en+'-sequence-window.txt'))
        sumfistab.append(['%s sequence window'%en,txt2link('summary/%s-sequence-window.txt'%en)])
    for es in sitecls:
        if len(expdic.keys())>0:
            tmpstexpt,tmpstexps=countexp(modsitems2dic[es],ms2expdic,expdic.keys())
            tab2file(tmpstexpt,os.path.join(sdir,'maxreport','summary',es+'-exp-assignment.txt'))
            loginfo('%s-experiment assignment file created: %s' %(es,'summary/'+es+'-exp-assignment.txt'))
            sumfistab.append(['%s site-experiment assignment'%es,txt2link('summary/%s-exp-assignment.txt'%es)])
            tmpsitelis=[es+' sites']
            for ets in expdic.keys():
                tmpsitelis.append(tmpstexps[ets])
            expsumtab.append(tmpsitelis)
        tmpopf.writelines(joinlis([es,len(stsitedic[es]),len(stpepdic[es]),len(stprodic[es]),len(stms2dic[es])]))
        tab2file([['Protein','Count']]+dic2tab(prositedic[es]),os.path.join(sdir,'maxreport','summary','protein-'+es+'-sites-count.txt'))
        loginfo('Protein-%s sites count file created: %s' %(es,'summary/protein-'+es+'-sites-count.txt'))
        sumfistab.append(['Protein-%s sites count'%es,txt2link('summary/protein-%s-sites-count.txt'%es)])
        tab2file(siteseqwdic[es],os.path.join(sdir,'maxreport','summary',es+'-sequence-window.txt'))
        loginfo('%s sequence window file created: %s' %(es,'summary/'+es+'-sequence-window.txt'))
        sumfistab.append(['%s sequence window'%es,txt2link('summary/%s-sequence-window.txt'%es)])
    tmpopf.close()
    loginfo('Site summary file created: %s' %('summary/sites-summary.txt'))
    sumfistab.append(['Site summary',txt2link('summary/sites-summary.txt')])
    if len(expdic.keys())>0:
        tab2file(expsumtab,os.path.join(sdir,'maxreport','summary','exp-specific-summary.txt'))
        loginfo('Experiment specific count file created: %s' %('summary/exp-specific-summary.txt'))
        sumfistab.append(['Experiment specific count',txt2link('summary/exp-specific-summary.txt')])
    #separate fasta files
    rseqopf=open(os.path.join(sdir,'maxreport','misc','mapped.fasta'),'w')
    #oseqopf=open(os.path.join(sdir,'maxreport','misc','others.fasta'),'w')
    for em in fseqdic.keys():
        tmpcont=formatseq(fseqdic[em])
        if em in profilis:
            rseqopf.writelines(tmpcont)
        #else:
        #    oseqopf.writelines(tmpcont)
    rseqopf.close()
    #oseqopf.close()
    loginfo('Mapped protein sequences file created: %s' %('misc/mapped.fasta'))
    #generate MS/MS IDs for further separation
    tab2file([[sortid,ms2rawfile,ms2scannu]]+dic2tab(ms2idxlis),os.path.join(sdir,'maxreport','misc','mapped-MS2-scans.txt'))
    loginfo('Mapped MS/MS scans file created: %s' %('misc/mapped-MS2-scans.txt'))
    if isquant:
        #preparation for quantification
        loginfo('Initial isobaric labeling based quantification ...')
        #posisible correction
        if args.cormtx:
            ms2rpions=datcorrct(ms2rpions,cormtx)
            loginfo('Tag impurity corrected.')
        #check experiments:
        ms2exprpdic={}
        expms2rpdic={}
        if len(expdic.keys())==0:
            for em in ms2idxlis.keys():
                ms2exprpdic[em]=''
            expms2rpdic['']=ms2idxlis.keys()
        else:
            ms2exprpdic=ms2expdic
            for em in ms2exprpdic.keys():
                if ms2exprpdic[em] not in expms2rpdic.keys():
                    expms2rpdic[ms2exprpdic[em]]=[em]
                else:
                    expms2rpdic[ms2exprpdic[em]].append(em)
        explis=expms2rpdic.keys()
        #overall calibration
        #exp: ms2: [reportor intensities ...]
        expms2rpions={}
        #exp: calibration factors
        expcalb={}
        for ee in explis:
            expms2rpions[ee]={}
            expcalb[ee]=[]
            for ems in expms2rpdic[ee]:
                expms2rpions[ee][ems]=ms2rpions[ems]
            expms2rpions[ee],expcalb[ee]=datcalibra(expms2rpions[ee])
        loginfo('Overall intensities calibrated.')
        #protein quantification
        quantctdic={}
        #exp: pgid: [[reportor intensities ...]]
        exppgrpions={}
        exppgctuni={}
        exppgctlbr={}
        exppgqt={}
        for ee in explis:
            exppgrpions[ee]={}
            exppgctuni[ee]={}
            exppgctlbr[ee]={}
            exppgqt[ee]={}
            for epg in pgantdic.keys():
                exppgrpions[ee][epg]=[]
                exppgctuni[ee][epg]=[0,0]
                for ems in pgantdic[epg][0]:
                    if ms2exprpdic[ems]==ee:
                        exppgctuni[ee][epg][0]+=1
                        if ms2unipg[ems]=='YES':
                            exppgctuni[ee][epg][1]+=1
                            exppgrpions[ee][epg].append(expms2rpions[ee][ems])
                        else:
                            if allpep:
                                exppgrpions[ee][epg].append(expms2rpions[ee][ems])
            exppgqt[ee],exppgctlbr[ee]=datlibra(exppgrpions[ee],minith,1.96)
        #write results
        rplen=len(ms2rpions.values()[0])
        rptits=[]
        for ee in explis:
            for esc in ['Total spectra','Unique spectra','Filtered spectra','Libra sepctra']:
                if ee:
                    rptits.append(esc+(' (%s)'%ee))
                else:
                    rptits.append(esc)
            for est in ['Median','Mean','STDEV']:
                for erp in range(rplen):
                    if ee:
                        rptits.append(est+(' reporter %d' %erp)+(' (%s)'%ee))
                    else:
                        rptits.append(est+(' reporter %d' %erp))
        tmppgiddic=file2dic(os.path.join(sdir,'maxreport','report','proteinGroups.txt'))
        tmpopf=open(os.path.join(sdir,'maxreport','quant','proteinGroups.txt'),'w')
        tmpopf.writelines(joinlis(['id','First protein ID','Gene name','Description']+rptits+tmppgiddic['id'][:-2]))
        quantctdic['Protein groups']=0
        for epg in pgantdic.keys():
            tmplis=[epg]+pgantdic[epg][1:]
            islibra=0
            for ee in explis:
                tmplis+=exppgctuni[ee][epg]+exppgctlbr[ee][epg]
                tmpqts=exppgqt[ee][epg]
                if tmpqts:
                    tmplis+=tmpqts[0]+tmpqts[1]+tmpqts[2]
                    islibra=1
                else:
                    tmplis+=['']*rplen+['']*rplen+['']*rplen
            quantctdic['Protein groups']+=islibra
            tmpopf.writelines(joinlis(tmplis+tmppgiddic[epg][:-2]))
        tmpopf.close()
        loginfo('Protein group quantification file created: %s' %('quant/proteinGroups.txt'))
        #write calibration factors
        tmpopf=open(os.path.join(sdir,'maxreport','quant','calibration-factors.txt'),'w')
        tmptits=['Reporter ions']
        for erp in range(rplen):
            tmptits.append('Reporter %d' %erp)
        tmpopf.writelines(joinlis(tmptits))
        for ee in explis:
            tmpexptit=ee
            if ee=='':
                tmpexptit='Factors'
            tmpopf.writelines(joinlis([tmpexptit]+list(expcalb[ee])))
        tmpopf.close()
        loginfo('Calibration factors file created: %s' %('quant/calibration-factors.txt'))
        #peptide quantification
        #exp: pepid: [[reportor intensities ...]]
        exppeprpions={}
        exppepctuni={}
        exppepctlbr={}
        exppepqt={}
        for ee in explis:
            exppeprpions[ee]={}
            exppepctuni[ee]={}
            exppepctlbr[ee]={}
            exppepqt[ee]={}
            for epe in pepantdic.keys():
                exppeprpions[ee][epe]=[]
                exppepctuni[ee][epe]=[0,0]
                for ems in pepantdic[epe][0]:
                    if ms2exprpdic[ems]==ee:
                        exppepctuni[ee][epe][0]+=1
                        if ms2unipg[ems]=='YES':
                            exppepctuni[ee][epe][1]+=1
                            exppeprpions[ee][epe].append(expms2rpions[ee][ems])
                        else:
                            if allpep:
                                exppeprpions[ee][epe].append(expms2rpions[ee][ems])
            exppepqt[ee],exppepctlbr[ee]=datlibra(exppeprpions[ee],minith,1.96)
        #write results
        #using the same rptits
        tmppepiddic=file2dic(os.path.join(sdir,'maxreport','report','peptides.txt'))
        tmpopf=open(os.path.join(sdir,'maxreport','quant','peptides.txt'),'w')
        tmpopf.writelines(joinlis(['id','Best protein ID','Gene name','Description','Peptide sequence']+rptits+tmppepiddic['id'][:-2]))
        quantctdic['Peptides']=0
        for epe in pepantdic.keys():
            tmplis=[epe]+pepantdic[epe][1:]
            islibra=0
            for ee in explis:
                tmplis+=exppepctuni[ee][epe]+exppepctlbr[ee][epe]
                tmpqts=exppepqt[ee][epe]
                if tmpqts:
                    tmplis+=tmpqts[0]+tmpqts[1]+tmpqts[2]
                    islibra=1
                else:
                    tmplis+=['']*rplen+['']*rplen+['']*rplen
            quantctdic['Peptides']+=islibra
            tmpopf.writelines(joinlis(tmplis+tmppepiddic[epe][:-2]))
        tmpopf.close()
        loginfo('Peptides quantification file created: %s' %('quant/peptides.txt'))
        #N-term mods quantification
        for en in nterms:
            expntsrpions={}
            expntsctuni={}
            expntsctlbr={}
            expntsqt={}
            for ee in explis:
                expntsrpions[ee]={}
                expntsctuni[ee]={}
                expntsctlbr[ee]={}
                expntsqt[ee]={}
                for ent in ntsantdic[en].keys():
                    expntsrpions[ee][ent]=[]
                    expntsctuni[ee][ent]=[0,0]
                    for ems in ntsantdic[en][ent][0]:
                        if ms2exprpdic[ems]==ee:
                            expntsctuni[ee][ent][0]+=1
                            if ms2unipg[ems]=='YES':
                                expntsctuni[ee][ent][1]+=1
                                expntsrpions[ee][ent].append(expms2rpions[ee][ems])
                            else:
                                if allpep:
                                    expntsrpions[ee][ent].append(expms2rpions[ee][ems])
                expntsqt[ee],expntsctlbr[ee]=datlibra(expntsrpions[ee],minith,1.96)
            #write results
            #using the same rptits
            tmpntiddic=file2dic(os.path.join(sdir,'maxreport','report',en+'Sites.txt'))
            tmpopf=open(os.path.join(sdir,'maxreport','quant',en+'Sites.txt'),'w')
            tmpopf.writelines(joinlis(['id','Best protein ID','Gene name','Description','Position','Site','Sequence window']+rptits+tmpntiddic['id'][0:2]+tmpntiddic['id'][8:]))
            quantctdic[en]=0
            for ent in ntsantdic[en].keys():
                tmplis=[ent]+ntsantdic[en][ent][1:]
                islibra=0
                for ee in explis:
                    tmplis+=expntsctuni[ee][ent]+expntsctlbr[ee][ent]
                    tmpqts=expntsqt[ee][ent]
                    if tmpqts:
                        tmplis+=tmpqts[0]+tmpqts[1]+tmpqts[2]
                        islibra=1
                    else:
                        tmplis+=['']*rplen+['']*rplen+['']*rplen
                quantctdic[en]+=islibra
                tmpopf.writelines(joinlis(tmplis+tmpntiddic[ent][0:2]+tmpntiddic[ent][8:]))
            tmpopf.close()
            loginfo('%s quantification file created: %s' %(en,'quant/'+en+'Sites.txt'))
        #Sites quantification
        for es in sitecls:
            expstsrpions={}
            expstsctuni={}
            expstsctlbr={}
            expstsqt={}
            for ee in explis:
                expstsrpions[ee]={}
                expstsctuni[ee]={}
                expstsctlbr[ee]={}
                expstsqt[ee]={}
                for ens in stsantdic[es].keys():
                    expstsrpions[ee][ens]=[]
                    expstsctuni[ee][ens]=[0,0]
                    for ems in stsantdic[es][ens][0]:
                        if ms2exprpdic[ems]==ee:
                            expstsctuni[ee][ens][0]+=1
                            if ms2unipg[ems]=='YES':
                                expstsctuni[ee][ens][1]+=1
                                expstsrpions[ee][ens].append(expms2rpions[ee][ems])
                            else:
                                if allpep:
                                    expstsrpions[ee][ens].append(expms2rpions[ee][ems])
                expstsqt[ee],expstsctlbr[ee]=datlibra(expstsrpions[ee],minith,1.96)
            #write results
            #using the same rptits
            tmpstiddic=file2dic(os.path.join(sdir,'maxreport','report',es+'Sites.txt'))
            tmpopf=open(os.path.join(sdir,'maxreport','quant',es+'Sites.txt'),'w')
            tmpopf.writelines(joinlis(['id','Best protein ID','Gene name','Description','Position','Site','Sequence window']+rptits+tmpstiddic['id'][:-3]+[tmpstiddic['id'][-1]]))
            quantctdic[es]=0
            for ens in stsantdic[es].keys():
                tmplis=[ens]+stsantdic[es][ens][1:]
                islibra=0
                for ee in explis:
                    tmplis+=expstsctuni[ee][ens]+expstsctlbr[ee][ens]
                    tmpqts=expstsqt[ee][ens]
                    if tmpqts:
                        tmplis+=tmpqts[0]+tmpqts[1]+tmpqts[2]
                        islibra=1
                    else:
                        tmplis+=['']*rplen+['']*rplen+['']*rplen
                quantctdic[es]+=islibra
                tmpopf.writelines(joinlis(tmplis+tmpstiddic[ens][:-3]+[tmpstiddic[ens][-1]]))
            tmpopf.close()
            loginfo('%s quantification file created: %s' %(es,'quant/'+es+'Sites.txt'))
        #write quantification summary
        tab2file([['Class','Count']]+dic2tab(quantctdic),os.path.join(sdir,'maxreport','quant','quantification-summary.txt'),1)
        loginfo('Quantification summary file created: %s' %('quant/quantification-summary.txt'))
    #write combined xlsx files
    cbdident=''
    cbdquant=''
    if xlsxop:
        tmpstnames=['proteinGroups','peptides']
        for en in nterms:
            tmpstnames.append(en+'Sites')
        for es in sitecls:
            tmpstnames.append(es+'Sites')
        #tmpstnames.append('msms')
        tmpsttabs=[]
        for est in tmpstnames:
            try:
                tmpsttabs.append(file2tab(os.path.join(sdir,'maxreport','report',est+'.txt')))
            except Exception,em:
                sys.stderr.write('Too large to add to a xlsx file: %s\n'%os.path.join(sdir,'maxreport','report',est+'.txt'))
        tabs2xls(os.path.join(sdir,'maxreport',prona+'-identification.xlsx'),tmpstnames,tmpsttabs)
        loginfo('Combined identification xlsx file created: %s' %(prona+'-identification.xlsx'))
        cbdident='''                    <tr>
                      <td width="30%%" align="left" valign="top">Combined identification file</td>
                      <td align="left" valign="top">%s</td>
                    </tr>'''%txt2link(prona+'-identification.xlsx')
        if isquant:
            tmpstnames=['proteinGroups','peptides']
            for en in nterms:
                tmpstnames.append(en+'Sites')
            for es in sitecls:
                tmpstnames.append(es+'Sites')
            tmpsttabs=[]
            for est in tmpstnames:
                tmpsttabs.append(file2tab(os.path.join(sdir,'maxreport','quant',est+'.txt')))
            tabs2xls(os.path.join(sdir,'maxreport',prona+'-quantification.xlsx'),tmpstnames,tmpsttabs)
            loginfo('Combined quantification xlsx file created: %s' %(prona+'-quantification.xlsx'))
            cbdquant='''                    <tr>
                      <td width="30%%" align="left" valign="top">Combined quantification file</td>
                      <td align="left" valign="top">%s</td>
                    </tr>'''%txt2link(prona+'-quantification.xlsx')
    #generate html index file
    htmlis=[fmtime()+'<br>MaxQuant version: %s'%mqver]
    htmlis.append(joinlis(mrcmds,'&nbsp;',''))
    htmlis.append(cbdident)
    htmlis.append(txt2link('report/proteinGroups.txt'))
    htmlis.append(txt2link('report/peptides.txt'))
    tmpntlinks=[]
    for en in nterms:
        tmpntlinks.append(txt2link('report/'+en+'Sites.txt'))
    htmlis.append(joinlis(tmpntlinks,';&nbsp;',''))
    tmpstlinks=[]
    for es in sitecls:
        tmpstlinks.append(txt2link('report/'+es+'Sites.txt'))
    htmlis.append(joinlis(tmpstlinks,';&nbsp;',''))
    htmlis.append(txt2link('report/msms.txt'))
    tab2barsvg(file2tab(os.path.join(sdir,'maxreport','summary','proteome-summary.txt')),os.path.join(sdir,'maxreport','figures','proteome-summary.svg'))
    loginfo('Figure 3.1 (proteome summary) created: %s' %('figures/proteome-summary.svg'))
    htmlis.append(addsvg('figures/proteome-summary.svg','Proteome summary'))
    propepcts=ctpepdis(propeptab)
    dic2piesvg(propepcts['Peptides'],os.path.join(sdir,'maxreport','figures','pro-all-pep-count.svg'))
    loginfo('Figure 3.2 (protein-all peptides count) created: %s' %('figures/pro-all-pep-count.svg'))
    htmlis.append(addsvg('figures/pro-all-pep-count.svg','Protein-all peptides count'))
    dic2piesvg(propepcts['Unique peptides'],os.path.join(sdir,'maxreport','figures','pro-uni-pep-count.svg'))
    loginfo('Figure 3.3 (protein-unique peptides count) created: %s' %('figures/pro-uni-pep-count.svg'))
    htmlis.append(addsvg('figures/pro-uni-pep-count.svg','Protein-unique peptides count'))
    tab2barsvg(file2tab(os.path.join(sdir,'maxreport','summary','sites-summary.txt')),os.path.join(sdir,'maxreport','figures','sites-summary.svg'))
    loginfo('Figure 3.4 (sites summary) created: %s' %('figures/sites-summary.svg'))
    htmlis.append(addsvg('figures/sites-summary.svg','Sites summary'))
    cn=5
    prositehtm=''
    prositetmp='''                <tr>
                  <td align="center" valign="middle">%s</td>
                </tr>
                <tr>
                  <td align="center" valign="middle"><strong>Figure 3.%s Distribution of %s site counts</strong></td>
                </tr>
'''
    for en in nterms:
        dic2piesvg(ctsitedis([['Protein','Count']]+dic2tab(prontsdic[en])),os.path.join(sdir,'maxreport','figures','pro-%s-count.svg'%en))
        loginfo('Figure 3.%s (protein-%s site count) created: %s' %(cn,en,'figures/pro-%s-count.svg'%en))
        prositehtm+=prositetmp %(addsvg('figures/pro-%s-count.svg'%en,'Protein-%s sites count'%en),cn,en)
        cn+=1
    for es in sitecls:
        dic2piesvg(ctsitedis([['Protein','Count']]+dic2tab(prositedic[es])),os.path.join(sdir,'maxreport','figures','pro-%s-count.svg'%es))
        loginfo('Figure 3.%s (protein-%s site count) created: %s' %(cn,es,'figures/pro-%s-count.svg'%es))
        prositehtm+=prositetmp %(addsvg('figures/pro-%s-count.svg'%es,'Protein-%s sites count'%es),cn,es)
        cn+=1
    htmlis.append(prositehtm)
    expsumhtml=''
    if len(expdic.keys())>0:
        tab2barsvg(expsumtab,os.path.join(sdir,'maxreport','figures','exps-summary.svg'))
        loginfo('Figure 3.%s (experiment specific summary) created: %s' %(cn,'figures/exps-summary.svg'))
        expsumhtml='''<table width="95%%" border="0" align="center" cellpadding="0" cellspacing="10">
                <tr>
                  <td align="left" valign="top"><strong>3.5 Experimental replications assignment and count</strong></td>
                </tr>
                <tr>
                  <td align="center" valign="middle">%s</td>
                </tr>
                <tr>
                  <td align="center" valign="middle"><strong>Figure 3.%s Summary of experiment specific identification</strong></td>
                </tr>
                </table>''' %(addsvg('figures/exps-summary.svg','Experiment specific summary'),cn)
    htmlis.append(expsumhtml)
    sumfishtml=''
    for esf in sumfistab[:-1]:
        sumfishtml+='''                    <tr>
                      <td align="left" valign="top">%s</td>
                      <td align="left" valign="top">%s</td>
                      </tr>
'''%(esf[0],esf[1])
    sumfishtml+='''                    <tr>
                      <td align="left" valign="top" class="tltbbt">%s</td>
                      <td align="left" valign="top" class="tltbbt">%s</td>
                      </tr>'''%(sumfistab[-1][0],sumfistab[-1][1])
    htmlis.append(sumfishtml)
    if isquant:
        tmpntlinks=[]
        for en in nterms:
            tmpntlinks.append(txt2link('quant/'+en+'Sites.txt'))
        tmpstlinks=[]
        for es in sitecls:
            tmpstlinks.append(txt2link('quant/'+es+'Sites.txt'))
        tab2barsvg(file2tab(os.path.join(sdir,'maxreport','quant','calibration-factors.txt')),os.path.join(sdir,'maxreport','figures','calibration-factors.svg'),valtyp='f')
        loginfo('Figure 4.1 (calibration factors) created: %s' %('figures/calibration-factors.svg'))
        tab2barsvg(file2tab(os.path.join(sdir,'maxreport','quant','quantification-summary.txt')),os.path.join(sdir,'maxreport','figures','quantification-summary.svg'))
        loginfo('Figure 4.2 (quantification-summary) created: %s' %('figures/quantification-summary.svg'))
        quanthtml='''<table width="95%%" border="0" cellspacing="10" cellpadding="0">
              <tr>
                <td align="center" valign="top">%s</td>
              </tr>
              <tr>
                <td align="center" valign="top"><strong>Figure 4.1 Overall calibration factors</strong></td>
              </tr>
              <tr>
                <td align="center" valign="top">%s</td>
              </tr>
              <tr>
                <td align="center" valign="top"><strong>Figure 4.2 Summary of quantification results</strong></td>
              </tr>
              <tr>
                <td align="center" valign="top"><strong>Table 4.1 Location of quantification result files</strong></td>
              </tr>
              <tr>
                <td align="center" valign="top"><table width="100%%" border="0" align="center" cellpadding="0" cellspacing="0">
                  <tr>
                    <td width="30%%" align="left" valign="top" class="tltbhd">Description</td>
                    <td align="left" valign="top" class="tltbhd">File location</td>
                  </tr>
                  %s
                  <tr>
                    <td width="30%%" align="left" valign="top">Summary count</td>
                    <td align="left" valign="top">%s</td>
                  </tr>
                  <tr>
                    <td align="left" valign="top">Calibration factors</td>
                    <td align="left" valign="top">%s</td>
                  </tr>
                  <tr>
                    <td width="30%%" align="left" valign="top">Protein groups</td>
                    <td align="left" valign="top">%s</td>
                  </tr>
                  <tr>
                    <td align="left" valign="top">Peptides</td>
                    <td align="left" valign="top">%s</td>
                  </tr>
                  <tr>
                    <td align="left" valign="top">N-terminal modifications</td>
                    <td align="left" valign="top">%s</td>
                  </tr>
                  <tr>
                    <td align="left" valign="top" class="tltbbt">Modifications</td>
                    <td align="left" valign="top" class="tltbbt">%s</td>
                  </tr>
                </table></td>
              </tr>
            </table>'''% (addsvg('figures/calibration-factors.svg','Calibration factors'),addsvg('figures/quantification-summary.svg','Quantification-summary')\
                          ,cbdquant\
                          ,txt2link('quant/quantification-summary.txt')\
                          ,txt2link('quant/calibration-factors.txt')\
                          ,txt2link('quant/proteinGroups.txt')\
                          ,txt2link('quant/peptides.txt')\
                          ,joinlis(tmpntlinks,';&nbsp;','')\
                          ,joinlis(tmpstlinks,';&nbsp;','')\
                          )
    else:
        quanthtml='''<table width="95%" border="0" cellspacing="10" cellpadding="0">
                <tr>
                  <td align="center" valign="middle">None</td>
                </tr>
              </table>'''
    htmlis.append(quanthtml)
    htmlis.append(txt2link('misc/mapped.fasta'))
    htmlis.append(txt2link('misc/mapped-MS2-scans.txt'))
    htmlis.append(__version__)
    htmltmp='''<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<title>MaxReport - Enhanced proteomic result reporting tool for MaxQuant</title>
<style type="text/css">
body {
	margin-left: 0px;
	margin-top: 0px;
	margin-right: 0px;
	margin-bottom: 0px;
	background-repeat: repeat-x;
	background-color: #399;
}
a:link {
	color: #CC0099;
	text-decoration: none;
}
a:visited {
	text-decoration: none;
	color:#F09;
}
a:hover {
	position: relative; left:1px; top:1px; color: #F09;text-decoration:underline;
}
a:active {
	text-decoration: none;
	color:#F09;
}
.tit {
font-size: 24px;
font-weight: bold;
color: #000;
text-transform: capitalize;
} 
.subtit {
font-size: 18px;
font-weight: bold;
color: #333;
}
.ctbg {background-color: #CFF;}
.solidbox {background-color: #CFF;
BORDER-RIGHT:#099 1px solid; 
BORDER-TOP:#099 1px solid; 
BORDER-LEFT:#099 1px solid; 
BORDER-BOTTOM:#099 1px solid;	
}
.tibg {background-color: #099;
color:#FFF;
font-weight: bold;
text-align: left;
}
.dashedbox {
	background-color: #FC9;
BORDER-RIGHT: #FC9 2px solid; 
BORDER-TOP: #FC9 2px solid; 
BORDER-LEFT: #FC9 2px solid; 
BORDER-BOTTOM: #FC9 2px solid;
}
.tltbhd
{
font-weight: bold;
BORDER-TOP:#000000 1px solid; 
BORDER-BOTTOM:#000000 1px solid;
	}
.tltbbt
{BORDER-BOTTOM: #000000 1px solid;	}
body,td,th {
	font-family: "Times New Roman", Times, serif;
	font-size: 18px;
	color: #000;
}
</style>
      <script language="javascript"> 
function dh(TagName){
	var obj = document.getElementById(TagName);
	if(obj.style.display==""){
		obj.style.display = "none";
	}else{
		obj.style.display = "";
	}
}
      </script>
</head>
<body>
<table width="969" border="0" align="center" cellpadding="0" cellspacing="0" bgcolor="#66CCCC" style="word-break:break-all">
  <tr>
    <td width="969" align="center" valign="top"><table width="95%%" border="0" cellspacing="0" cellpadding="0">
      <tr>
        <td width="180" height="120" align="center" valign="middle">
          <table width="120" height="90" border="0" cellspacing="0" cellpadding="0" style="width:120px;height:100px;border-radius: 15px 5px;background-color:#CFF;";>
            <tr>
              <td align="center" valign="middle"><span style="font-size:64px;color:#CC0099;text-decoration:underline;"><a href="http://websdoor.net/bioinfo/maxreport/" target="_blank" title="Visit MaxReport Homepage">MR</a></span></td>
            </tr>
          </table>
        </td>
        <td height="120" align="left" valign="middle"><span class="tit">Index page for proteomic result reporting</span><br>
          <span class="subtit">Generated by MaxReport on %s</span></td>
      </tr>
    </table></td>
  </tr>
  <tr>
    <td align="center" valign="top"><table width="95%%" border="0" cellpadding="0" cellspacing="5" class="dashedbox">
      <tr>
        <td><table width="100%%" border="0" align="center" cellpadding="0" cellspacing="0" class="solidbox">
          <tr>
            <td height="27" class="tibg">&nbsp;1. Input parameters</td>
            <td width="25%%" height="27" class="tibg"><table width="95%%" border="0" align="center" cellpadding="0" cellspacing="0">
              <tr>
                <td width="50%%" align="center" valign="middle" class="tibg" id="divth0" onclick="dh('divtc0');dh('divts0');dh('divth0');">- Hide Detail</td>
                <td width="50%%" align="center" valign="middle" class="tibg" id="divts0" style="display:none;" onclick="dh('divtc0');dh('divth0');dh('divts0');">+ Show Detail</td>
              </tr>
            </table></td>
          </tr>
          <tr id="divtc0">
            <td colspan="2" align="center" valign="top" class="ctbg"><table width="95%%" border="0" cellspacing="10" cellpadding="0">
    <tr>
      <td align="left" valign="top"><strong>Commands:</strong></td>
    </tr>
    <tr>
      <td align="left" valign="middle">%s</td>
    </tr>
  </table></td>
          </tr>
        </table></td>
      </tr>
      <tr>
        <td><table width="100%%" border="0" align="center" cellpadding="0" cellspacing="0" class="solidbox">
          <tr>
            <td height="27" class="tibg">&nbsp;2. Identification results</td>
            <td width="25%%" height="27" class="tibg"><table width="95%%" border="0" align="center" cellpadding="0" cellspacing="0">
              <tr>
                <td width="50%%" align="center" valign="middle" class="tibg" id="divth1" onclick="dh('divtc1');dh('divts1');dh('divth1');">- Hide Detail</td>
                <td width="50%%" align="center" valign="middle" class="tibg" id="divts1" style="display:none;" onclick="dh('divtc1');dh('divth1');dh('divts1');">+ Show Detail</td>
              </tr>
            </table></td>
          </tr>
          <tr id="divtc1">
            <td colspan="2" align="center" valign="top" class="ctbg"><table width="95%%" border="0" cellspacing="10" cellpadding="0">
              <tr>
                  <td align="center" valign="top"><strong>Table 2.1 Location of result files</strong></td>
                </tr>
                <tr>
                  <td><table width="100%%" border="0" align="center" cellpadding="0" cellspacing="0">
                    <tr>
                      <td width="30%%" align="left" valign="top" class="tltbhd">Description</td>
                      <td align="left" valign="top" class="tltbhd">File location</td>
                    </tr>
                    %s
                    <tr>
                      <td width="30%%" align="left" valign="top">Protein groups</td>
                      <td align="left" valign="top">%s</td>
                    </tr>
                    <tr>
                      <td align="left" valign="top">Peptides</td>
                      <td align="left" valign="top">%s</td>
                    </tr>
                    <tr>
                      <td align="left" valign="top">N-terminal modifications</td>
                      <td align="left" valign="top">%s</td>
                    </tr>
                    <tr>
                      <td align="left" valign="top">Modifications</td>
                      <td align="left" valign="top">%s</td>
                    </tr>
                    <tr>
                      <td align="left" valign="top" class="tltbbt">MS/MS spectra</td>
                      <td align="left" valign="top" class="tltbbt">%s</td>
                    </tr>
                  </table></td>
                </tr>
            </table></td>
          </tr>
        </table></td>
      </tr>
      <tr>
        <td align="center" valign="middle"><table width="100%%" border="0" align="center" cellpadding="0" cellspacing="0" class="solidbox">
          <tr>
            <td height="27" class="tibg">&nbsp;3. Summary and descriptive statistics</td>
            <td width="25%%" height="27" class="tibg"><table width="95%%" border="0" align="center" cellpadding="0" cellspacing="0">
              <tr>
                <td width="50%%" align="center" valign="middle" class="tibg" id="divth2" onclick="dh('divtc2');dh('divts2');dh('divth2');">- Hide Detail</td>
                <td width="50%%" align="center" valign="middle" class="tibg" id="divts2" style="display:none;" onclick="dh('divtc2');dh('divth2');dh('divts2');">+ Show Detail</td>
                </tr>
              </table></td>
            </tr>
          <tr id="divtc2">
            <td colspan="2" align="center" valign="middle" class="ctbg"><table width="95%%" border="0" align="center" cellpadding="0" cellspacing="10">
              <tr>
                <td align="left" valign="top"><strong>3.1 General proteome count</strong></td>
              </tr>
              <tr>
                <td align="center" valign="middle">%s</td>
              </tr>
              <tr>
                <td align="center" valign="middle"><strong>Figure 3.1 Summary of general protein identification</strong></td>
              </tr>
              </table>
              <table width="95%%" border="0" align="center" cellpadding="0" cellspacing="10">
                <tr>
                  <td align="left" valign="top"><strong>3.2 Protein-peptides count</strong></td>
                </tr>
                <tr>
                  <td align="center" valign="middle">
   %s</td>
                </tr>
                <tr>
                  <td align="center" valign="middle"><strong>Figure 3.2 Distribution of all peptide counts</strong></td>
                </tr>
                <tr>
                  <td align="center" valign="middle">%s</td>
                </tr>
                <tr>
                  <td align="center" valign="middle"><strong>Figure 3.3 Distribution of unique peptide counts</strong></td>
                </tr>
              </table>
              <table width="95%%" border="0" align="center" cellpadding="0" cellspacing="10">
                <tr>
                  <td align="left" valign="top"><strong>3.3 Post-translational modification (PTM) specific proteome count</strong></td>
                </tr>
                <tr>
                  <td align="center" valign="middle">%s</td>
                </tr>
                <tr>
                  <td align="center" valign="middle"><strong>Figure 3.4 Summary of PTM identification</strong></td>
                </tr>
              </table>
              <table width="95%%" border="0" align="center" cellpadding="0" cellspacing="10">
                <tr>
                  <td align="left" valign="top"><strong>3.4 Protein-sites count</strong></td>
                </tr>
%s
              </table>
              %s
              <table width="95%%" border="0" align="center" cellpadding="0" cellspacing="10">
                <tr>
                  <td align="center" valign="top"><strong>Table 3.1 Location of summary files</strong></td>
                </tr>
                <tr>
                  <td align="center" valign="middle"><table width="100%%" border="0" align="center" cellpadding="0" cellspacing="0">
                    <tr>
                      <td width="30%%" align="left" valign="top" class="tltbhd">Description</td>
                      <td align="left" valign="top" class="tltbhd">File location</td>
                      </tr>
%s
                    </table></td>
                </tr>
                </table></td>
            </tr>
          </table></td>
      </tr>
      <tr>
        <td><table width="100%%" border="0" align="center" cellpadding="0" cellspacing="0" class="solidbox">
          <tr>
            <td height="27" class="tibg">&nbsp;4. Quantification based on isobaric labeling</td>
            <td width="25%%" height="27" class="tibg"><table width="95%%" border="0" align="center" cellpadding="0" cellspacing="0">
              <tr>
                <td width="50%%" align="center" valign="middle" class="tibg" id="divth4" onclick="dh('divtc4');dh('divts4');dh('divth4');">- Hide Detail</td>
                <td width="50%%" align="center" valign="middle" class="tibg" id="divts4" style="display:none;" onclick="dh('divtc4');dh('divth4');dh('divts4');">+ Show Detail</td>
                </tr>
              </table></td>
            </tr>
          <tr id="divtc4">
            <td colspan="2" align="center" class="ctbg">%s</td>
            </tr>
          </table></td>
      </tr>
      <tr>
        <td><table width="100%%" border="0" align="center" cellpadding="0" cellspacing="0" class="solidbox">
          <tr>
            <td height="27" class="tibg">&nbsp;5. Miscellaneous</td>
            <td width="25%%" height="27" class="tibg"><table width="95%%" border="0" align="center" cellpadding="0" cellspacing="0">
              <tr>
                <td width="50%%" align="center" valign="middle" class="tibg" id="divth5" onclick="dh('divtc5');dh('divts5');dh('divth5');">- Hide Detail</td>
                <td width="50%%" align="center" valign="middle" class="tibg" id="divts5" style="display:none;" onclick="dh('divtc5');dh('divth5');dh('divts5');">+ Show Detail</td>
              </tr>
            </table></td>
          </tr>
          <tr id="divtc5">
            <td colspan="2" align="center" class="ctbg"><table width="95%%" border="0" align="center" cellpadding="0" cellspacing="10">
              <tr>
                <td align="center" valign="top"><strong>Table 5.1 Location of miscellaneous files</strong></td>
              </tr>
                <tr>
                  <td align="center" valign="middle"><table width="100%%" border="0" align="center" cellpadding="0" cellspacing="0">
                    <tr>
                      <td width="30%%" align="left" valign="top" class="tltbhd">Description</td>
                      <td align="left" valign="top" class="tltbhd">File location</td>
                    </tr>
                    <tr>
                      <td width="30%%" align="left" valign="top">Mapped protein sequences</td>
                      <td align="left" valign="top">%s</td>
                    </tr>
                    <tr>
                      <td align="left" valign="top" class="tltbbt">Mapped MS/MS scans</td>
                      <td align="left" valign="top" class="tltbbt">%s</td>
                    </tr>
                  </table></td>
                </tr>
            </table></td>
          </tr>
        </table></td>
      </tr>
      <tr>
        <td><table width="100%%" border="0" align="center" cellpadding="0" cellspacing="0" class="solidbox">
          <tr>
            <td height="27" class="tibg">&nbsp;6. Log file</td>
            <td width="25%%" height="27" class="tibg"><table width="95%%" border="0" align="center" cellpadding="0" cellspacing="0">
              <tr>
                <td width="50%%" align="center" valign="middle" class="tibg" id="divth6" onclick="dh('divtc6');dh('divts6');dh('divth6');">- Hide Detail</td>
                <td width="50%%" align="center" valign="middle" class="tibg" id="divts6" style="display:none;" onclick="dh('divtc6');dh('divth6');dh('divts6');">+ Show Detail</td>
              </tr>
            </table></td>
          </tr>
          <tr id="divtc6">
            <td colspan="2" align="center" class="ctbg"><table width="95%%" border="0" cellspacing="10" cellpadding="0">
              <tr>
                <td align="left" valign="top"><strong>File location:</strong></td>
              </tr>
              <tr>
                <td align="left" valign="middle"><a href="log.txt" target="_blank">log.txt</a></td>
              </tr>
            </table></td>
          </tr>
        </table></td>
      </tr>
    </table></td>
  </tr>
  <tr>
    <td align="center" valign="middle">      <span class="subtit">MaxReport: Enhanced proteomic result reporting tool for MaxQuant <br>
Version: %s<br>
      Designed and developped by Dr. Tao Zhou<br>
Contact: <a href="mailto:zhoutao@njmu.edu.cn">zhoutao@njmu.edu.cn</a></span></td>
  </tr>
</table>
</body>
</html>
'''%tuple(htmlis)
    topf=open(os.path.join(sdir,'maxreport','index.html'),'w')
    topf.writelines(htmltmp)
    topf.close()
    loginfo('Index file created: %s' %('index.html'))
    webbrowser.open(os.path.join(sdir,'maxreport','index.html'))
    loginfo('Index file opened!')
    loginfo('All finished!')
    lgopf.close()













