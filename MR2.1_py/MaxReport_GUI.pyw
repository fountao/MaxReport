# -*- coding: UTF-8 -*-
'''
MaxReport_GUI for MaxReport_CMD
'''
import wx
import wx.gizmos as gizmos
import sys
import subprocess
import os

def resource_path(relative_path):
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

class MainFrame(wx.MiniFrame):
    def __init__(
        self,parent,title,pos=wx.DefaultPosition,size=wx.DefaultSize,
        style=wx.DEFAULT_FRAME_STYLE^wx.RESIZE_BORDER|wx.TINY_CAPTION_HORIZ|wx.STAY_ON_TOP
        #style=wx.NO_BORDER| wx.STAY_ON_TOP | wx.FRAME_NO_TASKBAR | wx.FRAME_SHAPED
        ):
        self.title=title
        wx.MiniFrame.__init__(self, parent, -1, title, pos, size, style)
        self.title=title
        #set frame style
        self.SetBackgroundColour('#66CCCC')
        self.SetForegroundColour('#000000')
        self.SetFont(wx.FFont(12,wx.FONTFAMILY_DEFAULT))
        #create main box
        self.mainsizer = wx.BoxSizer(wx.VERTICAL)
        #load logo image
        self.logo=wx.StaticBitmap(self,-1,wx.Image(resource_path('mr_logo.png'), wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        self.logo.SetToolTipString('MaxReport GUI!')
        #required arguments box
        rqdbox=wx.StaticBox(self, -1, 'Required arguments')
        self.rqdsizer = wx.StaticBoxSizer(rqdbox, wx.VERTICAL)
        #create sdir sizer
        self.sdirsizer=wx.BoxSizer(wx.HORIZONTAL)
        sdirtit=wx.StaticText(self,wx.ID_ANY,'Source folder:')
        self.sdir=wx.TextCtrl(self)
        self.sdir.SetToolTipString("Source folder for a MaxQuant project.")
        sdirbut=wx.Button(self,label='Select')
        sdirbut.Bind(wx.EVT_BUTTON,self.opensdir)
        self.sdirsizer.Add(sdirtit,0,wx.LEFT,0)
        self.sdirsizer.Add(self.sdir,1,wx.LEFT|wx.EXPAND,0)
        self.sdirsizer.Add(sdirbut,0,wx.LEFT,0)
        #create rcfg sizer
        self.rcfgsizer=wx.BoxSizer(wx.HORIZONTAL)
        rcfgtit=wx.StaticText(self,wx.ID_ANY,'Configuration file:')
        self.rcfg=wx.TextCtrl(self)
        self.rcfg.SetToolTipString("Configuration file for result reporting.")
        rcfgbut=wx.Button(self,label='Select')
        rcfgbut.Bind(wx.EVT_BUTTON,self.openrcfg)
        self.rcfgsizer.Add(rcfgtit,0,wx.LEFT,0)
        self.rcfgsizer.Add(self.rcfg,1,wx.LEFT|wx.EXPAND,0)
        self.rcfgsizer.Add(rcfgbut,0,wx.LEFT,0)
        #create frules sizer
        self.frssizer=wx.BoxSizer(wx.HORIZONTAL)
        frstit=wx.StaticText(self,wx.ID_ANY,'Rule files:')
        self.frs=gizmos.EditableListBox(self, -1, "List of files", size=(240,120))
        frsbut=wx.Button(self,label='Select')
        frsbut.SetToolTipString("Rule files for parsing fasta sequences.\nkeep file order the same with sequence files.")
        frsbut.Bind(wx.EVT_BUTTON,self.openfrs)
        self.frssizer.Add(frstit,0,wx.LEFT,0)
        self.frssizer.Add(self.frs,1,wx.LEFT|wx.EXPAND,0)
        self.frssizer.Add(frsbut,0,wx.LEFT,0)
        #create fseqs sizer
        self.fsssizer=wx.BoxSizer(wx.HORIZONTAL)
        fsstit=wx.StaticText(self,wx.ID_ANY,'Sequence files:')
        self.fss=gizmos.EditableListBox(self, -1, "List of files", size=(240,120))
        fssbut=wx.Button(self,label='Select')
        fssbut.SetToolTipString("Protein sequence files.\nkeep file order the same with rule files.")
        fssbut.Bind(wx.EVT_BUTTON,self.openfss)
        self.fsssizer.Add(fsstit,0,wx.LEFT,0)
        self.fsssizer.Add(self.fss,1,wx.LEFT|wx.EXPAND,0)
        self.fsssizer.Add(fssbut,0,wx.LEFT,0)
        #assemble required box
        self.rqdsizer.Add(self.sdirsizer,0,wx.ALL|wx.EXPAND,0)
        self.rqdsizer.Add(self.rcfgsizer,0,wx.ALL|wx.EXPAND,0)
        self.rqdsizer.Add(self.frssizer,0,wx.ALL|wx.EXPAND,0)
        self.rqdsizer.Add(self.fsssizer,0,wx.ALL|wx.EXPAND,0)
        #create optional box
        optbox=wx.StaticBox(self, -1, 'Optional arguments')
        self.optsizer = wx.StaticBoxSizer(optbox, wx.VERTICAL)
        #create exp sizer
        self.expsizer=wx.BoxSizer(wx.HORIZONTAL)
        exptit=wx.StaticText(self,wx.ID_ANY,'Experimental design:')
        self.exp=wx.TextCtrl(self)
        self.exp.SetToolTipString("File for experimental design.")
        expbut=wx.Button(self,label='Select')
        expbut.Bind(wx.EVT_BUTTON,self.openexp)
        self.expsizer.Add(exptit,0,wx.LEFT,0)
        self.expsizer.Add(self.exp,1,wx.LEFT|wx.EXPAND,0)
        self.expsizer.Add(expbut,0,wx.LEFT,0)
        #create cor sizer
        self.corsizer=wx.BoxSizer(wx.HORIZONTAL)
        cortit=wx.StaticText(self,wx.ID_ANY,'Correction matrix:')
        self.cor=wx.TextCtrl(self)
        self.cor.SetToolTipString("Correction matrix file for isobaric quantification.")
        corbut=wx.Button(self,label='Select')
        corbut.Bind(wx.EVT_BUTTON,self.opencor)
        self.corsizer.Add(cortit,0,wx.LEFT,0)
        self.corsizer.Add(self.cor,1,wx.LEFT|wx.EXPAND,0)
        self.corsizer.Add(corbut,0,wx.LEFT,0)
        #create mini sizer
        self.minisizer=wx.BoxSizer(wx.HORIZONTAL)
        minitit=wx.StaticText(self,wx.ID_ANY,'Threshold value:')
        self.mini=wx.TextCtrl(self,size=(100, -1))
        self.mini.SetToolTipString("Minimum threshold value for isobaric quantification.")
        self.minisizer.Add(minitit,0,wx.LEFT,0)
        self.minisizer.Add(self.mini,0,wx.LEFT|wx.EXPAND,0)
        #create pep checkbox
        self.pep=wx.CheckBox(self, -1, "Use all peptides for quantification.")
        #create xls checkbox
        self.xls=wx.CheckBox(self, -1, "Write results to xlsx files.")
        self.xls.SetValue(True)
        #assemble optional box
        self.optsizer.Add(self.expsizer,0,wx.ALL|wx.EXPAND,0)
        self.optsizer.Add(self.corsizer,0,wx.ALL|wx.EXPAND,0)
        self.optsizer.Add(self.minisizer,0,wx.ALL|wx.EXPAND,0)
        self.optsizer.Add(self.pep,0,wx.ALL|wx.EXPAND,0)
        self.optsizer.Add(self.xls,0,wx.ALL|wx.EXPAND,0)
        #create control sizer
        self.ctrlsizer=wx.BoxSizer(wx.HORIZONTAL)
        runbut=wx.Button(self,label='Run')
        runbut.SetToolTipString('Run MaxReport.')
        runbut.Bind(wx.EVT_BUTTON,self.runmr)
        rstbut=wx.Button(self,label='Reset')
        rstbut.SetToolTipString('Clear all the settings.')
        rstbut.Bind(wx.EVT_BUTTON,self.reset)
        closebut=wx.Button(self,label='Exit')
        self.ctrlsizer.Add(runbut,proportion=0,flag=wx.ALL,border=0)
        self.ctrlsizer.Add(rstbut,proportion=0,flag=wx.ALL,border=0)
        self.ctrlsizer.Add(closebut,proportion=0,flag=wx.ALL,border=0)
        closebut.Bind(wx.EVT_BUTTON,self.OnClose)
        #assemble main sizer
        self.mainsizer.Add(self.logo,0,wx.ALL|wx.CENTER,0)
        self.mainsizer.Add(self.rqdsizer,0,wx.ALL|wx.EXPAND,0)
        self.mainsizer.Add(self.optsizer,0,wx.ALL|wx.EXPAND,0)
        self.mainsizer.Add(wx.StaticLine(self,-1),0,wx.ALL|wx.EXPAND,0)
        self.mainsizer.Add(self.ctrlsizer,0,wx.ALL|wx.ALIGN_CENTER,0)
        self.SetSizer(self.mainsizer)
        #set the correct size of window
        self.Fit()
        #mainsizer.Fit(self)
        self.Center()
        #creat fade in and fade out feature
        self.amount = 20
        self.tincrs = 20
        self.stimer = wx.Timer(self,id=1)
        self.SetTransparent(self.amount)
        self.stimer.Start(50)
        self.Bind(wx.EVT_TIMER, self.fadein,id=1)
        self.ctimer = wx.Timer(self,id=2)
        self.Bind(wx.EVT_TIMER, self.fadeout,id=2)
        #rewrite close event
        self.Bind(wx.EVT_CLOSE,self.OnClose)
        self.Show()
        #check MaxReport CMD program:
        self.crtdir=os.path.split(os.path.realpath(__file__))[0]
        if not (os.path.isfile(os.path.join(self.crtdir,'MaxReport_CMD.exe')) or os.path.isfile(os.path.join(self.crtdir,'MaxReport_CMD.py'))):
            self.notedlg('MaxReport CMD is missing!')
            wx.Exit()

    #internal class functions
    def opensdir(self,event):
        dlg=wx.DirDialog(
            self,message="Choose result directory of a MaxReport project:",
            style=wx.DD_DEFAULT_STYLE|wx.DD_DIR_MUST_EXIST|wx.DD_CHANGE_DIR
            )
        if dlg.ShowModal()==wx.ID_OK:
            apath=dlg.GetPath()
            self.sdir.SetValue(apath)
        dlg.Destroy()

    def openrcfg(self,event):
        dlg=wx.FileDialog(
            self,message="Choose configuration file for result reporting:",
            defaultFile="",
            wildcard="CFG file (*.cfg)|*.cfg|Text file (*.txt)|*.txt|All files (*.*)|*.*",
            style=wx.OPEN
            )
        if dlg.ShowModal()==wx.ID_OK:
            apath=dlg.GetPath()
            self.rcfg.SetValue(apath)
        dlg.Destroy()

    def openexp(self,event):
        dlg=wx.FileDialog(
            self,message="Choose file for experimental design:",
            defaultFile="",
            wildcard="Text file (*.txt)|*.txt|All files (*.*)|*.*",
            style=wx.OPEN
            )
        if dlg.ShowModal()==wx.ID_OK:
            apath=dlg.GetPath()
            self.exp.SetValue(apath)
        dlg.Destroy()

    def opencor(self,event):
        dlg=wx.FileDialog(
            self,message="Choose correction matrix file for isobaric quantification:",
            defaultFile="",
            wildcard="CFG file (*.cfg)|*.cfg|Text file (*.txt)|*.txt|All files (*.*)|*.*",
            style=wx.OPEN
            )
        if dlg.ShowModal()==wx.ID_OK:
            apath=dlg.GetPath()
            self.cor.SetValue(apath)
        dlg.Destroy()

    def openfrs(self,event):
        dlg=wx.FileDialog(
            self,message="Choose rule files for parsing fasta sequences:",
            defaultFile="",
            wildcard="CFG file (*.cfg)|*.cfg|Text file (*.txt)|*.txt|All files (*.*)|*.*",
            style=wx.OPEN|wx.FD_MULTIPLE
            )
        if dlg.ShowModal()==wx.ID_OK:
            mpaths=dlg.GetPaths()
            self.frs.SetStrings(self.frs.GetStrings()+mpaths)
        dlg.Destroy()

    def openfss(self,event):
        dlg=wx.FileDialog(
            self,message="Choose protein sequence files:",
            defaultFile="",
            wildcard="FASTA file (*.fasta)|*.fasta|FA file (*.fa)|*.fa|All files (*.*)|*.*",
            style=wx.OPEN|wx.FD_MULTIPLE
            )
        if dlg.ShowModal()==wx.ID_OK:
            mpaths=dlg.GetPaths()
            self.fss.SetStrings(self.fss.GetStrings()+mpaths)
        dlg.Destroy()

    def OnClose(self, event):
        if not self.ctimer.IsRunning() and (not self.stimer.IsRunning()):
            self.ctimer.Start(50)

    def reset(self, event):
        self.sdir.SetValue('')
        self.rcfg.SetValue('')
        self.exp.SetValue('')
        self.cor.SetValue('')
        self.frs.SetStrings([])
        self.fss.SetStrings([])
        self.mini.SetValue('')
        self.pep.SetValue(False)
        self.xls.SetValue(True)

    def runmr(self, event):
        if os.path.isfile(os.path.join(self.crtdir,'MaxReport_CMD.py')):
            cmdlis=['python',os.path.join(self.crtdir,'MaxReport_CMD.py')]
        else:
             cmdlis=[os.path.join(self.crtdir,'MaxReport_CMD.exe')]
        sdir=self.sdir.GetValue()
        rcfg=self.rcfg.GetValue()
        frlis=self.frs.GetStrings()
        fslis=self.fss.GetStrings()
        if not (sdir and rcfg and frlis and fslis):
            self.notedlg('Please ensure all required fields are filled!!')
            return
        cmdlis.append(sdir)
        cmdlis.append(rcfg)
        cmdlis.append('-r')
        cmdlis+=frlis
        cmdlis.append('-s')
        cmdlis+=fslis
        expt=self.exp.GetValue()
        if expt:
            cmdlis.append('-e')
            cmdlis.append(expt)
        corr=self.cor.GetValue()
        if corr:
            cmdlis.append('-c')
            cmdlis.append(corr)
        mini=self.mini.GetValue()
        if mini:
            cmdlis.append('-m')
            cmdlis.append(mini)
        peps=self.pep.GetValue()
        if peps:
            cmdlis.append('-a')
        xlsx=self.xls.GetValue()
        if xlsx:
            cmdlis.append('-x')
        os.chdir(self.crtdir)
        self.Hide()
        self.s2p=subprocess.Popen(cmdlis,stdout=None,stderr=None,shell=False)
        #self.s2p=subprocess.Popen(cmdlis,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=False)
        self.notedlg('MaxReport started successfully!\nPlease wait for the response ...','Note',wx.OK|wx.ICON_INFORMATION)
        self.s2p.communicate()
        self.notedlg('MaxReport finished!\nPlease check the index file or log file!','Note',wx.OK|wx.ICON_INFORMATION)
        #self.Show()
        wx.Exit()

    def fadein(self, evt):
        self.amount += self.tincrs
        if self.amount >= 240:
            self.amount = 240
            self.stimer.Stop()
        self.SetTransparent(self.amount)

    def fadeout(self, evt):
        self.amount -= self.tincrs
        if self.amount <= 0:
            self.amount=0
            self.Destroy()
        self.SetTransparent(self.amount)

    def autosize(self,*vsizes):
        tmpwidth=0
        tmpheight=0
        for es in vsizes:
            if es[0]>tmpwidth:
                tmpwidth=es[0]
            tmpheight+=es[1]
        return (tmpwidth,tmpheight)

    def notedlg(self,ntmes,ntit='Warning',ntsty=wx.OK|wx.ICON_WARNING):
        dlg = wx.MessageDialog(self, ntmes, ntit,ntsty)
        dlg.ShowModal()
        dlg.Destroy()
        

# Run the program
if __name__ == '__main__':
    app=wx.App(False)
    frame=MainFrame(None,title='GUI for MaxReport 2.0')
    app.MainLoop()
