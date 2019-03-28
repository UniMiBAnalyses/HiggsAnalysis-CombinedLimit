#!/usr/bin/env python

import ROOT
import CombineHarvester.CombineTools.plotting as plot


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(0)

plot.ModTDRStyle(width=700, height=1000)

def MakeGraphs(func, var, points):
    fd1 = func.derivative(var, 1)
    fd2 = func.derivative(var, 2)
    xmin = var.getMin()
    xmax = var.getMax()
    gr = ROOT.TGraph(points)
    grd1 = ROOT.TGraph(points)
    grd2 = ROOT.TGraph(points)
    gr.SetName(var.GetName())
    grd1.SetName(var.GetName()+"_d1")
    grd2.SetName(var.GetName()+"_d2")
    w = (xmax - xmin) / float(points)
    for i in xrange(points):
        x = xmin + (float(i) + 0.5) * w
        var.setVal(x)
        gr.SetPoint(i, x, func.getVal())
        grd1.SetPoint(i, x, fd1.getVal())
        grd2.SetPoint(i, x, fd2.getVal())
    return (gr, grd1, grd2)


def MakePlots(name, graphs):
    canv = ROOT.TCanvas(name, name)
    pads = plot.MultiRatioSplit([0.4, 0.3], [0.005, 0.005], [0.005, 0.005])
    pads[0].cd()
    plot.Set(graphs[0], MarkerSize=0.5)
    graphs[0].Draw('AL')
    axis1 = plot.GetAxisHist(pads[0])
    axis1.GetYaxis().SetTitle('fn')
    pads[1].cd()
    plot.Set(graphs[1], MarkerSize=0.5)
    graphs[1].Draw('AL')
    graphs[1].Print()
    axis2 = plot.GetAxisHist(pads[1])
    axis2.GetYaxis().SetTitle('fn\'')
    pads[2].cd()
    plot.Set(graphs[2], MarkerSize=0.5)
    graphs[2].Draw('AL')
    graphs[2].Print()
    axis3 = plot.GetAxisHist(pads[2])
    axis3.GetYaxis().SetTitle('fn\'\'')
    plot.Set(axis3.GetXaxis(),Title=x.GetName(),
        TitleSize=axis3.GetXaxis().GetTitleSize()*0.5,
        TitleOffset=axis3.GetXaxis().GetTitleOffset()*2,
        )
    canv.Print('.pdf')
    canv.Print('.png')



x = ROOT.RooRealVar('x', '', 0, -1, 1)
mu = ROOT.RooRealVar('mu', '', 10.)
N = ROOT.RooRealVar('N', '', 10.5)

lo = -0.1
hi = 1.0
r1 = (hi - lo) / 2.
r2 = (hi + lo) * (15. / 8.)


iCurrent = ROOT.RooCMSInterpVar('current', '', x, lo, hi, 1, 0.93, 0.65)

nllCurrent = ROOT.RooFormulaVar('nll', '-2. * (@0*log(@1+@2) - (@1+@2))', ROOT.RooArgList(N, mu, iCurrent))

gCurrent = MakeGraphs(iCurrent, x, 100)
gNLLCurrent = MakeGraphs(nllCurrent, x, 100)
MakePlots('testInterVar', gCurrent)
MakePlots('testNLLVar', gNLLCurrent)

x.setVal(-0.4)
minim = ROOT.RooMinimizer(nllCurrent)
minim.setVerbose(True)
minim.setStrategy(1)
minim.minimize('Minuit2', 'migrad')
# minim.setPrintLevel(-1)

print (r1, r2)

