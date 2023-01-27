# bash command: python /location/of/code/histograms.py -a /path/rootfile1.root -b /path/rootfile2.root -c /path/rootfile3.root -d /path/rootfile4.root -e /path/rootfile5.root -f/path/rootfile6.root -m Mq -p Nf -q Nc1 -r Nc2 -s Nc3 -t Nc4
# Takes rootfiles as input and prints different histograms on separate canvas.
# Per plot, four different histograms for foud different colors and one given flavour are plotted

## Version 1

# - TO DO LIST:
# - Dynamic ranging for x,y axis
# - Will be adding histograms for cumulative distributions and invariant masses later (this week)
# - Will be adding back in multiple root files later (future)
# - Later will be adding histograms of METs and tracks back into the fold (far future)


# bash command: python /location/of/code/histograms.py -a /path/rootfile1.root -m Mq -p Nf -q Nc1

from calendar import c
import os, sys
from re import M
import string
from tkinter import E
import numpy as np
from numpy import *
from pylab import *
from scipy.interpolate import interp1d
import ROOT
from ROOT import gROOT
from array import array
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c' , help = "Nf")
parser.add_argument('-d' , help = "pion_decay")
args = parser.parse_args()


Nf=args.c
pion_decay=args.d

eventcounter=0
store={}
Nc_array=[]

#Nc_array=array(3,7,11,15,19,24)
#####
leadingpt={}
subleadingpt={}
jetmjj={}
jetdeta={}
jetpiPT={}
jetpinum={}
jetmet={}

for entry in os.listdir('/media/data/lockyer/Nf_variation_analysis/'):
        for i in range(3,25):
            if entry == f"Nc{i}_Nf{Nf}_pp_{pion_decay}pi_decay_lam_5.root":
                print(entry)
                string = f'/media/data/lockyer/Nf_variation_analysis/{entry}'
                store[i] = ROOT.TFile(string)
                Nc_array.append(i)
                eventcounter=eventcounter+1
print(Nc_array)
print(store)
print(eventcounter)

for i in Nc_array:         
    leadingpt[i] = (store[i].Get('jet1')).GetMean()
    subleadingpt[i] = (store[i].Get('jet2')).GetMean()
    jetmjj[i] = (store[i].Get('jet_mJJ')).GetMean()
    jetdeta[i] = (store[i].Get('delta_eta')).GetMean()
    jetpiPT[i] = (store[i].Get('pion_pt')).GetMean()
    jetpinum[i] = (store[i].Get('pion_number')).GetMean()
    jetmet[i] = (store[i].Get('jet_met')).GetMean()
        
colourguide_lambda={1: 8, 2: 9, 3: 46, 4: 8, 5: 9,6: 46, 7: 8,8: 9,9: 46}
colourguide_pinum={1: 8, 2: 8, 3: 8, 4: 9, 5: 9,6: 9, 7: 46,8: 46,9: 46}

########################################################################
c1 = ROOT.TCanvas( 'c1', 'Leading pT against Nc', 200, 10, 700, 500 )
c1.SetFillColor(0)
Nc, leadingptarray = array('d'), array('d')
for i in Nc_array:
    Nc.append(i)
    leadingptarray.append(leadingpt[i])
gr = ROOT.TGraph(eventcounter, Nc, leadingptarray)
gr.SetLineColor(1)
gr.SetLineWidth(4)
gr.SetMarkerColor(1)
gr.SetMarkerStyle(8)
gr.SetTitle( "Mean leading jet p_{T} for N_{f}=%s"%(Nf))
gr.GetXaxis().SetTitle( 'N_{c}' )
gr.GetYaxis().SetTitle("Mean leading p_{T}[GeV]")
gr.GetYaxis().SetRangeUser(140,180)
gr.Draw( 'AP' )

t1 = ROOT.TLatex()
t1.SetNDC(ROOT.kTRUE)
t1.SetTextSize(0.03)
t1.DrawLatex(0.18 ,0.86 ,"N_{f} = %s , %s pion decay , \Lambda = 5 GeV"%(Nf,pion_decay))
t1.SetTextSize(0.03)

c1.Update()
c1.GetFrame().SetFillColor(0)
c1.GetFrame().SetBorderSize(5)
c1.Modified()
c1.Update()
c1.Print("MeanleadpT_%s_%s.pdf"%(Nf,pion_decay))
########################################################################
c1 = ROOT.TCanvas( 'c2', 'subleading against Nc', 200, 10, 700, 500 )
c1.SetFillColor(0)
Nc, subleadingptarray = array('d'), array('d')
for i in Nc_array:
    Nc.append(i)
    subleadingptarray.append(subleadingpt[i])
gr = ROOT.TGraph(eventcounter, Nc, subleadingptarray)
gr.SetLineColor(1)
gr.SetLineWidth(4)
gr.SetMarkerColor(1)
gr.SetMarkerStyle(8)
gr.SetTitle( "Mean subleading jet p_{T} for N_{f}=%s"%(Nf))
gr.GetXaxis().SetTitle( 'N_{c}' )
gr.GetYaxis().SetTitle("Mean subleading p_{T}[GeV]")
gr.GetYaxis().SetRangeUser(75,100)
gr.Draw( 'AP' )

t1 = ROOT.TLatex()
t1.SetNDC(ROOT.kTRUE)
t1.SetTextSize(0.03)
t1.DrawLatex(0.18 ,0.86 ,"N_{f} = %s , %s pion decay , \Lambda = 5 GeV"%(Nf,pion_decay))
t1.SetTextSize(0.03)

c1.Update()
c1.GetFrame().SetFillColor(0)
c1.GetFrame().SetBorderSize(5)
c1.Modified()
c1.Update()
c1.Print("MeansubleadpT_%s_%s.pdf"%(Nf,pion_decay))
########################################################################
c1 = ROOT.TCanvas( 'c3', 'minv against Nc', 200, 10, 700, 500 )
c1.SetFillColor(0)
Nc, mjjarray = array('d'), array('d')
for i in Nc_array:
    Nc.append(i)
    mjjarray.append(jetmjj[i])
gr = ROOT.TGraph(eventcounter, Nc, mjjarray)
gr.SetLineColor(1)
gr.SetLineWidth(4)
gr.SetMarkerColor(1)
gr.SetMarkerStyle(8)
gr.SetTitle( "Mean invariant mass for N_{f}=%s"%(Nf))
gr.GetXaxis().SetTitle( 'N_{c}' )
gr.GetYaxis().SetTitle("Mean invariant mass [GeV]")
gr.GetYaxis().SetRangeUser(260,310)
gr.Draw( 'AP' )

t1 = ROOT.TLatex()
t1.SetNDC(ROOT.kTRUE)
t1.SetTextSize(0.03)
t1.DrawLatex(0.18 ,0.86 ,"N_{f} = %s , %s pion decay , \Lambda = 5 GeV"%(Nf,pion_decay))
t1.SetTextSize(0.03)

c1.Update()
c1.GetFrame().SetFillColor(0)
c1.GetFrame().SetBorderSize(5)
c1.Modified()
c1.Update()
c1.Print("Invariant_mass_%s_%s.pdf"%(Nf,pion_decay))
########################################################################
c1 = ROOT.TCanvas( 'c4', 'deleta against Nc', 200, 10, 700, 500 )
c1.SetFillColor(0)
Nc, detaarray = array('d'), array('d')
for i in Nc_array:
    Nc.append(i)
    detaarray.append(jetdeta[i])
gr = ROOT.TGraph(eventcounter, Nc, detaarray)
gr.SetLineColor(1)
gr.SetLineWidth(4)
gr.SetMarkerColor(1)
gr.SetMarkerStyle(8)
gr.SetTitle( "Mean \Delta\eta for N_{f}=%s"%(Nf))
gr.GetXaxis().SetTitle( 'N_{c}' )
gr.GetYaxis().SetTitle("\Delta\eta")
gr.GetYaxis().SetRangeUser(1,2)
gr.Draw( 'AP' )

t1 = ROOT.TLatex()
t1.SetNDC(ROOT.kTRUE)
t1.SetTextSize(0.03)
t1.DrawLatex(0.18 ,0.86 ,"N_{f} = %s , %s pion decay , \Lambda = 5 GeV"%(Nf,pion_decay))
t1.SetTextSize(0.03)

c1.Update()
c1.GetFrame().SetFillColor(0)
c1.GetFrame().SetBorderSize(5)
c1.Modified()
c1.Update()
c1.Print("delta_eta_%s_%s.pdf"%(Nf,pion_decay))
########################################################################
c1 = ROOT.TCanvas( 'c5', 'met', 200, 10, 700, 500 )
c1.SetFillColor(0)
Nc, metarray = array('d'), array('d')
for i in Nc_array:
    Nc.append(i)
    metarray.append(jetmet[i])
gr = ROOT.TGraph(eventcounter, Nc, metarray)
gr.SetLineColor(1)
gr.SetLineWidth(4)
gr.SetMarkerColor(1)
gr.SetMarkerStyle(8)
gr.SetTitle( "Mean #slash{E}_{T} for N_{f}=%s"%(Nf))
gr.GetXaxis().SetTitle( 'N_{c}' )
gr.GetYaxis().SetTitle("#slash{E}_{T} [GeV]")
gr.GetYaxis().SetRangeUser(80,120)
gr.Draw( 'AP' )

t1 = ROOT.TLatex()
t1.SetNDC(ROOT.kTRUE)
t1.SetTextSize(0.03)
t1.DrawLatex(0.18 ,0.86 ,"N_{f} = %s , %s pion decay , \Lambda = 5 GeV"%(Nf,pion_decay))
t1.SetTextSize(0.03)

c1.Update()
c1.GetFrame().SetFillColor(0)
c1.GetFrame().SetBorderSize(5)
c1.Modified()
c1.Update()
c1.Print("jet_met_%s_%s.pdf"%(Nf,pion_decay))
########################################################################
c1 = ROOT.TCanvas( 'c6', 'pionpt', 200, 10, 700, 500 )
c1.SetFillColor(0)
Nc, pionptarray = array('d'), array('d')
for i in Nc_array:
    Nc.append(i)
    pionptarray.append(jetpiPT[i])
gr = ROOT.TGraph(eventcounter, Nc, pionptarray)
gr.SetLineColor(1)
gr.SetLineWidth(4)
gr.SetMarkerColor(1)
gr.SetMarkerStyle(8)
gr.SetTitle( "Mean pion p_{T} for N_{f}=%s"%(Nf))
gr.GetXaxis().SetTitle( 'N_{c}' )
gr.GetYaxis().SetTitle("Mean pion p_{T} [GeV]")
gr.GetYaxis().SetRangeUser(20,40)
gr.Draw( 'AP' )

t1 = ROOT.TLatex()
t1.SetNDC(ROOT.kTRUE)
t1.SetTextSize(0.03)
t1.DrawLatex(0.18 ,0.86 ,"N_{f} = %s , %s pion decay , \Lambda = 5 GeV"%(Nf,pion_decay))
t1.SetTextSize(0.03)

c1.Update()
c1.GetFrame().SetFillColor(0)
c1.GetFrame().SetBorderSize(5)
c1.Modified()
c1.Update()
c1.Print("pion_pt_%s_%s.pdf"%(Nf,pion_decay))
########################################################################
c1 = ROOT.TCanvas( 'c7', 'pionnumber', 200, 10, 700, 500 )
c1.SetFillColor(0)
Nc, pionnumarray = array('d'), array('d')
for i in Nc_array:
    Nc.append(i)
    pionnumarray.append(jetpinum[i])
gr = ROOT.TGraph(eventcounter, Nc, pionnumarray)
gr.SetLineColor(1)
gr.SetLineWidth(4)
gr.SetMarkerColor(1)
gr.SetMarkerStyle(8)
gr.SetTitle( "Mean pion number for N_{f}=%s"%(Nf))
gr.GetXaxis().SetTitle( 'N_{c}' )
gr.GetYaxis().SetTitle("Mean pion number")
gr.GetYaxis().SetRangeUser(20,40)
gr.Draw( 'AP' )

t1 = ROOT.TLatex()
t1.SetNDC(ROOT.kTRUE)
t1.SetTextSize(0.03)
t1.DrawLatex(0.18 ,0.86 ,"N_{f} = %s , %s pion decay , \Lambda = 5 GeV"%(Nf,pion_decay))
t1.SetTextSize(0.03)

c1.Update()
c1.GetFrame().SetFillColor(0)
c1.GetFrame().SetBorderSize(5)
c1.Modified()
c1.Update()
c1.Print("pion_number_%s_%s.pdf"%(Nf,pion_decay))





