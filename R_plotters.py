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
from lib2to3.pgen2.pgen import NFAState
import os, sys
from tkinter import E
from tokenize import Pointfloat
import numpy as np
from numpy import *
from pylab import *
from scipy.interpolate import interp1d
import ROOT
from array import array
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-a' , help = "rootfile1")
parser.add_argument('-b' , help = "rootfile2")
#parser.add_argument('-c' , help = "rootfile3")
#parser.add_argument('-d' , help = "rootfile4")
#parser.add_argument('-e' , help = "variable choice")
parser.add_argument('-f' , help = "const1")
parser.add_argument('-g' , help = "const2")
parser.add_argument('-i' , help = "const3")
parser.add_argument('-j' , help = "const4")
#parser.add_argument('-y' , help = "variable2")
#parser.add_argument('-z' , help = "variable3")

 


######VARIABLE CHOICE GUIDE######
### -e 1 is Lambda Variation
### -f pion decay
### -g Nc
### -h Nf
######
### -e 2  is Pion decay Variation
### -f lambda
### -g Nc
### -h Nf
######
### -e 3  is Nc Variation
### -f lambda
### -g pion decay 
### -h Nf
######
### -e 4  is Nf Variation
### -f lambda
### -g pion decay [1 for single pion, all for all pions]
### -h Nc


######PARSE ARGUMENTS######
args = parser.parse_args()
rootfile1 = args.a
rootfile2 = args.b
#rootfile3 = args.c
#rootfile4 = args.d
#variable_choice=args.e
const1=args.f
const2=args.g
const3=args.i
const4=args.j
#variable1=args.x
#variable2=args.y
#variable3=args.z
print("Input files :")
print("rootfile1: ",rootfile1)
print("rootfile2: ",rootfile2)
#print("rootfile3: ",rootfile3)
#print("rootfile4: ",rootfile4)



legend_dict={
    "1": "\Lambda", 
    "2": "pion decay", 
    "3": "N_{c}",
    "4": "N_{f}"}

if const2 == "all":
   const2 = "All pions"
   print("All pions")
else:
   const2 = "1 pion"

#(0.18 ,0.86 ,"\Lambda = %s, %s decay, N_{f}=%s, N_{c}=%s"%(const1,const2,const3,const4))



#Dictionary full of jet histograms
jet1hist={}
jet2hist={}
#jet3hist={}
#jet4hist={}
#Dictionary full of canvasses
canvas_dictionary={}
legend={}
meancount={}
#Counts the number of jets
jet1counter=0
jet2counter=0
#jet3counter=0
#jet4counter=0
#Dummy number for event loop
dummy=0

#hist_dict_p={}
#hist_dict_n={}
x_axis_dict={
    "1": "\Delta\eta", 
    "2": "\Delta\phi", 
    "3": "\Delta\eta",
    "4": "\Delta\phi"}
title_dict={
    "1": "Asymmetry distribution of \Delta\eta for leading and subleading jets", 
    "2": "Asymmetry distribution of \Delta\phi for leading and subleading jets", 
    "3": "Asymmetry distribution of \Delta\eta of two SM quarks",
    "4": "Asymmetry distribution of \Delta\phi of two SM quarks"}

######IMPORTING FILES######
f1 = ROOT.TFile(rootfile1)
hist1mjj = f1.Get('jet_mJJ')
hist1delEta = f1.Get("delta_eta")
hist1piPT = f1.Get("pion_pt")
hist1pinum=f1.Get("pion_number")
hist1met=f1.Get("jet_met")
hist1dPhi=f1.Get("dPhi")
hist1_dE_q=f1.Get("delta_eta_quark")
hist1_dP_q=f1.Get("dPhi_quark")
hist1_pT_q=f1.Get("quark_pt")
hist1_num_q=f1.Get("number_of_quarks")
while dummy == 0:
   jet1hist[jet1counter+1] = f1.Get(f'jet{jet1counter+1}')
   if not jet1hist[jet1counter+1]:
      dummy = 1
      print("File 1 has ",jet1counter," jets")
   else:
      meancount[jet1counter+1]=jet1hist[jet1counter+1].GetMean()
      jet1counter = jet1counter + 1  
dummy=0
##############
f2 = ROOT.TFile(rootfile2)
hist2mjj = f2.Get('jet_mJJ')
hist2delEta = f2.Get("delta_eta")
hist2piPT = f2.Get("pion_pt")
hist2pinum=f2.Get("pion_number")
hist2met=f2.Get("jet_met")
hist2dPhi=f2.Get("dPhi")
hist2_dE_q=f2.Get("delta_eta_quark")
hist2_dP_q=f2.Get("dPhi_quark")
hist2_pT_q=f2.Get("quark_pt")
hist2_num_q=f2.Get("number_of_quarks")
while dummy == 0:
   jet2hist[jet2counter+1] = f2.Get(f'jet{jet2counter+1}')
   if not jet2hist[jet2counter+1]:
      dummy = 1
      print("File 2 has ",jet2counter," jets")
   else:
      jet2counter = jet2counter + 1
dummy=0
##############
#f4 = ROOT.TFile(rootfile4)
#hist4mjj = f4.Get('jet_mJJ')
#hist4delEta = f4.Get("delta_eta")
#hist4piPT = f4.Get("pion_pt")
#hist4pinum=f4.Get("pion_number")
#hist4met=f4.Get("jet_met")
#while dummy == 0:
#   jet4hist[jet4counter+1] = f4.Get(f'jet{jet4counter+1}')
#   if not jet4hist[jet4counter+1]:
#      dummy = 1
#      print("File 4 has ",jet4counter," jets")
#   else:
#      jet4counter = jet4counter + 1
#dummy=0
##############

#jetcounter = min(jet1counter,jet2counter,jet3counter,jet4counter)
jetcounter = min(jet1counter,jet2counter)
print("Number of jets containing all files is ", jetcounter)

###### LEADING JET DISTRIBUTION ######
canvas1 = ROOT.TCanvas("canvas1")
canvas1.cd()
jet1hist[1].SetTitle("p_{T} distribution of leading jet")
jet1hist[1].GetXaxis().SetTitle("p_{T}[GeV]")
jet1hist[1].GetYaxis().SetTitle("A.U.")
jet1hist[1].GetXaxis().SetLabelSize(0.04)
jet1hist[1].GetYaxis().SetLabelSize(0.04)
jet1hist[1].GetXaxis().SetTitleSize(0.04)
jet1hist[1].GetYaxis().SetTitleSize(0.04)
jet1hist[1].SetStats(0)
jet1hist[1].SetLineColor(2)
jet1hist[1].SetLineWidth(2)
jet1hist[1].Draw("hist1 E")
jet1hist[1].GetXaxis().SetRangeUser(0,3*jet1hist[1].GetMean())
jet1hist[1].GetYaxis().SetRangeUser(0,0.15)
jet1hist[1].GetXaxis().SetTitleOffset(1.4)
jet1hist[1].GetYaxis().SetTitleOffset(1.4)
jet2hist[1].SetLineColor(1)
jet2hist[1].SetLineWidth(2)
jet2hist[1].SetStats(0)
jet2hist[1].Draw("hist1 E same")
#jet4hist[1].SetLineColor(1)
#jet4hist[1].SetLineWidth(2)
#jet4hist[1].SetStats(0)
#jet4hist[1].Draw("hist1 E same")
t1 = ROOT.TLatex()
t1.SetNDC(ROOT.kTRUE)
t1.SetTextSize(0.03)
t1.DrawLatex(0.18 ,0.86 ,"\Lambda = %s, %s decay, N_{f}=%s, N_{c}=%s"%(const1,const2,const3,const4))
t1.SetTextSize(0.03)
legend1 = ROOT.TLegend(0.50,0.65,0.9,0.9)
legend1.SetTextSize(0.03)
legend1.AddEntry(jet1hist[1],"R=0.4, Mean = %.1f GeV"%(meancount[1]),"l")
legend1.AddEntry(jet2hist[1],"R=1.0, Mean = %.1f GeV"%(jet2hist[1].GetMean()),"l")
#legend1.AddEntry(jet4hist[1],"N_{f} = 8 , Mean = %.1f GeV"%(jet4hist[1].GetMean()),"l")
legend1.SetLineWidth(1)
legend1.Draw()
canvas1.SetRightMargin(0.09)
canvas1.SetLeftMargin(0.15)
canvas1.SetBottomMargin(0.15)
canvas1.Update()
ROOT.gSystem.ProcessEvents()
canvas1.Print("pT_leading_jet_varyingradius.pdf")

###### SUBLEADING JET DISTRIBUTION ######
canvas2 = ROOT.TCanvas("canvas2")
canvas2.cd()
jet1hist[2].SetTitle("p_{T} distribution of subleading jet")
jet1hist[2].GetXaxis().SetTitle("p_{T}[GeV]")
jet1hist[2].GetYaxis().SetTitle("A.U.")
jet1hist[2].GetXaxis().SetLabelSize(0.04)
jet1hist[2].GetYaxis().SetLabelSize(0.04)
jet1hist[2].GetXaxis().SetTitleSize(0.04)
jet1hist[2].GetYaxis().SetTitleSize(0.04)
jet1hist[2].SetStats(0)
jet1hist[2].SetLineColor(2)
jet1hist[2].SetLineWidth(2)
jet1hist[2].Draw("hist2 E")
jet1hist[2].GetXaxis().SetRangeUser(0,3*jet1hist[2].GetMean())
jet1hist[2].GetYaxis().SetRangeUser(0,0.3)
jet1hist[2].GetXaxis().SetTitleOffset(1.4)
jet1hist[2].GetYaxis().SetTitleOffset(1.4)
jet2hist[2].SetLineColor(1)
jet2hist[2].SetLineWidth(2)
jet2hist[2].SetStats(0)
jet2hist[2].Draw("hist2 E same")
#jet4hist[2].SetLineColor(1)
#jet4hist[2].SetLineWidth(2)
#jet4hist[2].SetStats(0)
#jet4hist[2].Draw("hist2 E same")
t2 = ROOT.TLatex()
t2.SetNDC(ROOT.kTRUE)
t2.SetTextSize(0.03)
t2.DrawLatex(0.18 ,0.86 ,"\Lambda = %s, %s decay, N_{f}=%s, N_{c}=%s"%(const1,const2,const3,const4))
t2.SetTextSize(0.03)
legend2 = ROOT.TLegend(0.55,0.65,0.9,0.9)
legend2.SetTextSize(0.03)
legend2.AddEntry(jet1hist[2],"R=0.4, Mean = %.1f GeV"%(meancount[2]),"l")
legend2.AddEntry(jet2hist[2],"R=1.0, Mean = %.1f GeV"%(jet2hist[2].GetMean()),"l")
#legend2.AddEntry(jet4hist[2],"N_{f} = 8 , Mean = %.1f GeV"%(jet4hist[2].GetMean()),"l")
legend2.SetLineWidth(1)
legend2.Draw()
canvas2.SetRightMargin(0.09)
canvas2.SetLeftMargin(0.15)
canvas2.SetBottomMargin(0.15)
canvas2.Update()
ROOT.gSystem.ProcessEvents()
canvas2.Print("pT_subleading_jet_varyingradius.pdf")


###### ETA DISTRIBUTION ######
meaneta = hist1delEta.GetMean()
canvas3 = ROOT.TCanvas("canvas3")
canvas3.cd()
hist1delEta.SetTitle("Distribution of \Delta\eta for leading and sub-leading jets")
hist1delEta.GetXaxis().SetTitle("\Delta\eta")
hist1delEta.GetYaxis().SetTitle("A.U.")
hist1delEta.GetXaxis().SetLabelSize(0.04)
hist1delEta.GetYaxis().SetLabelSize(0.04)
hist1delEta.GetXaxis().SetTitleSize(0.04)
hist1delEta.GetYaxis().SetTitleSize(0.04)
hist1delEta.SetStats(0)
hist1delEta.SetLineColor(2)
hist1delEta.SetLineWidth(2)
hist1delEta.Draw("hist3 E")
hist1delEta.GetXaxis().SetRangeUser(-5,5)
hist1delEta.GetYaxis().SetRangeUser(0,0.1)
hist1delEta.GetXaxis().SetTitleOffset(1.4)
hist1delEta.GetYaxis().SetTitleOffset(1.4)
hist2delEta.SetLineColor(1)
hist2delEta.SetLineWidth(2)
hist2delEta.SetStats(0)
hist2delEta.Draw("hist3 E same")
#hist4delEta.SetLineColor(1)
#hist4delEta.SetLineWidth(2)
#hist4delEta.SetStats(0)
#hist4delEta.Draw("hist3 E same")
canvas3.SetRightMargin(0.09)
canvas3.SetLeftMargin(0.15)
canvas3.SetBottomMargin(0.15)
t3 = ROOT.TLatex()
t3.SetNDC(ROOT.kTRUE)
t3.SetTextSize(0.03)
t3.DrawLatex(0.18 ,0.86 ,"\Lambda = %s, %s decay, N_{f}=%s, N_{c}=%s"%(const1,const2,const3,const4))
t3.SetTextSize(0.03)
legend3 = ROOT.TLegend(0.50,0.65,0.9,0.9)
legend3.SetTextSize(0.03)
legend3.AddEntry(hist1delEta,"R=0.4, Mean = %.3f"%(meaneta),"l")
legend3.AddEntry(hist2delEta,"R=1.0, Mean = %.3f"%(hist2delEta.GetMean()),"l")
#legend3.AddEntry(hist4delEta,"N_{f} = 8 , Mean = %.1f"%(hist4delEta.GetMean()),"l")
legend3.SetLineWidth(1)
legend3.Draw()
canvas3.Update()
ROOT.gSystem.ProcessEvents()
canvas3.Print("delta_eta_varyingradius.pdf")


###### INVARIANT MASS DISTRIBUTION ######
canvas4 = ROOT.TCanvas("canvas4")
canvas4.cd()
meanmjj = hist1mjj.GetMean()
hist1mjj.SetTitle("Invariant mass distribution")
hist1mjj.GetXaxis().SetTitle("Invariant mass [GeV]")
hist1mjj.GetYaxis().SetTitle("A.U.")
hist1mjj.GetXaxis().SetLabelSize(0.04)
hist1mjj.GetYaxis().SetLabelSize(0.04)
hist1mjj.GetXaxis().SetTitleSize(0.04)
hist1mjj.GetYaxis().SetTitleSize(0.04)
hist1mjj.SetStats(0)
hist1mjj.SetLineColor(2)
hist1mjj.SetLineWidth(2)
hist1mjj.Draw("hist4 E")
hist2mjj.SetLineColor(1)
hist2mjj.SetLineWidth(2)
hist2mjj.SetStats(0)
hist2mjj.Draw("hist4 E same")
hist1mjj.GetXaxis().SetRangeUser(0,600)
hist1mjj.GetYaxis().SetRangeUser(0,0.12)
hist1mjj.GetXaxis().SetTitleOffset(1.4)
hist1mjj.GetYaxis().SetTitleOffset(1.4)
#hist4mjj.SetLineColor(1)
#hist4mjj.SetLineWidth(2)
#hist4mjj.SetStats(0)
#hist4mjj.Draw("hist4 E same")
canvas4.SetRightMargin(0.09)
canvas4.SetLeftMargin(0.15)
canvas4.SetBottomMargin(0.15)
t4 = ROOT.TLatex()
t4.SetNDC(ROOT.kTRUE)
t4.SetTextSize(0.03)
t4.DrawLatex(0.18 ,0.86 ,"\Lambda = %s, %s decay, N_{f}=%s, N_{c}=%s"%(const1,const2,const3,const4))
t4.SetTextSize(0.03)
legend4 = ROOT.TLegend(0.50,0.65,0.9,0.9)
legend4.SetTextSize(0.03)
legend4.AddEntry(hist1mjj,"R=0.4, Mean = %.1f GeV"%(meanmjj),"l")
legend4.AddEntry(hist2mjj,"R=1.0, Mean = %.1f GeV"%(hist2mjj.GetMean()),"l")
#legend4.AddEntry(hist4mjj,"N_{f} = 8 , Mean = %.1f GeV"%(hist4mjj.GetMean()),"l")
legend4.SetLineWidth(1)
legend4.Draw()
canvas4.Update()
ROOT.gSystem.ProcessEvents()
canvas4.Print("invariant_mass_varyingradius.pdf")


###### MISSING TRANSVERSE ENERGY (MET) DISTRIBUTION ######
canvas5 = ROOT.TCanvas("canvas5")
canvas5.cd()
metmean = hist1met.GetMean()
hist1met.SetTitle("#slash{E}_{T} distribution")
hist1met.GetXaxis().SetTitle("#slash{E}_{T} [GeV]")
hist1met.GetYaxis().SetTitle("A.U.")
hist1met.GetXaxis().SetLabelSize(0.04)
hist1met.GetYaxis().SetLabelSize(0.04)
hist1met.GetXaxis().SetTitleSize(0.04)
hist1met.GetYaxis().SetTitleSize(0.04)
hist1met.SetStats(0)
hist1met.SetLineColor(2)
hist1met.SetLineWidth(2)
hist1met.Draw("hist5 E")
hist2met.SetLineColor(1)
hist2met.SetLineWidth(2)
hist2met.SetStats(0)
hist2met.Draw("hist5 E same")
hist1met.GetXaxis().SetRangeUser(0,300)
hist1met.GetYaxis().SetRangeUser(0,0.15)
hist1met.GetXaxis().SetTitleOffset(1.4)
hist1met.GetYaxis().SetTitleOffset(1.4)
#hist4met.SetLineColor(1)
#hist4met.SetLineWidth(2)
#hist4met.SetStats(0)
#hist4met.Draw("hist5 E same")
canvas5.SetRightMargin(0.09)
canvas5.SetLeftMargin(0.15)
canvas5.SetBottomMargin(0.15)
t5 = ROOT.TLatex()
t5.SetNDC(ROOT.kTRUE)
t5.SetTextSize(0.03)
t5.DrawLatex(0.18 ,0.86 ,"\Lambda = %s, %s decay, N_{f}=%s, N_{c}=%s"%(const1,const2,const3,const4))
t5.SetTextSize(0.03)
legend5 = ROOT.TLegend(0.50,0.65,0.9,0.9)
legend5.SetTextSize(0.03)
legend5.AddEntry(hist1met,"R=0.4, Mean = %.1f GeV"%(metmean),"l")
legend5.AddEntry(hist2met,"R=1.0, Mean = %.1f GeV"%(hist2met.GetMean()),"l")
#legend5.AddEntry(hist4met,"N_{f} = 8 , Mean = %.1f GeV"%(hist4met.GetMean()),"l")
legend5.SetLineWidth(1)
legend5.Draw()
canvas5.Update()
ROOT.gSystem.ProcessEvents()
canvas5.Print("jet_met_varyingradius.pdf")

###### PION PT DISTRIBUTION ######
canvas6 = ROOT.TCanvas("canvas6")
canvas6.cd()
meanpiPT = hist1piPT.GetMean()
hist1piPT.SetTitle("Dark pion p_{T} distribution")
hist1piPT.GetXaxis().SetTitle("p_{T} [GeV]")
hist1piPT.GetYaxis().SetTitle("A.U.")
hist1piPT.GetXaxis().SetLabelSize(0.04)
hist1piPT.GetYaxis().SetLabelSize(0.04)
hist1piPT.GetXaxis().SetTitleSize(0.04)
hist1piPT.GetYaxis().SetTitleSize(0.04)
hist1piPT.SetStats(0)
hist1piPT.SetLineColor(2)
hist1piPT.SetLineWidth(2)
hist1piPT.Draw("hist6 E")
hist2piPT.SetLineColor(1)
hist2piPT.SetLineWidth(2)
hist2piPT.SetStats(0)
hist2piPT.Draw("hist6 E same")
hist1piPT.GetXaxis().SetRangeUser(0,120)
hist1piPT.GetYaxis().SetRangeUser(0,0.4)
hist1piPT.GetXaxis().SetTitleOffset(1.4)
hist1piPT.GetYaxis().SetTitleOffset(1.4)
#hist4piPT.SetLineColor(1)
#hist4piPT.SetLineWidth(2)
#hist4piPT.SetStats(0)
#hist4piPT.Draw("hist6 E same")
canvas6.SetRightMargin(0.09)
canvas6.SetLeftMargin(0.15)
canvas6.SetBottomMargin(0.15)
t6 = ROOT.TLatex()
t6.SetNDC(ROOT.kTRUE)
t6.SetTextSize(0.03)
t6.DrawLatex(0.18 ,0.86 ,"\Lambda = %s, %s decay, N_{f}=%s, N_{c}=%s"%(const1,const2,const3,const4))
t6.SetTextSize(0.03)
legend6 = ROOT.TLegend(0.55,0.65,0.9,0.9)
legend6.SetTextSize(0.03)
legend6.AddEntry(hist1piPT,"R=0.4, Mean = %.1f GeV"%(meanpiPT),"l")
legend6.AddEntry(hist2piPT,"R=1.0, Mean = %.1f GeV"%(hist2piPT.GetMean()),"l")
#legend6.AddEntry(hist4piPT,"N_{f} = 8 , Mean = %.1f GeV"%(hist4piPT.GetMean()),"l")
legend6.SetLineWidth(1)
legend6.Draw()
canvas6.Update()
ROOT.gSystem.ProcessEvents()
canvas6.Print("pion_pT_varyingradius.pdf")


###### PION NUMBER DISTRIBUTION ######
canvas7 = ROOT.TCanvas("canvas7")
canvas7.cd()
meanpinum = hist1pinum.GetMean()
hist1pinum.SetTitle("Dark pion number distribution")
hist1pinum.GetXaxis().SetTitle("Number of dark pions")
hist1pinum.GetYaxis().SetTitle("A.U.")
hist1pinum.GetXaxis().SetLabelSize(0.04)
hist1pinum.GetYaxis().SetLabelSize(0.04)
hist1pinum.GetXaxis().SetTitleSize(0.04)
hist1pinum.GetYaxis().SetTitleSize(0.04)
hist1pinum.SetStats(0)
hist1pinum.SetLineColor(2)
hist1pinum.SetLineWidth(2)
hist1pinum.Draw("hist7 E")
hist2pinum.SetLineColor(1)
hist2pinum.SetLineWidth(2)
hist2pinum.SetStats(0)
hist2pinum.Draw("hist7 E same")
hist1pinum.GetXaxis().SetRangeUser(0,60)
#hist1pinum.GetYaxis().SetRangeUser(0,0.30)######0.15 for nf,nc
hist1pinum.GetYaxis().SetRangeUser(0,0.2)
hist1pinum.GetXaxis().SetTitleOffset(1.4)
hist1pinum.GetYaxis().SetTitleOffset(1.4)
#hist4pinum.SetLineColor(1)
#hist4pinum.SetLineWidth(2)
#hist4pinum.SetStats(0)
#hist4pinum.Draw("hist7 E same")
canvas7.SetRightMargin(0.09)
canvas7.SetLeftMargin(0.15)
canvas7.SetBottomMargin(0.15)
t7 = ROOT.TLatex()
t7.SetNDC(ROOT.kTRUE)
t7.SetTextSize(0.03)
t7.DrawLatex(0.18 ,0.86 ,"\Lambda = %s, %s decay, N_{f}=%s, N_{c}=%s"%(const1,const2,const3,const4))
t7.SetTextSize(0.03)
legend7 = ROOT.TLegend(0.55,0.65,0.9,0.9)
legend7.SetTextSize(0.03)
legend7.AddEntry(hist1pinum,"R=0.4, Mean = %.1f"%(meanpinum),"l")
legend7.AddEntry(hist2pinum,"R=1.0, Mean = %.1f"%(hist2pinum.GetMean()),"l")
#legend7.AddEntry(hist4pinum,"N_{f} = 8 , Mean = %.1f"%(hist4pinum.GetMean()),"l")
legend7.SetLineWidth(1)
legend7.Draw()
canvas7.Update()
ROOT.gSystem.ProcessEvents()
canvas7.Print("pion_number_varyingradius.pdf")


###### DELTA PHI DISTRIBUTION ######
canvas8 = ROOT.TCanvas("canvas8")
canvas8.cd()
meandelphi = hist1dPhi.GetMean()
hist1dPhi.SetTitle("Distribution of \Delta\phi for leading and sub-leading jets")
hist1dPhi.GetXaxis().SetTitle("\Delta\phi")
hist1dPhi.GetYaxis().SetTitle("A.U.")
hist1dPhi.GetXaxis().SetLabelSize(0.04)
hist1dPhi.GetYaxis().SetLabelSize(0.04)
hist1dPhi.GetXaxis().SetTitleSize(0.04)
hist1dPhi.GetYaxis().SetTitleSize(0.04)
hist1dPhi.SetStats(0)
hist1dPhi.SetLineColor(2)
hist1dPhi.SetLineWidth(2)
hist1dPhi.Draw("hist8 E")
hist2dPhi.SetLineColor(1)
hist2dPhi.SetLineWidth(2)
hist2dPhi.SetStats(0)
hist2dPhi.Draw("hist8 E same")
hist1dPhi.GetXaxis().SetRangeUser(-3.14,3.14)
hist1dPhi.GetYaxis().SetRangeUser(0,0.2)
hist1dPhi.GetXaxis().SetTitleOffset(1.4)
hist1dPhi.GetYaxis().SetTitleOffset(1.4)
#hist4pinum.SetLineColor(1)
#hist4pinum.SetLineWidth(2)
#hist4pinum.SetStats(0)
#hist4pinum.Draw("hist8 E same")
canvas8.SetRightMargin(0.09)
canvas8.SetLeftMargin(0.15)
canvas8.SetBottomMargin(0.15)
t8 = ROOT.TLatex()
t8.SetNDC(ROOT.kTRUE)
t8.SetTextSize(0.03)
t8.DrawLatex(0.18 ,0.86 ,"\Lambda = %s, %s decay, N_{f}=%s, N_{c}=%s"%(const1,const2,const3,const4))
t8.SetTextSize(0.03)
legend8 = ROOT.TLegend(0.55,0.65,0.9,0.9)
legend8.SetTextSize(0.03)
legend8.AddEntry(hist1dPhi,"R=0.4, Mean = %.3f"%(meandelphi),"l")
legend8.AddEntry(hist2dPhi,"R=1.0, Mean = %.3f"%(hist2dPhi.GetMean()),"l")
#legend7.AddEntry(hist4pinum,"N_{f} = 8 , Mean = %.1f"%(hist4pinum.GetMean()),"l")
legend8.SetLineWidth(1)
legend8.Draw()
canvas8.Update()
ROOT.gSystem.ProcessEvents()
canvas8.Print("del_phi_varyingradius.pdf")

###### SM QUARK NUMBER DISTRIBUTION ######
canvas9 = ROOT.TCanvas("canvas9")
canvas9.cd()
mean_num_q = hist1_num_q.GetMean()
hist1_num_q.SetTitle("SM quark number distribution")
hist1_num_q.GetXaxis().SetTitle("Number of SM quarks")
hist1_num_q.GetYaxis().SetTitle("A.U.")
hist1_num_q.GetXaxis().SetLabelSize(0.04)
hist1_num_q.GetYaxis().SetLabelSize(0.04)
hist1_num_q.GetXaxis().SetTitleSize(0.04)
hist1_num_q.GetYaxis().SetTitleSize(0.04)
hist1_num_q.SetStats(0)
hist1_num_q.SetLineColor(2)
hist1_num_q.SetLineWidth(2)
hist1_num_q.Draw("hist10 E")
hist2_num_q.SetLineColor(1)
hist2_num_q.SetLineWidth(2)
hist2_num_q.SetStats(0)
hist2_num_q.Draw("hist10 E same")
hist1_num_q.GetXaxis().SetRangeUser(0,30)
hist1_num_q.GetYaxis().SetRangeUser(0,1)
hist1_num_q.GetXaxis().SetTitleOffset(1.4)
hist1_num_q.GetYaxis().SetTitleOffset(1.4)
#hist4pinum.SetLineColor(1)
#hist4pinum.SetLineWidth(2)
#hist4pinum.SetStats(0)
#hist4pinum.Draw("hist10 E same")
canvas9.SetRightMargin(0.09)
canvas9.SetLeftMargin(0.15)
canvas9.SetBottomMargin(0.15)
t9 = ROOT.TLatex()
t9.SetNDC(ROOT.kTRUE)
t9.SetTextSize(0.03)
t9.DrawLatex(0.18 ,0.86 ,"\Lambda = %s, %s decay, N_{f}=%s, N_{c}=%s"%(const1,const2,const3,const4))
t9.SetTextSize(0.03)
legend9 = ROOT.TLegend(0.55,0.65,0.9,0.9)
legend9.SetTextSize(0.03)
legend9.AddEntry(hist1_num_q,"R=0.4, Mean = %.1f"%(mean_num_q),"l")
legend9.AddEntry(hist2_num_q,"R=1.0, Mean = %.1f"%(hist2_num_q.GetMean()),"l")
#legend7.AddEntry(hist4pinum,"N_{f} = 8 , Mean = %.1f"%(hist4pinum.GetMean()),"l")
legend9.SetLineWidth(1)
legend9.Draw()
canvas9.Update()
ROOT.gSystem.ProcessEvents()
canvas9.Print("quark_num_varyingradius.pdf")

###### SM QUARK PT DISTRIBUTION ######
canvas10 = ROOT.TCanvas("canvas10")
canvas10.cd()
mean_pt_q = hist1_pT_q.GetMean()
hist1_pT_q.SetTitle("SM quark p_{T} distribution")
hist1_pT_q.GetXaxis().SetTitle("p_{T}[GeV]")
hist1_pT_q.GetYaxis().SetTitle("A.U.")
hist1_pT_q.GetXaxis().SetLabelSize(0.04)
hist1_pT_q.GetYaxis().SetLabelSize(0.04)
hist1_pT_q.GetXaxis().SetTitleSize(0.04)
hist1_pT_q.GetYaxis().SetTitleSize(0.04)
hist1_pT_q.SetStats(0)
hist1_pT_q.SetLineColor(2)
hist1_pT_q.SetLineWidth(2)
hist1_pT_q.Draw("hist10 E")
hist2_pT_q.SetLineColor(1)
hist2_pT_q.SetLineWidth(2)
hist2_pT_q.SetStats(0)
hist2_pT_q.Draw("hist10 E same")
hist1_pT_q.GetXaxis().SetRangeUser(0,50)
hist1_pT_q.GetYaxis().SetRangeUser(0,1)
hist1_pT_q.GetXaxis().SetTitleOffset(1.4)
hist1_pT_q.GetYaxis().SetTitleOffset(1.4)
#hist4pinum.SetLineColor(1)
#hist4pinum.SetLineWidth(2)
#hist4pinum.SetStats(0)
#hist4pinum.Draw("hist10 E same")
canvas10.SetRightMargin(0.09)
canvas10.SetLeftMargin(0.15)
canvas10.SetBottomMargin(0.15)
t10 = ROOT.TLatex()
t10.SetNDC(ROOT.kTRUE)
t10.SetTextSize(0.03)
t10.DrawLatex(0.18 ,0.86 ,"\Lambda = %s, %s decay, N_{f}=%s, N_{c}=%s"%(const1,const2,const3,const4))
t10.SetTextSize(0.03)
legend10 = ROOT.TLegend(0.55,0.65,0.9,0.9)
legend10.SetTextSize(0.03)
legend10.AddEntry(hist1_pT_q,"R=0.4, Mean = %.1f GeV"%(mean_pt_q),"l")
legend10.AddEntry(hist2_pT_q,"R=1.0, Mean = %.1f GeV"%(hist2_pT_q.GetMean()),"l")
#legend7.AddEntry(hist4pinum,"N_{f} = 8 , Mean = %.1f"%(hist4pinum.GetMean()),"l")
legend10.SetLineWidth(1)
legend10.Draw()
canvas10.Update()
ROOT.gSystem.ProcessEvents()
canvas10.Print("quark_pt_varyingradius.pdf")

###### SM QUARK DELTA ETA DISTRIBUTION ######
canvas11 = ROOT.TCanvas("canvas11")
canvas11.cd()
mean_dE_q = hist1_dE_q.GetMean()
hist1_dE_q.SetTitle("Distribution of \Delta\eta of two SM quarks")
hist1_dE_q.GetXaxis().SetTitle("\Delta\eta")
hist1_dE_q.GetYaxis().SetTitle("A.U.")
hist1_dE_q.GetXaxis().SetLabelSize(0.04)
hist1_dE_q.GetYaxis().SetLabelSize(0.04)
hist1_dE_q.GetXaxis().SetTitleSize(0.04)
hist1_dE_q.GetYaxis().SetTitleSize(0.04)
hist1_dE_q.SetStats(0)
hist1_dE_q.SetLineColor(2)
hist1_dE_q.SetLineWidth(2)
hist1_dE_q.Draw("hist11 E")
hist2_dE_q.SetLineColor(1)
hist2_dE_q.SetLineWidth(2)
hist2_dE_q.SetStats(0)
hist2_dE_q.Draw("hist11 E same")
hist1_dE_q.GetXaxis().SetRangeUser(-5,5)
hist1_dE_q.GetYaxis().SetRangeUser(0,0.25)
hist1_dE_q.GetXaxis().SetTitleOffset(1.4)
hist1_dE_q.GetYaxis().SetTitleOffset(1.4)
#hist4pinum.SetLineColor(1)
#hist4pinum.SetLineWidth(2)
#hist4pinum.SetStats(0)
#hist4pinum.Draw("hist10 E same")
canvas11.SetRightMargin(0.09)
canvas11.SetLeftMargin(0.15)
canvas11.SetBottomMargin(0.15)
t11 = ROOT.TLatex()
t11.SetNDC(ROOT.kTRUE)
t11.SetTextSize(0.03)
t11.DrawLatex(0.18 ,0.86 ,"\Lambda = %s, %s decay, N_{f}=%s, N_{c}=%s"%(const1,const2,const3,const4))
t11.SetTextSize(0.03)
legend11 = ROOT.TLegend(0.55,0.65,0.9,0.9)
legend11.SetTextSize(0.03)
legend11.AddEntry(hist1_dE_q,"R=0.4, Mean = %.3f"%(mean_dE_q),"l")
legend11.AddEntry(hist2_dE_q,"R=1.0, Mean = %.3f"%(hist2_dE_q.GetMean()),"l")
#legend7.AddEntry(hist4pinum,"N_{f} = 8 , Mean = %.1f"%(hist4pinum.GetMean()),"l")
legend11.SetLineWidth(1)
legend11.Draw()
canvas11.Update()
ROOT.gSystem.ProcessEvents()
canvas11.Print("quark_del_eta_varyingradius.pdf")

###### SM QUARK DELTA PHI DISTRIBUTION ######
canvas12 = ROOT.TCanvas("canvas12")
canvas12.cd()
mean_dP_q = hist1_dP_q.GetMean()
hist1_dP_q.SetTitle("Distribution of \Delta\phi of two SM quarks")
hist1_dP_q.GetXaxis().SetTitle("\Delta\phi")
hist1_dP_q.GetYaxis().SetTitle("A.U.")
hist1_dP_q.GetXaxis().SetLabelSize(0.04)
hist1_dP_q.GetYaxis().SetLabelSize(0.04)
hist1_dP_q.GetXaxis().SetTitleSize(0.04)
hist1_dP_q.GetYaxis().SetTitleSize(0.04)
hist1_dP_q.SetStats(0)
hist1_dP_q.SetLineColor(2)
hist1_dP_q.SetLineWidth(2)
hist1_dP_q.Draw("hist12 E")
hist2_dP_q.SetLineColor(1)
hist2_dP_q.SetLineWidth(2)
hist2_dP_q.SetStats(0)
hist2_dP_q.Draw("hist12 E same")
hist1_dP_q.GetXaxis().SetRangeUser(-3.14,3.14)
hist1_dP_q.GetYaxis().SetRangeUser(0,0.2)
hist1_dP_q.GetXaxis().SetTitleOffset(1.4)
hist1_dP_q.GetYaxis().SetTitleOffset(1.4)
#hist4pinum.SetLineColor(1)
#hist4pinum.SetLineWidth(2)
#hist4pinum.SetStats(0)
#hist4pinum.Draw("hist10 E same")
canvas12.SetRightMargin(0.09)
canvas12.SetLeftMargin(0.15)
canvas12.SetBottomMargin(0.15)
t12 = ROOT.TLatex()
t12.SetNDC(ROOT.kTRUE)
t12.SetTextSize(0.03)
t12.DrawLatex(0.18 ,0.86 ,"\Lambda = %s, %s decay, N_{f}=%s, N_{c}=%s"%(const1,const2,const3,const4))
t12.SetTextSize(0.03)
legend12 = ROOT.TLegend(0.55,0.65,0.9,0.9)
legend12.SetTextSize(0.03)
legend12.AddEntry(hist1_dP_q,"R=0.4, Mean = %.3f"%(mean_dP_q),"l")
legend12.AddEntry(hist2_dP_q,"R=1.0, Mean = %.3f"%(hist2_dP_q.GetMean()),"l")
#legend7.AddEntry(hist4pinum,"N_{f} = 8 , Mean = %.1f"%(hist4pinum.GetMean()),"l")
legend12.SetLineWidth(1)
legend12.Draw()
canvas12.Update()
ROOT.gSystem.ProcessEvents()
canvas12.Print("quark_del_phi_varyingradius.pdf")
