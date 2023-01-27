# Authors and contacts:
# Guillaume Albouy: guillaume.albouy@etu.univ-grenoble-alpes.fr
# Akanksha Singh: akki153209@gmail.com
# Harikrishnan Nair: hunair1996@gmail.com

# This code needs helpers.py and histo_defs.py which are stored in the same folder as this file
# Following cuts on events are assumed
# jet1.PT > 500 and np.abs(jet1.Eta) < 2.5 and jet2.PT > 500 and np.abs(jet2.Eta) < 2.5
# We assume that different branches corresponding to different jet radii are defined in root files
# In this code jet clustering radius is fixed to 1.4, can be changed at line 52
# This code will analyse input sample and create following reconstructed level normalized distributions:
# 1) pt of leading/subleading jet
# 2) dijet invarint mass
# 3) missing energy
# 4) transverse mass of dijet and met system
# 5) transverse ratio
# 6) delta phi between missin energy and leading/subleading jet
# 7) 2D histo of track pT of leading jet
# 8) delta eta between leading and subleading jet
# 9) pt and invariant mass for trimmed and SoftDropped leading/subleading jets
#
# command: python /path_of_code/transverse_mass.py /path_of_rootfile/name_of_rootfile.root /path_of_rootfile/output_name.root
# Takes the rootfile as input and computes the defined variables and fills the respective histograms.


from math import sin
import sys
import numpy as np
import ROOT
from array import array
from helpers import *


try:
  input = raw_input
except:
  pass

if len(sys.argv) < 3:
  print(" Usage: Analysis_code/Jet_analysis.py /path/delphes_file.root /path/output.root")
  sys.exit(1)

ROOT.gSystem.Load("libDelphes")

try:
        ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
        ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
        pass

# Parameters
############################################
# Radius of jets (0.4, 1.0, 1.4) :
R = 1.4
# Events selection (pT in GeV)
pT_min_jet1 = 500
pT_min_jet2 = 500
eta_max = 2.5
M_PI = 3.14
############################################

inputFile = sys.argv[1]
print("Input file :")
print(inputFile)

outputFile = sys.argv[2]
print("output file :")
print(outputFile)

# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

# Get pointer to branches used in this analysis
# R-jet branches : 04, 08, 10, 12, 14
R_jet = str(int(R*10))
if R<1.0: R_jet = '0' + R_jet

# Getting the required branches from Delphes ROOT file.
#branchJet = treeReader.UseBranch("ParticleFlowJet%s"%R_jet)
#branchMET = treeReader.UseBranch("GenMissingET")
branchMET = treeReader.UseBranch("MissingET")
branchtrack = treeReader.UseBranch("Track")
branchJet = treeReader.UseBranch("Jet")
branchPtcl = treeReader.UseBranch("Particle")

####BOOKING ALL (NON-LOOPED) HISTOGRAMS HERE#####
###JetMet angles
Nbins=10
histquarknum = ROOT.TH1F("number_of_quarks", "Number of quark decay products", Nbins, 0.0, 30.0)
Nbins=20
hdphi1 = ROOT.TH1F("dPhi1", "dPhi_jet1_met", Nbins, 0.0, 10.0)
hdphi2 = ROOT.TH1F("dPhi2", "dPhi_jet2_met", Nbins, 0.0, 10.0)
###Phi/Eta overlay distributions
Nbins=25
histdelEta_positive=ROOT.TH1F("delta_eta_positive", "delta_eta_positive", Nbins,0,5.0)
histdelEta_negative=ROOT.TH1F("delta_eta_negative", "delta_eta_negative", Nbins,0,5.0)
histdelPhi_positive=ROOT.TH1F("delta_phi_positive", "delta_phi_positive", Nbins,0,3.15)
histdelPhi_negative=ROOT.TH1F("delta_phi_negative", "delta_phi_negative", Nbins,0,3.15)
histdelEta_positive_quark=ROOT.TH1F("delta_eta_positive_quark", "delta_eta_positive_quark", Nbins,0,5.0)
histdelEta_negative_quark=ROOT.TH1F("delta_eta_negative_quark", "delta_eta_negative_quark", Nbins,0,5.0)
histdelPhi_positive_quark=ROOT.TH1F("delta_phi_positive_quark", "delta_phi_positive_quark", Nbins,0,3.15)
histdelPhi_negative_quark=ROOT.TH1F("delta_phi_negative_quark", "delta_phi_negative_quark", Nbins,0,3.15)
histdelEta_positive_darkquark=ROOT.TH1F("delta_eta_positive_darkquark", "delta_eta_positive_quark", Nbins,0,5.0)
histdelEta_negative_darkquark=ROOT.TH1F("delta_eta_negative_darkquark", "delta_eta_negative_quark", Nbins,0,5.0)
histdelPhi_positive_darkquark=ROOT.TH1F("delta_phi_positive_darkquark", "delta_phi_positive_quark", Nbins,0,3.15)
histdelPhi_negative_darkquark=ROOT.TH1F("delta_phi_negative_darkquark", "delta_phi_negative_quark", Nbins,0,3.15)
histdelEta_positive_darkquark_initial=ROOT.TH1F("delta_eta_positive_darkquark_initial", "delta_eta_positive_quark", Nbins,0,5.0)
histdelEta_negative_darkquark_initial=ROOT.TH1F("delta_eta_negative_darkquark_initial", "delta_eta_negative_quark", Nbins,0,5.0)
histdelPhi_positive_darkquark_initial=ROOT.TH1F("delta_phi_positive_darkquark_initial", "delta_phi_positive_quark", Nbins,0,3.15)
histdelPhi_negative_darkquark_initial=ROOT.TH1F("delta_phi_negative_darkquark_initial", "delta_phi_negative_quark", Nbins,0,3.15)
histdarkquarknum = ROOT.TH1F("number_of_dark_quarks", "Number of quark decay products", Nbins, 0.0, 100.0)
histdarkquarknum_initial = ROOT.TH1F("number_of_initial_dark_quarks", "Number of quark decay products", Nbins, 0.0, 25.0)
#Phi/Eta full distributions
Nbins = 50
histdelEta = ROOT.TH1F("delta_eta", "delta_eta_jet", Nbins,-5.0,5.0)
histdelPhi = ROOT.TH1F("dPhi", "dPhi_jet1", Nbins, -3.15, 3.15)
histdelEta_quark = ROOT.TH1F("delta_eta_quark", "delta_eta_jet_quark", Nbins,-5.0,5.0)
histdelPhi_quark = ROOT.TH1F("dPhi_quark", "dPhi_jet1_quark", Nbins, -3.15, 3.15)
histdelEta_darkquark = ROOT.TH1F("delta_eta_darkquark", "delta_eta_jet_quark", Nbins,-5.0,5.0)
histdelPhi_darkquark = ROOT.TH1F("dPhi_darkquark", "dPhi_jet1_quark", Nbins, -3.15, 3.15)
histdelEta_darkquark_initial = ROOT.TH1F("delta_eta_darkquark_initial", "delta_eta_jet_quark", Nbins,-5.0,5.0)
histdelPhi_darkquark_initial = ROOT.TH1F("dPhi_darkquark_initial", "dPhi_jet1_quark", Nbins, -3.15, 3.15)
#Pion number distributions
Nbins= 60
histpionnum = ROOT.TH1F("pion_number", "Pion number", Nbins, 0.0, 60.0)
histrhonum = ROOT.TH1F("rho_number", "Pion number", Nbins, 0.0, 60.0)
histstablepionnum = ROOT.TH1F("stable_pion_number", "Stable pion number", Nbins, 0.0, 60.0)
histunstablepionnum = ROOT.TH1F("unstable_pion_number", "Unstable pion number", Nbins, 0.0, 60.0)
#GeV distributions
Nbins = 150
histMET = ROOT.TH1F("jet_met", "Missing transverse energy", Nbins, 0.0, 1500.0)
histpionPT = ROOT.TH1F("pion_pt", "Pion pT", Nbins, 0.0, 1000.0)
histrhoPT = ROOT.TH1F("rho_pt", "Pion pT", Nbins, 0.0, 1000.0)
histquarkPT = ROOT.TH1F("quark_pt", "Quark pT", Nbins, 0.0, 500.0)
histdarkquarkPT_initial = ROOT.TH1F("initial_dark_quark_pt", "Quark pT", Nbins, 0.0, 1000.0)
histdarkquarkPT = ROOT.TH1F("dark_quark_pt", "Quark pT", Nbins, 0.0, 500.0)
hist_initialdarkquark_invmass = ROOT.TH1F("initial_dark_quark_invariantmass", "Quark invariant mass", Nbins, 0.0, 1500.0)
hist_darkquark_invmass= ROOT.TH1F("dark_quark_invariantmass", "Quark invariant mass", Nbins, 0.0, 1000.0)
histjetmJJ = ROOT.TH1F("jet_mJJ", " Invariant mass", Nbins, 0.0, 2500.0)

######2D HISTOGRAM######
###two_D_histogram= ROOT.TH2F("2dhistogram", "2D histogram of pion number and pT distribution",350,0,2000,50,0,50)
#### x and y - REMEMBER####

####COUNTERS,INDEXES AND IDS STORED HERE######
#index of Particle ID's for stable pions
piPID = [4900111, 4900121, 4900131, 4900141, 4900151, 4900161, 4900171, 4900181,
        4900211, 4900221, 4900231, 4900241, 4900251, 4900261, 4900271, 4900281,
        4900311, 4900321, 4900331, 4900341, 4900351, 4900361, 4900371, 4900381,
        4900411, 4900421, 4900431, 4900441, 4900451, 4900461, 4900471, 4900481,
        4900511, 4900521, 4900531, 4900541, 4900551, 4900561, 4900571, 4900581,
        4900611, 4900621, 4900631, 4900641, 4900651, 4900661, 4900671, 4900681,
        4900711, 4900721, 4900731, 4900741, 4900751, 4900761, 4900771, 4900781,
        4900811, 4900821, 4900831, 4900841, 4900851, 4900861, 4900871, 4900881
        ]
rhoPID = [4900113, 4900123, 4900133, 4900143, 4900153, 4900163, 4900173, 4900183,
        4900213, 4900223, 4900233, 4900243, 4900253, 4900263, 4900273, 4900283,
        4900313, 4900323, 4900333, 4900343, 4900353, 4900363, 4900373, 4900383,
        4900413, 4900423, 4900433, 4900443, 4900453, 4900463, 4900473, 4900483,
        4900513, 4900523, 4900533, 4900543, 4900553, 4900563, 4900573, 4900583,
        4900613, 4900623, 4900633, 4900643, 4900653, 4900663, 4900673, 4900683,
        4900713, 4900723, 4900733, 4900743, 4900753, 4900763, 4900773, 4900783,
        4900813, 4900823, 4900833, 4900843, 4900853, 4900863, 4900873, 4900883
        ]
#index of Particle ID's for on-diagonal stable pions (not necessary yet)
diagonalpiPID = [4900111, 4900221, 4900331, 4900441, 4900551, 4900661, 4900771, 4900881]
#index of number of jets per event
jetindex={}
#counts the number of events
eventcounter=1
#Dictionary for the PT j-th jet histogram 
PTjetJ={}
#minimum number of jets (to avoid overcounting)
minimum_jet_number=0
#All histograms stored in here
histlist = ROOT.TList()

#######
#######
#######

# Loop over all events
print(numberOfEntries)
for entry in range(0, numberOfEntries):
    # Load selected branches with data from specified event
    treeReader.ReadEntry(entry)
    # If event contains at least 2 jets
    if branchJet.GetEntries() > 1:
      ###### Defining jet-1, jet-2 vectors ######
        vec1 = GetJetVector(branchJet.At(0).PT , branchJet.At(0).Eta , branchJet.At(0).Phi , branchJet.At(0).Mass)
        vec2 = GetJetVector(branchJet.At(1).PT , branchJet.At(1).Eta , branchJet.At(1).Phi , branchJet.At(1).Mass)

      #######PRODUCING N-JET DISTRIBUTIONS######
        jetindex[eventcounter]=branchJet.GetEntries()
        maximum_jet_number=max(jetindex.values())
        for n in range(minimum_jet_number, maximum_jet_number):
            PTjetJ[n+1] = ROOT.TH1F("jet"+str(n+1), "", Nbins, 0.0, 2000.0)
        minimum_jet_number=maximum_jet_number
        for n in range(0, branchJet.GetEntries()):
            PTjetJ[n+1].Fill((branchJet.At(n)).PT)

      ######PRODUCING D_ETA AND INVARIANT MASS DISTRIBUTIONS######
        delEta = GetdEta(branchJet.At(0).Eta , branchJet.At(1).Eta)
        delPhi = GetdPhi(branchJet.At(0).Phi, branchJet.At(1).Phi)
        histjetmJJ.Fill((vec1+vec2).M())
        histdelEta.Fill(delEta)
        if delEta >= 0:
          histdelEta_positive.Fill(delEta)
        else:
          histdelEta_negative.Fill((-1)*delEta)
        histdelPhi.Fill(delPhi)
        if delPhi >= 0:
          histdelPhi_positive.Fill(delPhi)
        else:
          histdelPhi_negative.Fill(-1*delPhi)
      
      ######PRODUCING MET DISTRIBUTION######
        met = branchMET.At(0)
        m = met.MET
        METPhi = met.Phi
        METx = m* np.cos(METPhi)
        METy = m* np.sin(METPhi)
        vecmet = ROOT.TLorentzVector()
        vecmet.SetPxPyPzE(METx , METy , 0 , m)
        histMET.Fill(branchMET.At(0).MET)
        dPhi1 = GetdPhi(branchJet.At(0).Phi, met.Phi)
        dPhi2 = GetdPhi(branchJet.At(1).Phi , met.Phi)
        
        #counts the total number of pions per event
        pionnum=0
        #counts the total number of rhos per event
        rhonum=0
        #counts the number of stable pions per event
        st_pionnum=0
        #counts the number of unstable pions per event
        un_pionnum=0
        #counts the number of standard model quarks per event
        quarknum=0
        #counts the number of initial dark quarks per event
        darkquarknum_initial=0
        #counts the number of initial dark quarks per event
        darkquarknum=0

        ######PRODUCING PARTON DISTRIBUTIONS######

        for iptcl in range(0, branchPtcl.GetEntries()):
          if branchPtcl.At(iptcl).PID == 4900023 and branchPtcl.At(iptcl).Status == 22: 
            ###INITIAL DARK QUARK DISTRIBUTIONS###
            D1 = branchPtcl.At(iptcl).D1
            D2 = branchPtcl.At(iptcl).D2
            momptcl = D1
            while abs(branchPtcl.At(momptcl).PID) == 4900023:
              D1 = branchPtcl.At(momptcl).D1
              D2 = branchPtcl.At(momptcl).D2
              momptcl = D1
            vec1 = ROOT.TLorentzVector()
            vec2 = ROOT.TLorentzVector()
            vec1 = GetJetVector(branchPtcl.At(D1).PT, branchPtcl.At(D1).Eta, branchPtcl.At(D1).Phi, branchPtcl.At(D1).Mass)
            vec2 = GetJetVector(branchPtcl.At(D2).PT, branchPtcl.At(D2).Eta, branchPtcl.At(D2).Phi, branchPtcl.At(D2).Mass)
            histdarkquarkPT_initial.Fill(branchPtcl.At(D1).PT)
            histdarkquarkPT_initial.Fill(branchPtcl.At(D2).PT)
            delEta_darkquark_initial = GetdEta(branchPtcl.At(D1).Eta , branchPtcl.At(D2).Eta)
            delPhi_darkquark_initial = GetdPhi(branchPtcl.At(D1).Phi, branchPtcl.At(D2).Phi)

            histdelEta_darkquark_initial.Fill(delEta_darkquark_initial)
            if delEta_darkquark_initial >= 0:
              histdelEta_positive_darkquark_initial.Fill(delEta_darkquark_initial)
            else:
              histdelEta_negative_darkquark_initial.Fill((-1)*delEta_darkquark_initial)
            histdelPhi_darkquark_initial.Fill(delPhi_darkquark_initial)
            if delPhi_darkquark_initial >= 0:
              histdelPhi_positive_darkquark_initial.Fill(delPhi_darkquark_initial)
            else:
              histdelPhi_negative_darkquark_initial.Fill((-1)*delPhi_darkquark_initial)
            darkquarknum_initial=darkquarknum_initial+2
            hist_initialdarkquark_invmass.Fill((vec1+vec2).M())

          if abs(branchPtcl.At(iptcl).PID) in piPID:
            pionnum=pionnum+1
            histpionPT.Fill(branchPtcl.At(iptcl).PT)

            ###DARK QUARK DISTRIBUTIONS###
            M1 = branchPtcl.At(iptcl).M1
            M2 = branchPtcl.At(iptcl).M2
            momptcl = M1
            while abs(branchPtcl.At(momptcl).PID) == 4900113 or abs(branchPtcl.At(momptcl).PID) == 4900213:
              M1 = branchPtcl.At(momptcl).M1
              M2 = branchPtcl.At(momptcl).M2
              momptcl = M1
            #print(branchPtcl.At(D1).PID, branchPtcl.At(D2).PID)
            vec1 = ROOT.TLorentzVector()
            vec2 = ROOT.TLorentzVector()
            vec1 = GetJetVector(branchPtcl.At(M1).PT, branchPtcl.At(M1).Eta, branchPtcl.At(M1).Phi, branchPtcl.At(M1).Mass)
            vec2 = GetJetVector(branchPtcl.At(M2).PT, branchPtcl.At(M2).Eta, branchPtcl.At(M2).Phi, branchPtcl.At(M2).Mass)
            histdarkquarkPT.Fill(branchPtcl.At(M1).PT)
            histdarkquarkPT.Fill(branchPtcl.At(M2).PT)
            delEta_darkquark = GetdEta(branchPtcl.At(M1).Eta , branchPtcl.At(M2).Eta)
            delPhi_darkquark = GetdPhi(branchPtcl.At(M1).Phi, branchPtcl.At(M2).Phi)

            histdelEta_darkquark.Fill(delEta_darkquark)
            if delEta_darkquark >= 0:
              histdelEta_positive_darkquark.Fill(delEta_darkquark)
            else:
              histdelEta_negative_darkquark.Fill((-1)*delEta_darkquark)
            histdelPhi_darkquark.Fill(delPhi_darkquark)
            if delPhi_darkquark >= 0:
              histdelPhi_positive_darkquark.Fill(delPhi_darkquark)
            else:
              histdelPhi_negative_darkquark.Fill((-1)*delPhi_darkquark)
            darkquarknum=darkquarknum+2
            hist_darkquark_invmass.Fill((vec1+vec2).M())

            ##############################################################################################################################


            ###STANDARD MODEL QUARK DISTRIBUTIONS
            if abs(branchPtcl.At(iptcl).PID) in diagonalpiPID and branchPtcl.At(iptcl).Status != 1:
              un_pionnum=un_pionnum+1
              D1 = branchPtcl.At(iptcl).D1
              D2 = branchPtcl.At(iptcl).D2
              vec1 = ROOT.TLorentzVector()
              vec2 = ROOT.TLorentzVector()
              vec1 = GetJetVector(branchPtcl.At(D1).PT, branchPtcl.At(D1).Eta, branchPtcl.At(D1).Phi, branchPtcl.At(D1).Mass)
              vec2 = GetJetVector(branchPtcl.At(D2).PT, branchPtcl.At(D2).Eta, branchPtcl.At(D2).Phi, branchPtcl.At(D2).Mass)
              if branchPtcl.At(D1).Status == 91:
                histquarkPT.Fill(branchPtcl.At(D1).PT)
                histquarkPT.Fill(branchPtcl.At(D2).PT)
                delEta_quark = GetdEta(branchPtcl.At(D1).Eta , branchPtcl.At(D2).Eta)
                delPhi_quark = GetdPhi(branchPtcl.At(D1).Phi, branchPtcl.At(D2).Phi)
                histdelEta_quark.Fill(delEta_quark)
                if delEta_quark >= 0:
                  histdelEta_positive_quark.Fill(delEta_quark)
                else:
                  histdelEta_negative_quark.Fill((-1)*delEta_quark)
                histdelPhi_quark.Fill(delPhi_quark)
                if delPhi_quark >= 0:
                  histdelPhi_positive_quark.Fill(delPhi_quark)
                else:
                  histdelPhi_negative_quark.Fill((-1)*delPhi_quark)
                quarknum=quarknum+2
            else:
              st_pionnum=st_pionnum+1

          if abs(branchPtcl.At(iptcl).PID) in rhoPID:
            rhonum=rhonum+1
            histrhoPT.Fill(branchPtcl.At(iptcl).PT)

        histpionnum.Fill(pionnum)
        histrhonum.Fill(pionnum)
        histunstablepionnum.Fill(un_pionnum)
        histstablepionnum.Fill(st_pionnum)
        histquarknum.Fill(quarknum)
        histdarkquarknum_initial.Fill(darkquarknum_initial)
        histdarkquarknum.Fill(darkquarknum)
        hdphi1.Fill(dPhi1)
        hdphi2.Fill(dPhi2)

        eventcounter = eventcounter + 1

### Normalising an adding all jet histograms in a loop ###
for n in range(0, minimum_jet_number):
              if PTjetJ[n+1].GetSumw2N()==0: PTjetJ[n+1].Sumw2(True)
              PTjetJ[n+1].Scale(1./PTjetJ[n+1].Integral())
              histlist.Add(PTjetJ[n+1])  
print("Accepted events: ", eventcounter-1)
print("Maximum number of jets: ", maximum_jet_number)


#######NORMALISE HISTOGRAMS
if histrhonum.GetSumw2N()==0: histrhonum.Sumw2(True)
if histrhoPT.GetSumw2N()==0: histrhoPT.Sumw2(True)
if histMET.GetSumw2N()==0: histMET.Sumw2(True)
if histjetmJJ.GetSumw2N()==0: histjetmJJ.Sumw2(True)
if histdelEta.GetSumw2N()==0: histdelEta.Sumw2(True)
if histpionnum.GetSumw2N()==0: histpionnum.Sumw2(True)
if histunstablepionnum.GetSumw2N()==0: histunstablepionnum.Sumw2(True)
if histstablepionnum.GetSumw2N()==0: histstablepionnum.Sumw2(True)
if histquarknum.GetSumw2N()==0: histquarknum.Sumw2(True)
if histquarkPT.GetSumw2N()==0: histquarkPT.Sumw2(True)
if histpionPT.GetSumw2N()==0: histpionPT.Sumw2(True)
if histdelPhi.GetSumw2N()==0: histdelPhi.Sumw2(True)
if hdphi1.GetSumw2N()==0: hdphi1.Sumw2(True)
if hdphi2.GetSumw2N()==0: hdphi2.Sumw2(True)
if histdelEta_positive.GetSumw2N()==0: histdelEta_positive.Sumw2(True)
if histdelEta_negative.GetSumw2N()==0: histdelEta_negative.Sumw2(True)
if histdelPhi_positive.GetSumw2N()==0: histdelPhi_positive.Sumw2(True)
if histdelPhi_negative.GetSumw2N()==0: histdelPhi_negative.Sumw2(True)
if histdelEta_positive_quark.GetSumw2N()==0: histdelEta_positive_quark.Sumw2(True)
if histdelEta_negative_quark.GetSumw2N()==0: histdelEta_negative_quark.Sumw2(True)
if histdelPhi_positive_quark.GetSumw2N()==0: histdelPhi_positive_quark.Sumw2(True)
if histdelPhi_negative_quark.GetSumw2N()==0: histdelPhi_negative_quark.Sumw2(True)
if histdelEta_positive_darkquark.GetSumw2N()==0: histdelEta_positive_darkquark.Sumw2(True)
if histdelEta_negative_darkquark.GetSumw2N()==0: histdelEta_negative_darkquark.Sumw2(True)
if histdelPhi_positive_darkquark.GetSumw2N()==0: histdelPhi_positive_darkquark.Sumw2(True)
if histdelPhi_negative_darkquark.GetSumw2N()==0: histdelPhi_negative_darkquark.Sumw2(True)
if histdelPhi_quark.GetSumw2N()==0: histdelPhi_quark.Sumw2(True)
if histdelEta_quark.GetSumw2N()==0: histdelEta_quark.Sumw2(True)
if histdelEta_darkquark.GetSumw2N()==0: histdelEta_darkquark.Sumw2(True)
if histdelPhi_darkquark.GetSumw2N()==0: histdelPhi_darkquark.Sumw2(True)
if histdarkquarknum.GetSumw2N()==0: histdarkquarknum.Sumw2(True)
if histdarkquarkPT.GetSumw2N()==0: histdarkquarkPT.Sumw2(True)
if histdelEta_positive_darkquark_initial.GetSumw2N()==0: histdelEta_positive_darkquark_initial.Sumw2(True)
if histdelEta_negative_darkquark_initial.GetSumw2N()==0: histdelEta_negative_darkquark_initial.Sumw2(True)
if histdelPhi_positive_darkquark_initial.GetSumw2N()==0: histdelPhi_positive_darkquark_initial.Sumw2(True)
if histdelPhi_negative_darkquark_initial.GetSumw2N()==0: histdelPhi_negative_darkquark_initial.Sumw2(True)
if histdarkquarknum_initial.GetSumw2N()==0: histdarkquarknum_initial.Sumw2(True)
if histdarkquarkPT_initial.GetSumw2N()==0: histdarkquarkPT_initial.Sumw2(True)
if histdelEta_darkquark_initial.GetSumw2N()==0: histdelEta_darkquark_initial.Sumw2(True)
if histdelPhi_darkquark_initial.GetSumw2N()==0: histdelPhi_darkquark_initial.Sumw2(True)
if hist_initialdarkquark_invmass.GetSumw2N()==0: hist_initialdarkquark_invmass.Sumw2(True)
if hist_darkquark_invmass.GetSumw2N()==0: hist_darkquark_invmass.Sumw2(True)

######SCALE HISTOGRAMS######
histrhonum.Scale(1./histrhonum.Integral())
histrhoPT.Scale(1./histrhoPT.Integral())
histMET.Scale(1./histMET.Integral())
histjetmJJ.Scale(1./histjetmJJ.Integral())
histdelEta.Scale(1./histdelEta.Integral())
histpionnum.Scale(1./histpionnum.Integral())
histunstablepionnum.Scale(1./histunstablepionnum.Integral())
histstablepionnum.Scale(1./histstablepionnum.Integral())
histquarknum.Scale(1./histquarknum.Integral())
histquarkPT.Scale(1./histquarkPT.Integral())
histpionPT.Scale(1./histpionPT.Integral())
histdelPhi.Scale(1./histdelPhi.Integral())
hdphi1.Scale(1./hdphi1.Integral())
hdphi2.Scale(1./hdphi2.Integral())
histdelEta_positive.Scale(1./histdelEta_positive.Integral())
histdelEta_negative.Scale(1./histdelEta_negative.Integral())
histdelPhi_positive.Scale(1./histdelPhi_positive.Integral())
histdelPhi_negative.Scale(1./histdelPhi_negative.Integral())
histdelEta_positive_quark.Scale(1./histdelEta_positive_quark.Integral())
histdelEta_negative_quark.Scale(1./histdelEta_negative_quark.Integral())
histdelPhi_positive_quark.Scale(1./histdelPhi_positive_quark.Integral())
histdelPhi_negative_quark.Scale(1./histdelPhi_negative_quark.Integral())
histdelEta_positive_darkquark.Scale(1./histdelEta_positive_darkquark.Integral())
histdelEta_negative_darkquark.Scale(1./histdelEta_negative_darkquark.Integral())
histdelPhi_positive_darkquark.Scale(1./histdelPhi_positive_darkquark.Integral())
histdelPhi_negative_darkquark.Scale(1./histdelPhi_negative_darkquark.Integral())
histdelEta_quark.Scale(1./histdelEta_quark.Integral())
histdelPhi_quark.Scale(1./histdelPhi_quark.Integral())
histdelEta_darkquark.Scale(1./histdelEta_darkquark.Integral())
histdelPhi_darkquark.Scale(1./histdelPhi_darkquark.Integral())
histdarkquarknum.Scale(1./histdarkquarknum.Integral())
histdarkquarkPT.Scale(1./histdarkquarkPT.Integral())
histdelEta_positive_darkquark_initial.Scale(1./histdelEta_positive_darkquark_initial.Integral())
histdelEta_negative_darkquark_initial.Scale(1./histdelEta_negative_darkquark_initial.Integral())
histdelPhi_positive_darkquark_initial.Scale(1./histdelPhi_positive_darkquark_initial.Integral())
histdelPhi_negative_darkquark_initial.Scale(1./histdelPhi_negative_darkquark_initial.Integral())
histdarkquarknum_initial.Scale(1./histdarkquarknum_initial.Integral())
histdarkquarkPT_initial.Scale(1./histdarkquarkPT_initial.Integral())
histdelEta_darkquark_initial.Scale(1./histdelEta_darkquark_initial.Integral())
histdelPhi_darkquark_initial.Scale(1./histdelPhi_darkquark_initial.Integral())
hist_initialdarkquark_invmass.Scale(1./hist_initialdarkquark_invmass.Integral())
hist_darkquark_invmass.Scale(1./hist_darkquark_invmass.Integral())

######ADD HISTOGRAMS TO LIST######
histlist.Add(histrhonum)
histlist.Add(histrhoPT)
histlist.Add(histMET)
histlist.Add(histjetmJJ)
histlist.Add(histdelEta)
histlist.Add(histpionnum)
histlist.Add(histunstablepionnum)
histlist.Add(histstablepionnum)
histlist.Add(histquarknum)
histlist.Add(histquarkPT)
histlist.Add(histpionPT)
histlist.Add(histdelPhi)
histlist.Add(hdphi1)
histlist.Add(hdphi2)
histlist.Add(histdelEta_positive)
histlist.Add(histdelEta_negative)
histlist.Add(histdelPhi_positive)
histlist.Add(histdelPhi_negative)
histlist.Add(histdelEta_positive_quark)
histlist.Add(histdelEta_negative_quark)
histlist.Add(histdelPhi_positive_quark)
histlist.Add(histdelPhi_negative_quark)
histlist.Add(histdelEta_positive_darkquark)
histlist.Add(histdelEta_negative_darkquark)
histlist.Add(histdelPhi_positive_darkquark)
histlist.Add(histdelPhi_negative_darkquark)
histlist.Add(histdelEta_quark)
histlist.Add(histdelPhi_quark)
histlist.Add(histdelEta_darkquark)
histlist.Add(histdelPhi_darkquark)
histlist.Add(histdarkquarknum)
histlist.Add(histdarkquarkPT)
histlist.Add(histdelEta_positive_darkquark_initial)
histlist.Add(histdelEta_negative_darkquark_initial)
histlist.Add(histdelPhi_positive_darkquark_initial)
histlist.Add(histdelPhi_negative_darkquark_initial)
histlist.Add(histdarkquarknum_initial)
histlist.Add(histdarkquarkPT_initial)
histlist.Add(histdelEta_darkquark_initial)
histlist.Add(histdelPhi_darkquark_initial)
histlist.Add(hist_initialdarkquark_invmass)
histlist.Add(hist_darkquark_invmass)

rootFile = ROOT.TFile(outputFile, "RECREATE")
histlist.Write()
rootFile.Close()





