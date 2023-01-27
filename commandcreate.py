
from calendar import c
import os, sys
from tkinter import E
import numpy as np
from numpy import *
from pylab import *
from scipy.interpolate import interp1d
import ROOT
from array import array
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-a' , help = "file-address")
parser.add_argument('-b' , help = "pion_decay", type = float)
parser.add_argument('-c' , help = "lam", type = float)
parser.add_argument('-d' , help = "Nf", type= float)
parser.add_argument('-e' , help = "Nc")
parser.add_argument('-f' , help = "m_pi / lam", type = float)
parser.add_argument('-g' , help = "loop_order")

args = parser.parse_args()

file_address = args.a
pion_decay = args.b
lam = args.c
Nf = args.d
Nc = args.e
ml = args.f
order=args.g

pion_mass=ml*lam
print("Input files :")
print("File address: ", file_address)
print("Number of pion decays: ", pion_decay)
print("Lambda: ", lam)
print("Number of colours: ", Nc)
print("Number of flavours: ", Nf)
print("Pion lambda ratio: ", ml)
print("Loop order: ", order)
if pion_mass<3.2:
    print("Strange decay")
    decaynumber=3
elif pion_mass<10:
    print("Charm decay")
    decaynumber=4
else:
    print("Bottom decay")
    decaynumber=5

decay_product_dictionary={}

pion_mass=ml*lam
rho_mass=((5.76+1.5*(ml**2))**(1/2))*lam
constit_mass=(1+(1/5.5**(2))*(ml**2))*lam
r_inv=pion_decay/Nf

string = '%s/Nc%sNf%.0f_sFoff_pp_%.0fpi_decay_lam_%.0f_probvec_75_order_%s.cmnd'%(file_address,Nc,Nf,pion_decay,lam,order)
print(string)
with open(string, 'w') as f:
    f.write(
    '! Produced with commandcreate.py'
    '! The value of r_inv approx 1 - 1/N_f * Br(pi_0->visible) \n'
    '! This is because we have Nf number of pions\n'
    '! Out of these any number of pions can decay to the SM\n'
    '\n'
    '\n'
    '! 1) Settings that will be used in a main program.\n'
    'Main:numberOfEvents = 50000          ! number of events to generate\n'
    'Main:timesAllowErrors = 3          ! abort run after this many flawed events\n'
    '\n'
    '!Random seed = 0 gives random seed based on time.\n'
    'Random:setSeed = on\n'
    'Random:seed = 0\n'
    '\n'
    '! 2) Settings related to output in init(), next() and stat().\n'
    'Init:showChangedSettings = on      ! list changed settings\n'
    'Init:showAllSettings = on         ! list all settings\n'
    'Init:showChangedParticleData = on  ! list changed particle data\n'
    'Init:showAllParticleData = on     ! list all particle data\n'
    'Next:numberCount = 1000            ! print message every n events\n'
    'Next:numberShowLHA = 1             ! print LHA information n times\n'
    'Next:numberShowInfo = 1            ! print event information n times\n'
    'Next:numberShowProcess = 1         ! print process record n times\n'
    'Next:numberShowEvent = 1           ! print event record n times\n'
    'Stat:showPartonLevel = on          ! additional statistics on MPI\n'
    '\n'
    '! 3) Beam parameter settings. Values below agree with default ones.\n'
    'Beams:idA = 2212                   ! first beam, p = 2212, pbar = -2212\n'
    'Beams:idB = 2212                   ! second beam, p = 2212, pbar = -2212\n'
    'Beams:eCM = 13000.                 ! CM energy of collision\n'
    '\n'
    '! 4) Turn on the production process.\n'
    'HiddenValley:ffbar2Zv = on\n'
    '4900023:m0 = 1000\n'
    '4900023:mMin = 800\n'
    '4900023:mMax = 1200\n'
    '4900023:mWidth = 20.0\n'
    '4900023:doForceWidth = on\n'
    '\n'
    '4900023:oneChannel = 1 0.0750000 100 1 -1\n'
    '4900023:addChannel = 1 0.0750000 100 2 -2\n'
    '4900023:addChannel = 1 0.0750000 100 3 -3\n'
    '4900023:addChannel = 1 0.0750000 100 4 -4\n'
    '4900023:addChannel = 1 0.0750000 100 5 -5\n'
    '4900023:addChannel = 1 0.0750000 100 6 -6\n'
    '4900023:addChannel = 1 0.0250000 100 11 -11\n'
    '4900023:addChannel = 1 0.0250000 100 12 -12\n'
    '4900023:addChannel = 1 0.0250000 100 13 -13\n'
    '4900023:addChannel = 1 0.0250000 100 14 -14\n'
    '4900023:addChannel = 1 0.0250000 100 15 -15\n'
    '4900023:addChannel = 1 0.0250000 100 16 -16\n'
    '4900023:addChannel = 1 0.4 100 4900101 -4900101\n'
    '\n'
    '! 5) For MC efficiency let zprime decay only to dark quarks.\n'
    '4900023:onMode = off\n'
    '4900023:onIfAny = 4900101\n'
    '\n'
    '\n'
    '! 6) Decouple default bifundamental quarks.\n'
    '4900001:m0 = 5000\n'
    '4900002:m0 = 5000\n'
    '4900003:m0 = 5000\n'
    '4900004:m0 = 5000\n'
    '4900005:m0 = 5000\n'
    '4900006:m0 = 5000\n'
    '4900011:m0 = 5000\n'
    '4900012:m0 = 5000\n'
    '4900013:m0 = 5000\n'
    '4900014:m0 = 5000\n'
    '4900015:m0 = 5000\n'
    '4900016:m0 = 5000\n'
    '\n'
    '\n'
    '! 7) Set mass of dark quark.\n'
    '! this is Lambda + current quark mass\n'
    f'4900101:m0 = {round(constit_mass, 3)}\n'
    '\n'
    '! 8) Set gauge group and number of flavours\n'
    '! caution: user should not use Nc=2 or Nf=1 without consulting theory experts.\n'
    f'HiddenValley:Ngauge  = {Nc}\n'
    f'HiddenValley:Nflav  = {round(Nf, 0)}\n'
    'HiddenValley:separateFlav  = off\n'
    '\n'
    '! 9) Set other darksector parmaters, must setHiddenValley:alphaOrder = 1.0\n'
    f'HiddenValley:Lambda = {lam}\n'
    'HiddenValley:FSR = on\n'
    'HiddenValley:fragment = on\n'
    f'HiddenValley:alphaOrder = {order}\n'
    'HiddenValley:probVector = 0.75\n'
    '\n'
    '! 10) Set masses and branching ratios of spin-1 and spin-0 mesons\n'
    '! For now this treats flavor singlets on the same footing as flavor adjoint multiplets, this is not entirely consistent.\n'
    '! adjustment to flavor singlet is possible but only with separateFlav=on\n'
    f'4900111:m0 = {round(pion_mass, 3)}\n'
    f'4900211:m0 = {round(pion_mass, 3)}\n'
    f'4900113:m0 = {round(rho_mass, 3)}\n'
    f'4900213:m0 = {round(rho_mass, 3)}\n'
    '\n'
    '! Declare additional stable state to emulate invisible fraction\n'
    '53:m0 = 0.0\n'
    '53:isResonance = false\n'
    '53:mayDecay = off\n'
    '\n'
    '! All spin-1 mesons decay to spin-0 mesons\n'
    '\n'
    '4900113:mayDecay = on\n'
    '4900113:oneChannel = 1 1 100 4900211 -4900211\n'
    '\n'
    '4900213:mayDecay = on\n'
    '4900213:oneChannel = 1 1 100 4900111 4900211\n'
    '\n'
    '\n'
    '! All off-diagonal pions are stable and all diagonal pions decay\n'
    '! Being pseudo-scalar, the pion should decay almost exclusively to heaviest available state\n'
    '4900211:mayDecay = off\n')
    if pion_decay == 0:
        f.write(
            '4900111:mayDecay = off\n')
    else:
        f.write(
        '4900111:mayDecay = on\n'
        f'4900111:onechannel = 1  {round(r_inv, 3)}  91 -{decaynumber} {decaynumber}\n'
        f'4900111:addchannel = 1  {round(1-r_inv, 3)}  100 53 -53\n')
