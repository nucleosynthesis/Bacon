#!/usr/bin/env python
# make the mono-X workspaces for statistical analysis
import sys, array
import ROOT

if len(sys.argv)<2: sys.exit("Not enough arguments, pass configuration file")

# usual stuff of which config to import :)
config = sys.argv[1]
x = __import__(config)

# Make an output file 
fout = ROOT.TFile(x.outFileName,"RECREATE")
fin  = ROOT.TFile.Open(x.inFileName)

# load the library to get a model builder
#ROOT.gROOT.ProcessLine(".L ./libModelBuilder.so")
ROOT.gROOT.ProcessLine(".L ./ModelBuilder.cc+") # euch, need to create a makefile to avoid recompile

for catid,cat in enumerate(x.categories): # loop over categories
  fdir = fout.mkdir("category_%d"%catid)

  mb = ROOT.ModelBuilder(catid)
  mb.luminosity     = float(x.categories[catid]["luminosity"])
  mb.run_correction = int(x.categories[catid]["run_correction"])
  mb.run_bkgsys     = int(x.categories[catid]["run_bkgsys"])
  mb.cutstring      = str(x.categories[catid]["cutstring"])
  mb.varstring      = str(x.categories[catid]["varstring"])
  mb.xmin	    = float(x.categories[catid]["min"])
  mb.xmax	    = float(x.categories[catid]["max"])
  mb.runphotonmodel = int("runPhoton" in x.categories[catid].keys())
  try: 
    bins = x.categories[catid]["bins"]
    histobase = ROOT.TH1F("base_%d"%catid,"base"
        ,len(bins)-1
	,array.array('d',bins))

  except:
    histobase = ROOT.TH1F("base_%d"%catid,"base"
  	,x.categories[catid]["nbins"]
  	,x.categories[catid]["min"]
  	,x.categories[catid]["max"])

  mb.lMet = histobase.Clone()  # ???

  mb.buildModel(fin,fdir);
  cstr = ROOT.TNamed("cut_category_%d"%catid,x.categories[catid]["cutstring"])
  fdir.cd()
  cstr.Write()

print "done!, Model saved in -> ", fout.GetName()

