#!/bin/bash
text2workspace.py ttH_hbb_13TeV_sl_j5_tge4_high.txt
combine ttH_hbb_13TeV_sl_j5_tge4_high.root -M MaxLikelihoodFit --rMin -5 --rMax 5  
combine ttH_hbb_13TeV_sl_j5_tge4_high.root -M MaxLikelihoodFit --toysFrequentist –t 5000 --rMin -5 --rMax 5 --expectSignal 0 –n toys
