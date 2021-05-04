#!/bin/bash

hadd -f histograms_pelee_overlay_wgt.root ../Run1/histograms_pelee_overlay_wgt.root ../Run2/histograms_pelee_overlay_wgt.root ../Run3/histograms_pelee_overlay_wgt.root 
hadd -f histograms_pelee_dirt_wgt.root ../Run1/histograms_pelee_dirt_wgt.root ../Run2/histograms_pelee_dirt_wgt.root ../Run3/histograms_pelee_dirt_wgt.root 
hadd -f histograms_pelee_bnb.root ../Run1/histograms_pelee_bnb.root ../Run2/histograms_pelee_bnb.root ../Run3/histograms_pelee_bnb.root 
hadd -f histograms_pelee_ext.root ../Run1/histograms_pelee_ext.root ../Run2/histograms_pelee_ext.root ../Run3/histograms_pelee_ext.root
