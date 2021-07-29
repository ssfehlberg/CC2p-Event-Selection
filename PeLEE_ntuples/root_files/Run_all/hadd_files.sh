#!/bin/bash

#hadds all the original files together
hadd -f histograms_pelee_overlay_wgt.root ../Run1/histograms_pelee_overlay_wgt.root ../Run2/histograms_pelee_overlay_wgt.root ../Run3/histograms_pelee_overlay_wgt.root 
hadd -f histograms_pelee_dirt_wgt.root ../Run1/histograms_pelee_dirt_wgt.root ../Run2/histograms_pelee_dirt_wgt.root ../Run3/histograms_pelee_dirt_wgt.root 
hadd -f histograms_pelee_bnb.root ../Run1/histograms_pelee_bnb.root ../Run2/histograms_pelee_bnb.root ../Run3/histograms_pelee_bnb.root 
hadd -f histograms_pelee_ext.root ../Run1/histograms_pelee_ext.root ../Run2/histograms_pelee_ext.root ../Run3/histograms_pelee_ext.root
hadd -f histograms_mc_eff.root ../Run1/histograms_mc_eff.root ../Run2/histograms_mc_eff.root ../Run3/histograms_mc_eff.root

#hadds all the xsec binning files together
hadd -f histograms_pelee_xsec_overlay_wgt.root ../Run1/histograms_pelee_xsec_overlay_wgt.root ../Run2/histograms_pelee_xsec_overlay_wgt.root ../Run3/histograms_pelee_xsec_overlay_wgt.root
hadd -f histograms_pelee_xsec_dirt_wgt.root ../Run1/histograms_pelee_xsec_dirt_wgt.root ../Run2/histograms_pelee_xsec_dirt_wgt.root ../Run3/histograms_pelee_xsec_dirt_wgt.root
hadd -f histograms_pelee_xsec_bnb.root ../Run1/histograms_pelee_xsec_bnb.root ../Run2/histograms_pelee_xsec_bnb.root ../Run3/histograms_pelee_xsec_bnb.root
hadd -f histograms_pelee_xsec_ext.root ../Run1/histograms_pelee_xsec_ext.root ../Run2/histograms_pelee_xsec_ext.root ../Run3/histograms_pelee_xsec_ext.root
