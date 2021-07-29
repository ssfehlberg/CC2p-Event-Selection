//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 18 15:10:23 2020 by ROOT version 6.12/06
// from TTree NeutrinoSelectionFilter/Neutrino Selection TTree
// found on file: /pnfs/uboone/persistent/users/davidc/searchingfornues/v08_00_00_43/0928/prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root
//////////////////////////////////////////////////////////

#ifndef twoproton_nuwro_overlay_h
#define twoproton_nuwro_overlay_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "string"
#include "map"
#include "vector"
#include "variables.h"
using namespace Constants;

class twoproton_nuwro_overlay {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           selected;
   Int_t           run;
   Int_t           sub;
   Int_t           evt;
   UInt_t          trk_id;
   UInt_t          shr_id;
   UInt_t          trk2_id;
   UInt_t          shr2_id;
   Float_t         shr_energy_tot;
   Float_t         shr_energy;
   Float_t         shr_energy_tot_cali;
   Float_t         shr_energy_cali;
   Float_t         shr_theta;
   Float_t         shr_phi;
   Float_t         shr_pca_0;
   Float_t         shr_pca_1;
   Float_t         shr_pca_2;
   Float_t         shr_px;
   Float_t         shr_py;
   Float_t         shr_pz;
   Float_t         shr_openangle;
   Float_t         shr_tkfit_start_x;
   Float_t         shr_tkfit_start_y;
   Float_t         shr_tkfit_start_z;
   Float_t         shr_tkfit_theta;
   Float_t         shr_tkfit_phi;
   Float_t         shr_start_x;
   Float_t         shr_start_y;
   Float_t         shr_start_z;
   Float_t         shr_dedx_Y;
   Float_t         shr_dedx_V;
   Float_t         shr_dedx_U;
   Float_t         shr_dedx_Y_cali;
   Float_t         shr_dedx_V_cali;
   Float_t         shr_dedx_U_cali;
   Float_t         shr_tkfit_dedx_Y;
   Float_t         shr_tkfit_dedx_V;
   Float_t         shr_tkfit_dedx_U;
   Float_t         shr_tkfit_dedx_max;
   UInt_t          shr_tkfit_nhits_Y;
   UInt_t          shr_tkfit_nhits_V;
   UInt_t          shr_tkfit_nhits_U;
   Float_t         shr_llrpid_dedx_Y;
   Float_t         shr_llrpid_dedx_V;
   Float_t         shr_llrpid_dedx_U;
   Float_t         shr_llrpid_dedx;
   Float_t         shr_tkfit_dedx_Y_alt;
   Float_t         shr_tkfit_dedx_V_alt;
   Float_t         shr_tkfit_dedx_U_alt;
   UInt_t          shr_tkfit_nhits_Y_alt;
   UInt_t          shr_tkfit_nhits_V_alt;
   UInt_t          shr_tkfit_nhits_U_alt;
   Float_t         trkfit;
   UInt_t          shr_tkfit_npoints;
   UInt_t          shr_tkfit_npointsvalid;
   Float_t         shr_trkfitmedangle;
   Float_t         shrmoliereavg;
   Float_t         shrmoliererms;
   Float_t         shr1shr2moliereavg;
   Float_t         shr1shr2moliererms;
   Float_t         shr1trk1moliereavg;
   Float_t         shr1trk1moliererms;
   Float_t         shr1trk2moliereavg;
   Float_t         shr1trk2moliererms;
   UChar_t         ismerged;
   Float_t         merge_bestdot;
   Float_t         merge_bestdist;
   Float_t         merge_vtx_x;
   Float_t         merge_vtx_y;
   Float_t         merge_vtx_z;
   UInt_t          merge_tk_ipfp;
   Float_t         shr_tkfit_2cm_dedx_Y;
   Float_t         shr_tkfit_2cm_dedx_V;
   Float_t         shr_tkfit_2cm_dedx_U;
   UInt_t          shr_tkfit_2cm_nhits_Y;
   UInt_t          shr_tkfit_2cm_nhits_V;
   UInt_t          shr_tkfit_2cm_nhits_U;
   Float_t         shr_tkfit_gap05_dedx_Y;
   Float_t         shr_tkfit_gap05_dedx_V;
   Float_t         shr_tkfit_gap05_dedx_U;
   UInt_t          shr_tkfit_gap05_nhits_Y;
   UInt_t          shr_tkfit_gap05_nhits_V;
   UInt_t          shr_tkfit_gap05_nhits_U;
   Float_t         shr_tkfit_gap10_dedx_Y;
   Float_t         shr_tkfit_gap10_dedx_V;
   Float_t         shr_tkfit_gap10_dedx_U;
   UInt_t          shr_tkfit_gap10_nhits_Y;
   UInt_t          shr_tkfit_gap10_nhits_V;
   UInt_t          shr_tkfit_gap10_nhits_U;
   Float_t         shr_chipr;
   Float_t         shr_chimu;
   Float_t         shr_bragg_p;
   Float_t         shr_bragg_mu;
   Float_t         shr_bragg_mip;
   Float_t         shr_bragg_kaon;
   Float_t         shr_bragg_pion;
   Float_t         tksh_distance;
   Float_t         tksh_angle;
   Float_t         shr_distance;
   Float_t         shr_score;
   Int_t           shr_bkt_pdg;
   Float_t         shr_bkt_purity;
   Float_t         shr_bkt_completeness;
   Float_t         shr_bkt_E;
   Float_t         trk_len;
   Float_t         trk_theta;
   Float_t         trk_phi;
   Float_t         trk_energy;
   Float_t         trk_energy_muon;
   Float_t         trk_energy_muon_mcs;
   Float_t         trk_energy_tot;
   Float_t         trk_energy_muon_tot;
   Float_t         trk_distance;
   Float_t         trk_score;
   Int_t           trk_bkt_pdg;
   Float_t         trk_bkt_purity;
   Float_t         trk_bkt_completeness;
   Float_t         trk_bkt_E;
   Float_t         trk_chipr_best;
   Float_t         trk_chipr_worst;
   Float_t         trk_chimu_best;
   Float_t         trk_chimu_worst;
   Float_t         trk_chipr;
   Float_t         trk_chimu;
   Float_t         trk_pida;
   Float_t         trk_bragg_p;
   Float_t         trk_bragg_mu;
   Float_t         trk_bragg_mip;
   Float_t         trk_bragg_kaon;
   Float_t         trk_bragg_pion;
   UInt_t          trk_hits_max;
   UInt_t          shr_hits_max;
   UInt_t          trk_hits_2nd;
   UInt_t          shr_hits_2nd;
   Float_t         trkshrhitdist0;
   Float_t         trkshrhitdist1;
   Float_t         trkshrhitdist2;
   Float_t         trk2shrhitdist0;
   Float_t         trk2shrhitdist1;
   Float_t         trk2shrhitdist2;
   Float_t         trk1trk2hitdist0;
   Float_t         trk1trk2hitdist1;
   Float_t         trk1trk2hitdist2;
   UInt_t          total_hits_y;
   Float_t         extra_energy_y;
   Float_t         trk_energy_hits_tot;
   UInt_t          subcluster;
   UInt_t          shrsubclusters0;
   UInt_t          shrsubclusters1;
   UInt_t          shrsubclusters2;
   Float_t         shrclusfrac0;
   Float_t         shrclusfrac1;
   Float_t         shrclusfrac2;
   Float_t         shrclusdir0;
   Float_t         shrclusdir1;
   Float_t         shrclusdir2;
   UInt_t          shr_hits_tot;
   UInt_t          shr_hits_y_tot;
   UInt_t          shr_hits_u_tot;
   UInt_t          shr_hits_v_tot;
   UInt_t          trk_hits_tot;
   UInt_t          trk_hits_y_tot;
   UInt_t          trk_hits_u_tot;
   UInt_t          trk_hits_v_tot;
   Float_t         _elecclusters_U_charge;
   Float_t         _elecclusters_V_charge;
   Float_t         _elecclusters_Y_charge;
   Int_t           _elecclusters_U_N;
   Int_t           _elecclusters_V_N;
   Int_t           _elecclusters_Y_N;
   UInt_t          n_tracks_contained;
   UInt_t          n_showers_contained;
   Float_t         matched_E;
   Float_t         hits_ratio;
   Float_t         contained_fraction;
   Float_t         sps_contained_fraction;
   Float_t         pt;
   Float_t         p;
   Float_t         pt_assume_muon;
   Float_t         p_assume_muon;
   Float_t         reco_e;
   Float_t         dvtx;
   Float_t         dtrk;
   Float_t         contained_sps_ratio;
   vector<double>  *dtrk_x_boundary;
   vector<double>  *dtrk_y_boundary;
   vector<double>  *dtrk_z_boundary;
   vector<double>  *dshr_x_boundary;
   vector<double>  *dshr_y_boundary;
   vector<double>  *dshr_z_boundary;
   vector<double>  *dvtx_x_boundary;
   vector<double>  *dvtx_y_boundary;
   vector<double>  *dvtx_z_boundary;
   vector<vector<double> > *dtrk_boundary;
   vector<vector<double> > *dvtx_boundary;
   vector<vector<double> > *dshr_boundary;
   vector<vector<double> > *dmc_boundary;
   Float_t         CosmicIP;
   Float_t         CosmicIPAll3D;
   Float_t         CosmicDirAll3D;
   Float_t         CosmicIPAll2DEnds;
   Float_t         CosmicDirAll2DEnds;
   Float_t         CosmicIPAll2DOvlp;
   Float_t         CosmicDirAll2DOvlp;
   Float_t         leeweight;
   Float_t         true_pt;
   Float_t         true_pt_visible;
   Float_t         true_p;
   Float_t         true_p_visible;
   Float_t         true_e_visible;
   Float_t         _opfilter_pe_beam;
   Float_t         _opfilter_pe_veto;
   Int_t           nu_pdg;
   Int_t           ccnc;
   Int_t           nu_parent_pdg;
   Int_t           nu_hadron_pdg;
   Int_t           nu_decay_mode;
   Int_t           interaction;
   Float_t         nu_e;
   Float_t         nu_l;
   Float_t         nu_pt;
   Float_t         theta;
   Bool_t          isVtxInFiducial;
   Bool_t          truthFiducial;
   Float_t         true_nu_vtx_t;
   Float_t         true_nu_vtx_x;
   Float_t         true_nu_vtx_y;
   Float_t         true_nu_vtx_z;
   Float_t         true_nu_vtx_sce_x;
   Float_t         true_nu_vtx_sce_y;
   Float_t         true_nu_vtx_sce_z;
   Float_t         reco_nu_vtx_x;
   Float_t         reco_nu_vtx_y;
   Float_t         reco_nu_vtx_z;
   Float_t         reco_nu_vtx_sce_x;
   Float_t         reco_nu_vtx_sce_y;
   Float_t         reco_nu_vtx_sce_z;
   Int_t           nmuon;
   Float_t         muon_e;
   Float_t         muon_c;
   Float_t         muon_p;
   Int_t           nelec;
   Float_t         elec_e;
   Float_t         elec_c;
   Float_t         elec_p;
   Float_t         elec_vx;
   Float_t         elec_vy;
   Float_t         elec_vz;
   Float_t         elec_px;
   Float_t         elec_py;
   Float_t         elec_pz;
   Int_t           npi0;
   Float_t         pi0_e;
   Float_t         pi0_c;
   Float_t         pi0_p;
   Int_t           nneutron;
   Int_t           nproton;
   Float_t         proton_e;
   Float_t         proton_c;
   Float_t         proton_p;
   Int_t           npion;
   Float_t         pion_e;
   Float_t         pion_c;
   Float_t         pion_p;
   Int_t           neta;
   Float_t         eta_e;
   Int_t           nslice;
   Int_t           crtveto;
   Float_t         crthitpe;
   vector<int>     *pfp_slice_idx;
   Int_t           category;
   vector<int>     *backtracked_pdg;
   vector<float>   *backtracked_e;
   vector<int>     *backtracked_tid;
   vector<float>   *backtracked_purity;
   vector<float>   *backtracked_completeness;
   vector<float>   *backtracked_overlay_purity;
   vector<float>   *backtracked_px;
   vector<float>   *backtracked_py;
   vector<float>   *backtracked_pz;
   vector<float>   *backtracked_start_x;
   vector<float>   *backtracked_start_y;
   vector<float>   *backtracked_start_z;
   vector<float>   *backtracked_start_t;
   vector<float>   *backtracked_start_U;
   vector<float>   *backtracked_start_V;
   vector<float>   *backtracked_start_Y;
   vector<float>   *backtracked_sce_start_x;
   vector<float>   *backtracked_sce_start_y;
   vector<float>   *backtracked_sce_start_z;
   vector<float>   *backtracked_sce_start_U;
   vector<float>   *backtracked_sce_start_V;
   vector<float>   *backtracked_sce_start_Y;
   Float_t         lep_e;
   Int_t           pass;
   Int_t           swtrig;
   Int_t           evnhits;
   Int_t           slpdg;
   Int_t           slnhits;
   Int_t           n_pfps;
   Int_t           n_tracks;
   Int_t           n_showers;
   vector<unsigned int> *pfp_generation_v;
   vector<unsigned int> *pfp_trk_daughters_v;
   vector<unsigned int> *pfp_shr_daughters_v;
   vector<float>   *trk_score_v;
   vector<int>     *pfpdg;
   vector<int>     *pfnhits;
   vector<int>     *pfnplanehits_U;
   vector<int>     *pfnplanehits_V;
   vector<int>     *pfnplanehits_Y;
   vector<int>     *pfpplanesubclusters_U;
   vector<int>     *pfpplanesubclusters_V;
   vector<int>     *pfpplanesubclusters_Y;
   vector<float>   *pfpplanesubhitfracmax_U;
   vector<float>   *pfpplanesubhitfracmax_V;
   vector<float>   *pfpplanesubhitfracmax_Y;
   UInt_t          hits_u;
   UInt_t          hits_v;
   UInt_t          hits_y;
   Float_t         topological_score;
   Float_t         slclustfrac;
   vector<int>     *mc_pdg;
   vector<float>   *mc_E;
   vector<float>   *mc_vx;
   vector<float>   *mc_vy;
   vector<float>   *mc_vz;
   vector<float>   *mc_endx;
   vector<float>   *mc_endy;
   vector<float>   *mc_endz;
   vector<float>   *mc_px;
   vector<float>   *mc_py;
   vector<float>   *mc_pz;
   vector<float>   *mc_completeness;
   vector<float>   *mc_purity;
   string          *endmuonprocess;
   Float_t         endmuonmichel;
   Int_t           filter_antibdt;
   Int_t           filter_ncpi0;
   Int_t           filter_pi0;
   Int_t           filter_ccinclusive;
   map<string,vector<double> > *weights;
   vector<unsigned short> *weightsFlux;
   vector<unsigned short> *weightsGenie;
   vector<unsigned short> *weightsReint;
   Float_t         weightSpline;
   Float_t         weightTune;
   Float_t         weightSplineTimesTune;
   Double_t        knobRPAup;
   Double_t        knobRPAdn;
   Double_t        knobCCMECup;
   Double_t        knobCCMECdn;
   Double_t        knobAxFFCCQEup;
   Double_t        knobAxFFCCQEdn;
   Double_t        knobVecFFCCQEup;
   Double_t        knobVecFFCCQEdn;
   Double_t        knobDecayAngMECup;
   Double_t        knobDecayAngMECdn;
   Double_t        knobThetaDelta2Npiup;
   Double_t        knobThetaDelta2Npidn;
   Float_t         flash_pe;
   Float_t         flash_time;
   Float_t         nu_flashmatch_score;
   Float_t         best_cosmic_flashmatch_score;
   Float_t         best_obviouscosmic_flashmatch_score;
   vector<float>   *cosmic_flashmatch_score_v;
   Float_t         mcf_nu_e;
   Float_t         mcf_lep_e;
   Int_t           mcf_actvol;
   Int_t           mcf_nmm;
   Int_t           mcf_nmp;
   Int_t           mcf_nem;
   Int_t           mcf_nep;
   Int_t           mcf_np0;
   Int_t           mcf_npp;
   Int_t           mcf_npm;
   Float_t         mcf_mcshr_elec_etot;
   Int_t           mcf_pass_ccpi0;
   Int_t           mcf_pass_ncpi0;
   Int_t           mcf_pass_ccnopi;
   Int_t           mcf_pass_ncnopi;
   Int_t           mcf_pass_cccpi;
   Int_t           mcf_pass_nccpi;
   vector<float>   *X_SpcPts_v;
   vector<float>   *Y_SpcPts_v;
   vector<float>   *Z_SpcPts_v;
   UInt_t          shr_id_MCStool;
   UInt_t          shr_hits_max_MCStool;
   UInt_t          n_showers_contained_MCStool;
   vector<float>   *trkshrscore_v;
   Float_t         shrPCA_1Cr;
   Float_t         shrPCA_2Cr;
   Float_t         shrPCA_3Cr;
   Float_t         shrPCA_1Ce;
   Float_t         shrPCA_2Ce;
   Float_t         shrPCA_3Ce;
   Float_t         shrPCA1CAS;
   Float_t         shrPCA2CAS;
   Float_t         shrPCA3CAS;
   Float_t         shrPCA_1Cr2h;
   Float_t         shrPCA_2Cr2h;
   Float_t         shrPCA_3Cr2h;
   Float_t         shrPCA_1Cr1h;
   Float_t         shrPCA_2Cr1h;
   Float_t         shrPCA_3Cr1h;
   Float_t         shrMCSMom;
   Float_t         shrMCSMom1h;
   Float_t         shrMCSMom2h;
   Float_t         shrPCALen;
   UInt_t          n_shrSpcPts;
   vector<float>   *PCAWin_1Cr_5cm;
   vector<float>   *PCAWin_2Cr_5cm;
   vector<float>   *PCAWin_3Cr_5cm;
   vector<float>   *PCAWin_dist_5cm;
   vector<int>     *PCAWin_npts_5cm;
   Float_t         shrStart_5cm;
   Float_t         shrStartMCS_5cm;
   Float_t         shrMCSAS_5cm;
   Float_t         shrPCA1CAS_5cm;
   Float_t         shrPCA2CAS_5cm;
   Float_t         shrPCA3CAS_5cm;
   Float_t         shrPCA1CMed_5cm;
   vector<float>   *PCAWin_1Cr_2_5cm;
   vector<float>   *PCAWin_2Cr_2_5cm;
   vector<float>   *PCAWin_3Cr_2_5cm;
   vector<float>   *PCAWin_dist_2_5cm;
   vector<int>     *PCAWin_npts_2_5cm;
   Float_t         shrStart_2_5cm;
   Float_t         shrStartMCS_2_5cm;
   Float_t         shrMCSAS_2_5cm;
   Float_t         shrPCA1CAS_2_5cm;
   Float_t         shrPCA2CAS_2_5cm;
   Float_t         shrPCA3CAS_2_5cm;
   Float_t         shrPCA1CMed_2_5cm;
   Float_t         DeltaMed;
   Float_t         DeltaMed1h;
   Float_t         DeltaMed2h;
   Float_t         DeltaRMS;
   Float_t         DeltaRMS1h;
   Float_t         DeltaRMS2h;
   Float_t         CylFrac_1cm;
   Float_t         CylFrac1h_1cm;
   Float_t         CylFrac2h_1cm;
   Float_t         CylFrac_2cm;
   Float_t         CylFrac1h_2cm;
   Float_t         CylFrac2h_2cm;
   Float_t         CylFrac_3cm;
   Float_t         CylFrac1h_3cm;
   Float_t         CylFrac2h_3cm;
   Float_t         CylFrac_4cm;
   Float_t         CylFrac1h_4cm;
   Float_t         CylFrac2h_4cm;
   Float_t         CylFrac_5cm;
   Float_t         CylFrac1h_5cm;
   Float_t         CylFrac2h_5cm;
   Float_t         NeutrinoEnergy0;
   Float_t         NeutrinoEnergy1;
   Float_t         NeutrinoEnergy2;
   Float_t         SliceCaloEnergy0;
   Float_t         SliceCaloEnergy1;
   Float_t         SliceCaloEnergy2;
   Float_t         pi0_mcgamma0_e;
   Float_t         pi0_mcgamma0_px;
   Float_t         pi0_mcgamma0_py;
   Float_t         pi0_mcgamma0_pz;
   Float_t         pi0_mcrcdot0;
   Float_t         pi0_mcrce0;
   Float_t         pi0_mcgamma1_e;
   Float_t         pi0_mcgamma1_px;
   Float_t         pi0_mcgamma1_py;
   Float_t         pi0_mcgamma1_pz;
   Float_t         pi0_mcrcdot1;
   Float_t         pi0_mcrce1;
   Int_t           pi0_nshower;
   Int_t           pi0_ntrack;
   Int_t           pi0_ngamma;
   Float_t         pi0_radlen1;
   Float_t         pi0_radlen2;
   Float_t         pi0_dot1;
   Float_t         pi0_dot2;
   Float_t         pi0_energy1_Y;
   Float_t         pi0_energy2_Y;
   Float_t         pi0_dir1_x;
   Float_t         pi0_dir1_y;
   Float_t         pi0_dir1_z;
   Float_t         pi0_dir2_x;
   Float_t         pi0_dir2_y;
   Float_t         pi0_dir2_z;
   Float_t         pi0_dedx1_Y;
   Float_t         pi0_dedx2_Y;
   Float_t         pi0_dedx1_fit_Y;
   Float_t         pi0_dedx2_fit_Y;
   Float_t         pi0_energy1_V;
   Float_t         pi0_energy2_V;
   Float_t         pi0_dedx1_V;
   Float_t         pi0_dedx2_V;
   Float_t         pi0_dedx1_fit_V;
   Float_t         pi0_dedx2_fit_V;
   Float_t         pi0_energy1_U;
   Float_t         pi0_energy2_U;
   Float_t         pi0_dedx1_U;
   Float_t         pi0_dedx2_U;
   Float_t         pi0_dedx1_fit_U;
   Float_t         pi0_dedx2_fit_U;
   Float_t         pi0_shrscore1;
   Float_t         pi0_shrscore2;
   Float_t         pi0_gammadot;
   Float_t         pi0_mass_Y;
   Float_t         pi0_mass_V;
   Float_t         pi0_mass_U;
   Float_t         pi0_rc_vtx_x;
   Float_t         pi0_rc_vtx_y;
   Float_t         pi0_rc_vtx_z;
   Int_t           pi0truth_gamma_parent;
   Float_t         pi0truth_elec_edep;
   Float_t         pi0truth_elec_etot;
   Float_t         pi0truth_elec_dist;
   Int_t           pi0truth_elec_parent;
   Int_t           pi0truth_gamma1_tid;
   Float_t         pi0truth_gamma1_edep;
   Float_t         pi0truth_gamma1_etot;
   Float_t         pi0truth_gamma1_dist;
   Float_t         pi0truth_gamma1_elec1;
   Float_t         pi0truth_gamma1_elec2;
   Float_t         pi0truth_gamma1_xpos;
   Float_t         pi0truth_gamma1_ypos;
   Float_t         pi0truth_gamma1_zpos;
   Int_t           pi0truth_gamma2_tid;
   Float_t         pi0truth_gamma2_edep;
   Float_t         pi0truth_gamma2_etot;
   Float_t         pi0truth_gamma2_dist;
   Float_t         pi0truth_gamma2_elec1;
   Float_t         pi0truth_gamma2_elec2;
   Float_t         pi0truth_gamma2_xpos;
   Float_t         pi0truth_gamma2_ypos;
   Float_t         pi0truth_gamma2_zpos;
   Float_t         pi0truth_gammadot;
   Int_t           pi0truth_run;
   Int_t           pi0truth_sub;
   Int_t           pi0truth_evt;
   Int_t           nflag_pl1;
   Int_t           nnoise_pl1;
   Int_t           nslhits_pl1;
   Int_t           nslnoise_pl1;
   Int_t           nhits_pl1;
   Float_t         frac_slnoise_pl1;
   Float_t         secondshower_U_charge;
   Int_t           secondshower_U_nhit;
   Float_t         secondshower_U_vtxdist;
   Float_t         secondshower_U_eigenratio;
   Float_t         secondshower_U_dot;
   Float_t         secondshower_U_dir;
   Float_t         secondshower_V_charge;
   Int_t           secondshower_V_nhit;
   Float_t         secondshower_V_vtxdist;
   Float_t         secondshower_V_eigenratio;
   Float_t         secondshower_V_dot;
   Float_t         secondshower_V_dir;
   Float_t         secondshower_Y_charge;
   Int_t           secondshower_Y_nhit;
   Float_t         secondshower_Y_vtxdist;
   Float_t         secondshower_Y_eigenratio;
   Float_t         secondshower_Y_dot;
   Float_t         secondshower_Y_dir;
   vector<float>   *shr_dedx_u_v;
   vector<float>   *shr_dedx_v_v;
   vector<float>   *shr_dedx_y_v;
   vector<float>   *shr_energy_u_v;
   vector<float>   *shr_energy_v_v;
   vector<float>   *shr_energy_y_v;
   vector<unsigned long> *shr_pfp_id_v;
   vector<float>   *shr_start_x_v;
   vector<float>   *shr_start_y_v;
   vector<float>   *shr_start_z_v;
   vector<float>   *shr_dist_v;
   vector<float>   *shr_start_U_v;
   vector<float>   *shr_start_V_v;
   vector<float>   *shr_px_v;
   vector<float>   *shr_py_v;
   vector<float>   *shr_pz_v;
   vector<float>   *shr_openangle_v;
   vector<float>   *shr_theta_v;
   vector<float>   *shr_phi_v;
   vector<float>   *shr_pitch_u_v;
   vector<float>   *shr_pitch_v_v;
   vector<float>   *shr_pitch_y_v;
   vector<int>     *shr_tkfit_nhits_v;
   vector<float>   *shr_tkfit_start_x_v;
   vector<float>   *shr_tkfit_start_y_v;
   vector<float>   *shr_tkfit_start_z_v;
   vector<float>   *shr_tkfit_start_U_v;
   vector<float>   *shr_tkfit_start_V_v;
   vector<float>   *shr_tkfit_theta_v;
   vector<float>   *shr_tkfit_phi_v;
   vector<float>   *shr_tkfit_pitch_u_v;
   vector<float>   *shr_tkfit_pitch_v_v;
   vector<float>   *shr_tkfit_pitch_y_v;
   vector<float>   *shr_tkfit_dedx_u_v;
   vector<float>   *shr_tkfit_dedx_v_v;
   vector<float>   *shr_tkfit_dedx_y_v;
   vector<float>   *shr_tkfit_gap10_dedx_u_v;
   vector<float>   *shr_tkfit_gap10_dedx_v_v;
   vector<float>   *shr_tkfit_gap10_dedx_y_v;
   vector<int>     *shr_tkfit_dedx_nhits_u_v;
   vector<int>     *shr_tkfit_dedx_nhits_v_v;
   vector<int>     *shr_tkfit_dedx_nhits_y_v;
   vector<float>   *shr_llr_pid_u_v;
   vector<float>   *shr_llr_pid_v_v;
   vector<float>   *shr_llr_pid_y_v;
   vector<float>   *shr_llr_pid_v;
   vector<float>   *shr_llr_pid_score_v;
   vector<float>   *shr_moliere_avg_v;
   vector<float>   *shr_moliere_rms_v;
   Int_t           evnunhits;
   Int_t           evlepnhits;
   Int_t           evpronhits;
   Int_t           evpi1nhits;
   Int_t           evpi0nhits;
   Int_t           evneunhits;
   Int_t           evgamnhits;
   Int_t           evothnhits;
   Int_t           slnunhits;
   Int_t           sllepnhits;
   Int_t           slpronhits;
   Int_t           slpi1nhits;
   Int_t           slpi0nhits;
   Int_t           slneunhits;
   Int_t           slgamnhits;
   Int_t           slothnhits;
   vector<int>     *pfnunhits;
   vector<int>     *pflepnhits;
   vector<int>     *pfpronhits;
   vector<int>     *pfpi1nhits;
   vector<int>     *pfpi0nhits;
   vector<int>     *pfneunhits;
   vector<int>     *pfgamnhits;
   vector<int>     *pfothnhits;
   Float_t         nu_completeness_from_pfp;
   Float_t         nu_purity_from_pfp;
   vector<float>   *trk_bragg_p_v;
   vector<float>   *trk_bragg_mu_v;
   vector<float>   *trk_bragg_mip_v;
   vector<float>   *trk_pida_v;
   vector<float>   *trk_pid_chipr_v;
   vector<float>   *trk_pid_chipi_v;
   vector<float>   *trk_pid_chika_v;
   vector<float>   *trk_pid_chimu_v;
   vector<float>   *trk_bragg_p_u_v;
   vector<float>   *trk_bragg_mu_u_v;
   vector<float>   *trk_bragg_mip_u_v;
   vector<float>   *trk_pida_u_v;
   vector<float>   *trk_pid_chipr_u_v;
   vector<float>   *trk_pid_chipi_u_v;
   vector<float>   *trk_pid_chika_u_v;
   vector<float>   *trk_pid_chimu_u_v;
   vector<float>   *trk_bragg_p_v_v;
   vector<float>   *trk_bragg_mu_v_v;
   vector<float>   *trk_bragg_mip_v_v;
   vector<float>   *trk_pida_v_v;
   vector<float>   *trk_pid_chipr_v_v;
   vector<float>   *trk_pid_chipi_v_v;
   vector<float>   *trk_pid_chika_v_v;
   vector<float>   *trk_pid_chimu_v_v;
   vector<unsigned long> *trk_pfp_id_v;
   vector<float>   *trk_dir_x_v;
   vector<float>   *trk_dir_y_v;
   vector<float>   *trk_dir_z_v;
   vector<float>   *trk_start_x_v;
   vector<float>   *trk_start_y_v;
   vector<float>   *trk_start_z_v;
   vector<float>   *trk_sce_start_x_v;
   vector<float>   *trk_sce_start_y_v;
   vector<float>   *trk_sce_start_z_v;
   vector<float>   *trk_end_x_v;
   vector<float>   *trk_end_y_v;
   vector<float>   *trk_end_z_v;
   vector<float>   *trk_sce_end_x_v;
   vector<float>   *trk_sce_end_y_v;
   vector<float>   *trk_sce_end_z_v;
   vector<float>   *trk_distance_v;
   vector<float>   *trk_theta_v;
   vector<float>   *trk_phi_v;
   vector<float>   *trk_len_v;
   vector<float>   *trk_mcs_muon_mom_v;
   vector<float>   *trk_range_muon_mom_v;
   vector<float>   *trk_energy_proton_v;
   vector<float>   *trk_energy_muon_v;
   vector<float>   *trk_calo_energy_u_v;
   vector<float>   *trk_calo_energy_v_v;
   vector<float>   *trk_calo_energy_y_v;
   vector<float>   *trk_llr_pid_u_v;
   vector<float>   *trk_llr_pid_v_v;
   vector<float>   *trk_llr_pid_y_v;
   vector<float>   *trk_llr_pid_v;
   vector<float>   *trk_llr_pid_score_v;
   Float_t         bdt_nuNCpi0;
   Float_t         bdt_numuCCpi0;
   Float_t         bdt_numuCC;
   Float_t         bdt_ext;
   Float_t         bdt_cosmic;
   Float_t         bdt_global;
   Int_t           pass_antibdt_filter;
   Float_t         bdt_pi0_np;
   Float_t         bdt_nonpi0_np;
   Float_t         bdt_bkg_0p;
   Float_t         anglediff_Y;
   Float_t         anglediff_V;
   Float_t         anglediff_U;
   Float_t         trkpid;

   // List of branches
   TBranch        *b_selected;   //!
   TBranch        *b_run;   //!
   TBranch        *b_sub;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_trk_pfp_id;   //!
   TBranch        *b_shr_pfp_id;   //!
   TBranch        *b_trk2_pfp_id;   //!
   TBranch        *b_shr2_pfp_id;   //!
   TBranch        *b_shr_energy_tot;   //!
   TBranch        *b_shr_energy;   //!
   TBranch        *b_shr_energy_tot_cali;   //!
   TBranch        *b_shr_energy_cali;   //!
   TBranch        *b_shr_theta;   //!
   TBranch        *b_shr_phi;   //!
   TBranch        *b_shr_pca_0;   //!
   TBranch        *b_shr_pca_1;   //!
   TBranch        *b_shr_pca_2;   //!
   TBranch        *b_shr_px;   //!
   TBranch        *b_shr_py;   //!
   TBranch        *b_shr_pz;   //!
   TBranch        *b_shr_openangle;   //!
   TBranch        *b_shr_tkfit_start_x;   //!
   TBranch        *b_shr_tkfit_start_y;   //!
   TBranch        *b_shr_tkfit_start_z;   //!
   TBranch        *b_shr_tkfit_theta;   //!
   TBranch        *b_shr_tkfit_phi;   //!
   TBranch        *b_shr_start_x;   //!
   TBranch        *b_shr_start_y;   //!
   TBranch        *b_shr_start_z;   //!
   TBranch        *b_shr_dedx_Y;   //!
   TBranch        *b_shr_dedx_V;   //!
   TBranch        *b_shr_dedx_U;   //!
   TBranch        *b_shr_dedx_Y_cali;   //!
   TBranch        *b_shr_dedx_V_cali;   //!
   TBranch        *b_shr_dedx_U_cali;   //!
   TBranch        *b_shr_tkfit_dedx_Y;   //!
   TBranch        *b_shr_tkfit_dedx_V;   //!
   TBranch        *b_shr_tkfit_dedx_U;   //!
   TBranch        *b_shr_tkfit_dedx_max;   //!
   TBranch        *b_shr_tkfit_nhits_Y;   //!
   TBranch        *b_shr_tkfit_nhits_V;   //!
   TBranch        *b_shr_tkfit_nhits_U;   //!
   TBranch        *b_shr_llrpid_dedx_Y;   //!
   TBranch        *b_shr_llrpid_dedx_V;   //!
   TBranch        *b_shr_llrpid_dedx_U;   //!
   TBranch        *b_shr_llrpid_dedx;   //!
   TBranch        *b_shr_tkfit_dedx_Y_alt;   //!
   TBranch        *b_shr_tkfit_dedx_V_alt;   //!
   TBranch        *b_shr_tkfit_dedx_U_alt;   //!
   TBranch        *b_shr_tkfit_nhits_Y_alt;   //!
   TBranch        *b_shr_tkfit_nhits_V_alt;   //!
   TBranch        *b_shr_tkfit_nhits_U_alt;   //!
   TBranch        *b__trkfit;   //!
   TBranch        *b_shr_tkfit_npoints;   //!
   TBranch        *b_shr_tkfit_npointsvalid;   //!
   TBranch        *b_shr_trkfitmedangle;   //!
   TBranch        *b_shrmoliereavg;   //!
   TBranch        *b_shrmoliererms;   //!
   TBranch        *b_shr1shr2moliereavg;   //!
   TBranch        *b_shr1shr2moliererms;   //!
   TBranch        *b_shr1trk1moliereavg;   //!
   TBranch        *b_shr1trk1moliererms;   //!
   TBranch        *b_shr1trk2moliereavg;   //!
   TBranch        *b_shr1trk2moliererms;   //!
   TBranch        *b_ismerged;   //!
   TBranch        *b_merge_bestdot;   //!
   TBranch        *b_merge_bestdist;   //!
   TBranch        *b_merge_vtx_x;   //!
   TBranch        *b_merge_vtx_y;   //!
   TBranch        *b_merge_vtx_z;   //!
   TBranch        *b_merge_tk_ipfp;   //!
   TBranch        *b_shr_tkfit_2cm_dedx_Y;   //!
   TBranch        *b_shr_tkfit_2cm_dedx_V;   //!
   TBranch        *b_shr_tkfit_2cm_dedx_U;   //!
   TBranch        *b_shr_tkfit_2cm_nhits_Y;   //!
   TBranch        *b_shr_tkfit_2cm_nhits_V;   //!
   TBranch        *b_shr_tkfit_2cm_nhits_U;   //!
   TBranch        *b_shr_tkfit_gap05_dedx_Y;   //!
   TBranch        *b_shr_tkfit_gap05_dedx_V;   //!
   TBranch        *b_shr_tkfit_gap05_dedx_U;   //!
   TBranch        *b_shr_tkfit_gap05_nhits_Y;   //!
   TBranch        *b_shr_tkfit_gap05_nhits_V;   //!
   TBranch        *b_shr_tkfit_gap05_nhits_U;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_Y;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_V;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_U;   //!
   TBranch        *b_shr_tkfit_gap10_nhits_Y;   //!
   TBranch        *b_shr_tkfit_gap10_nhits_V;   //!
   TBranch        *b_shr_tkfit_gap10_nhits_U;   //!
   TBranch        *b_shr_chipr;   //!
   TBranch        *b_shr_chimu;   //!
   TBranch        *b_shr_bragg_p;   //!
   TBranch        *b_shr_bragg_mu;   //!
   TBranch        *b_shr_bragg_mip;   //!
   TBranch        *b_shr_bragg_kaon;   //!
   TBranch        *b_shr_bragg_pion;   //!
   TBranch        *b_tksh_distance;   //!
   TBranch        *b_tksh_angle;   //!
   TBranch        *b_shr_distance;   //!
   TBranch        *b_shr_score;   //!
   TBranch        *b_shr_bkt_pdg;   //!
   TBranch        *b_shr_bkt_purity;   //!
   TBranch        *b_shr_bkt_completeness;   //!
   TBranch        *b_shr_bkt_E;   //!
   TBranch        *b_trk_len;   //!
   TBranch        *b_trk_theta;   //!
   TBranch        *b_trk_phi;   //!
   TBranch        *b_trk_energy;   //!
   TBranch        *b_trk_energy_muon;   //!
   TBranch        *b_trk_energy_muon_mcs;   //!
   TBranch        *b_trk_energy_tot;   //!
   TBranch        *b_trk_energy_muon_tot;   //!
   TBranch        *b_trk_distance;   //!
   TBranch        *b_trk_score;   //!
   TBranch        *b_trk_bkt_pdg;   //!
   TBranch        *b_trk_bkt_purity;   //!
   TBranch        *b_trk_bkt_completeness;   //!
   TBranch        *b_trk_bkt_E;   //!
   TBranch        *b_trk_chipr_best;   //!
   TBranch        *b_trk_chipr_worst;   //!
   TBranch        *b_trk_chimu_best;   //!
   TBranch        *b_trk_chimu_worst;   //!
   TBranch        *b_trk_chipr;   //!
   TBranch        *b_trk_chimu;   //!
   TBranch        *b_trk_pida;   //!
   TBranch        *b_trk_bragg_p;   //!
   TBranch        *b_trk_bragg_mu;   //!
   TBranch        *b_trk_bragg_mip;   //!
   TBranch        *b_trk_bragg_kaon;   //!
   TBranch        *b_trk_bragg_pion;   //!
   TBranch        *b_trk_hits_max;   //!
   TBranch        *b_shr_hits_max;   //!
   TBranch        *b_trk_hits_2nd;   //!
   TBranch        *b_shr_hits_2nd;   //!
   TBranch        *b_trkshrhitdist0;   //!
   TBranch        *b_trkshrhitdist1;   //!
   TBranch        *b_trkshrhitdist2;   //!
   TBranch        *b_trk2shrhitdist0;   //!
   TBranch        *b_trk2shrhitdist1;   //!
   TBranch        *b_trk2shrhitdist2;   //!
   TBranch        *b_trk1trk2hitdist0;   //!
   TBranch        *b_trk1trk2hitdist1;   //!
   TBranch        *b_trk1trk2hitdist2;   //!
   TBranch        *b_total_hits_y;   //!
   TBranch        *b_extra_energy_y;   //!
   TBranch        *b_trk_energy_hits_tot;   //!
   TBranch        *b_subcluster;   //!
   TBranch        *b_shrsubclusters0;   //!
   TBranch        *b_shrsubclusters1;   //!
   TBranch        *b_shrsubclusters2;   //!
   TBranch        *b_shrclusfrac0;   //!
   TBranch        *b_shrclusfrac1;   //!
   TBranch        *b_shrclusfrac2;   //!
   TBranch        *b_shrclusdir0;   //!
   TBranch        *b_shrclusdir1;   //!
   TBranch        *b_shrclusdir2;   //!
   TBranch        *b_shr_hits_tot;   //!
   TBranch        *b_shr_hits_y_tot;   //!
   TBranch        *b_shr_hits_u_tot;   //!
   TBranch        *b_shr_hits_v_tot;   //!
   TBranch        *b_trk_hits_tot;   //!
   TBranch        *b_trk_hits_y_tot;   //!
   TBranch        *b_trk_hits_u_tot;   //!
   TBranch        *b_trk_hits_v_tot;   //!
   TBranch        *b_elecclusters_U_charge;   //!
   TBranch        *b_elecclusters_V_charge;   //!
   TBranch        *b_elecclusters_Y_charge;   //!
   TBranch        *b_elecclusters_U_N;   //!
   TBranch        *b_elecclusters_V_N;   //!
   TBranch        *b_elecclusters_Y_N;   //!
   TBranch        *b_n_tracks_contained;   //!
   TBranch        *b_n_showers_contained;   //!
   TBranch        *b_matched_E;   //!
   TBranch        *b_hits_ratio;   //!
   TBranch        *b_contained_fraction;   //!
   TBranch        *b_sps_contained_fraction;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_p;   //!
   TBranch        *b_pt_assume_muon;   //!
   TBranch        *b_p_assume_muon;   //!
   TBranch        *b_reco_e;   //!
   TBranch        *b_dvtx;   //!
   TBranch        *b_dtrk;   //!
   TBranch        *b_contained_sps_ratio;   //!
   TBranch        *b_dtrk_x_boundary;   //!
   TBranch        *b_dtrk_y_boundary;   //!
   TBranch        *b_dtrk_z_boundary;   //!
   TBranch        *b_dshr_x_boundary;   //!
   TBranch        *b_dshr_y_boundary;   //!
   TBranch        *b_dshr_z_boundary;   //!
   TBranch        *b_dvtx_x_boundary;   //!
   TBranch        *b_dvtx_y_boundary;   //!
   TBranch        *b_dvtx_z_boundary;   //!
   TBranch        *b_dtrk_boundary;   //!
   TBranch        *b_dvtx_boundary;   //!
   TBranch        *b_dshr_boundary;   //!
   TBranch        *b_dmc_boundary;   //!
   TBranch        *b_CosmicIP;   //!
   TBranch        *b_CosmicIPAll3D;   //!
   TBranch        *b_CosmicDirAll3D;   //!
   TBranch        *b_CosmicIPAll2DEnds;   //!
   TBranch        *b_CosmicDirAll2DEnds;   //!
   TBranch        *b_CosmicIPAll2DOvlp;   //!
   TBranch        *b_CosmicDirAll2DOvlp;   //!
   TBranch        *b_leeweight;   //!
   TBranch        *b_true_pt;   //!
   TBranch        *b_true_pt_visible;   //!
   TBranch        *b_true_p;   //!
   TBranch        *b_true_p_visible;   //!
   TBranch        *b_true_e_visible;   //!
   TBranch        *b_opfilter_pe_beam;   //!
   TBranch        *b_opfilter_pe_veto;   //!
   TBranch        *b_nu_pdg;   //!
   TBranch        *b_ccnc;   //!
   TBranch        *b_nu_parent_pdg;   //!
   TBranch        *b_nu_hadron_pdg;   //!
   TBranch        *b_nu_decay_mode;   //!
   TBranch        *b_interaction;   //!
   TBranch        *b_nu_e;   //!
   TBranch        *b_nu_l;   //!
   TBranch        *b_nu_pt;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_isVtxInFiducial;   //!
   TBranch        *b_truthFiducial;   //!
   TBranch        *b_true_nu_vtx_t;   //!
   TBranch        *b_true_nu_vtx_x;   //!
   TBranch        *b_true_nu_vtx_y;   //!
   TBranch        *b_true_nu_vtx_z;   //!
   TBranch        *b_true_nu_vtx_sce_x;   //!
   TBranch        *b_true_nu_vtx_sce_y;   //!
   TBranch        *b_true_nu_vtx_sce_z;   //!
   TBranch        *b_reco_nu_vtx_x;   //!
   TBranch        *b_reco_nu_vtx_y;   //!
   TBranch        *b_reco_nu_vtx_z;   //!
   TBranch        *b_reco_nu_vtx_sce_x;   //!
   TBranch        *b_reco_nu_vtx_sce_y;   //!
   TBranch        *b_reco_nu_vtx_sce_z;   //!
   TBranch        *b_nmuon;   //!
   TBranch        *b_muon_e;   //!
   TBranch        *b_muon_c;   //!
   TBranch        *b_muon_p;   //!
   TBranch        *b_nelec;   //!
   TBranch        *b_elec_e;   //!
   TBranch        *b_elec_c;   //!
   TBranch        *b_elec_p;   //!
   TBranch        *b_elec_vx;   //!
   TBranch        *b_elec_vy;   //!
   TBranch        *b_elec_vz;   //!
   TBranch        *b_elec_px;   //!
   TBranch        *b_elec_py;   //!
   TBranch        *b_elec_pz;   //!
   TBranch        *b_npi0;   //!
   TBranch        *b_pi0_e;   //!
   TBranch        *b_pi0_c;   //!
   TBranch        *b_pi0_p;   //!
   TBranch        *b_nneutron;   //!
   TBranch        *b_nproton;   //!
   TBranch        *b_proton_e;   //!
   TBranch        *b_proton_c;   //!
   TBranch        *b_proton_p;   //!
   TBranch        *b_npion;   //!
   TBranch        *b_pion_e;   //!
   TBranch        *b_pion_c;   //!
   TBranch        *b_pion_p;   //!
   TBranch        *b_neta;   //!
   TBranch        *b_eta_e;   //!
   TBranch        *b_nslice;   //!
   TBranch        *b_crtveto;   //!
   TBranch        *b_crthitpe;   //!
   TBranch        *b_pfp_slice_idx;   //!
   TBranch        *b_category;   //!
   TBranch        *b_backtracked_pdg;   //!
   TBranch        *b_backtracked_e;   //!
   TBranch        *b_backtracked_tid;   //!
   TBranch        *b_backtracked_purity;   //!
   TBranch        *b_backtracked_completeness;   //!
   TBranch        *b_backtracked_overlay_purity;   //!
   TBranch        *b_backtracked_px;   //!
   TBranch        *b_backtracked_py;   //!
   TBranch        *b_backtracked_pz;   //!
   TBranch        *b_backtracked_start_x;   //!
   TBranch        *b_backtracked_start_y;   //!
   TBranch        *b_backtracked_start_z;   //!
   TBranch        *b_backtracked_start_t;   //!
   TBranch        *b_backtracked_start_U;   //!
   TBranch        *b_backtracked_start_V;   //!
   TBranch        *b_backtracked_start_Y;   //!
   TBranch        *b_backtracked_sce_start_x;   //!
   TBranch        *b_backtracked_sce_start_y;   //!
   TBranch        *b_backtracked_sce_start_z;   //!
   TBranch        *b_backtracked_sce_start_U;   //!
   TBranch        *b_backtracked_sce_start_V;   //!
   TBranch        *b_backtracked_sce_start_Y;   //!
   TBranch        *b_lep_e;   //!
   TBranch        *b_pass;   //!
   TBranch        *b_swtrig;   //!
   TBranch        *b_evnhits;   //!
   TBranch        *b_slpdg;   //!
   TBranch        *b_slnhits;   //!
   TBranch        *b_n_pfps;   //!
   TBranch        *b_n_tracks;   //!
   TBranch        *b_n_showers;   //!
   TBranch        *b_pfp_generation_v;   //!
   TBranch        *b_pfp_trk_daughters_v;   //!
   TBranch        *b_pfp_shr_daughters_v;   //!
   TBranch        *b_trk_score_v;   //!
   TBranch        *b_pfpdg;   //!
   TBranch        *b_pfnhits;   //!
   TBranch        *b_pfnplanehits_U;   //!
   TBranch        *b_pfnplanehits_V;   //!
   TBranch        *b_pfnplanehits_Y;   //!
   TBranch        *b_pfpplanesubclusters_U;   //!
   TBranch        *b_pfpplanesubclusters_V;   //!
   TBranch        *b_pfpplanesubclusters_Y;   //!
   TBranch        *b_pfpplanesubhitfracmax_U;   //!
   TBranch        *b_pfpplanesubhitfracmax_V;   //!
   TBranch        *b_pfpplanesubhitfracmax_Y;   //!
   TBranch        *b_hits_u;   //!
   TBranch        *b_hits_v;   //!
   TBranch        *b_hits_y;   //!
   TBranch        *b_topological_score;   //!
   TBranch        *b_slclustfrac;   //!
   TBranch        *b_mc_pdg;   //!
   TBranch        *b_mc_E;   //!
   TBranch        *b_mc_vx;   //!
   TBranch        *b_mc_vy;   //!
   TBranch        *b_mc_vz;   //!
   TBranch        *b_mc_endx;   //!
   TBranch        *b_mc_endy;   //!
   TBranch        *b_mc_endz;   //!
   TBranch        *b_mc_px;   //!
   TBranch        *b_mc_py;   //!
   TBranch        *b_mc_pz;   //!
   TBranch        *b_mc_completeness;   //!
   TBranch        *b_mc_purity;   //!
   TBranch        *b_endmuonprocess;   //!
   TBranch        *b_endmuonmichel;   //!
   TBranch        *b_filter_antibdt;   //!
   TBranch        *b_filter_ncpi0;   //!
   TBranch        *b_filter_pi0;   //!
   TBranch        *b_filter_ccinclusive;   //!
   TBranch        *b_weights;   //!
   TBranch        *b_weightsFlux;   //!
   TBranch        *b_weightsGenie;   //!
   TBranch        *b_weightsReint;   //!
   TBranch        *b_weightSpline;   //!
   TBranch        *b_weightTune;   //!
   TBranch        *b_weightSplineTimesTune;   //!
   TBranch        *b_knobRPAup;   //!
   TBranch        *b_knobRPAdn;   //!
   TBranch        *b_knobCCMECup;   //!
   TBranch        *b_knobCCMECdn;   //!
   TBranch        *b_knobAxFFCCQEup;   //!
   TBranch        *b_knobAxFFCCQEdn;   //!
   TBranch        *b_knobVecFFCCQEup;   //!
   TBranch        *b_knobVecFFCCQEdn;   //!
   TBranch        *b_knobDecayAngMECup;   //!
   TBranch        *b_knobDecayAngMECdn;   //!
   TBranch        *b_knobThetaDelta2Npiup;   //!
   TBranch        *b_knobThetaDelta2Npidn;   //!
   TBranch        *b_flash_pe;   //!
   TBranch        *b_flash_time;   //!
   TBranch        *b_nu_flashmatch_score;   //!
   TBranch        *b_best_cosmic_flashmatch_score;   //!
   TBranch        *b_best_obviouscosmic_flashmatch_score;   //!
   TBranch        *b_cosmic_flashmatch_score_v;   //!
   TBranch        *b_mcf_nu_e;   //!
   TBranch        *b_mcf_lep_e;   //!
   TBranch        *b_mcf_actvol;   //!
   TBranch        *b_mcf_nmm;   //!
   TBranch        *b_mcf_nmp;   //!
   TBranch        *b_mcf_nem;   //!
   TBranch        *b_mcf_nep;   //!
   TBranch        *b_mcf_np0;   //!
   TBranch        *b_mcf_npp;   //!
   TBranch        *b_mcf_npm;   //!
   TBranch        *b_mcf_mcshr_elec_etot;   //!
   TBranch        *b_mcf_pass_ccpi0;   //!
   TBranch        *b_mcf_pass_ncpi0;   //!
   TBranch        *b_mcf_pass_ccnopi;   //!
   TBranch        *b_mcf_pass_ncnopi;   //!
   TBranch        *b_mcf_pass_cccpi;   //!
   TBranch        *b_mcf_pass_nccpi;   //!
   TBranch        *b_X_SpcPts_v;   //!
   TBranch        *b_Y_SpcPts_v;   //!
   TBranch        *b_Z_SpcPts_v;   //!
   //TBranch        *b_shr_pfp_id;   //!
   TBranch        *b_shr_hits_max_MCStool;   //!
   TBranch        *b_n_showers_contained_MCStool;   //!
   TBranch        *b_trkshrscore_v;   //!
   TBranch        *b_shrPCA_1Cr;   //!
   TBranch        *b_shrPCA_2Cr;   //!
   TBranch        *b_shrPCA_3Cr;   //!
   TBranch        *b_shrPCA_1Ce;   //!
   TBranch        *b_shrPCA_2Ce;   //!
   TBranch        *b_shrPCA_3Ce;   //!
   TBranch        *b_shrPCA1CAS;   //!
   TBranch        *b_shrPCA2CAS;   //!
   TBranch        *b_shrPCA3CAS;   //!
   TBranch        *b_shrPCA_1Cr2h;   //!
   TBranch        *b_shrPCA_2Cr2h;   //!
   TBranch        *b_shrPCA_3Cr2h;   //!
   TBranch        *b_shrPCA_1Cr1h;   //!
   TBranch        *b_shrPCA_2Cr1h;   //!
   TBranch        *b_shrPCA_3Cr1h;   //!
   TBranch        *b_shrMCSMom;   //!
   TBranch        *b_shrMCSMom1h;   //!
   TBranch        *b_shrMCSMom2h;   //!
   TBranch        *b_shrPCALen;   //!
   TBranch        *b_n_shrSpcPts;   //!
   TBranch        *b_PCAWin_1Cr_5cm;   //!
   TBranch        *b_PCAWin_2Cr_5cm;   //!
   TBranch        *b_PCAWin_3Cr_5cm;   //!
   TBranch        *b_PCAWin_dist_5cm;   //!
   TBranch        *b_PCAWin_npts_5cm;   //!
   TBranch        *b_shrStart_5cm;   //!
   TBranch        *b_shrStartMCS_5cm;   //!
   TBranch        *b_shrMCSAS_5cm;   //!
   TBranch        *b_shrPCA1CAS_5cm;   //!
   TBranch        *b_shrPCA2CAS_5cm;   //!
   TBranch        *b_shrPCA3CAS_5cm;   //!
   TBranch        *b__shrPCA1CMed_5cm;   //!
   TBranch        *b_PCAWin_1Cr_2_5cm;   //!
   TBranch        *b_PCAWin_2Cr_2_5cm;   //!
   TBranch        *b_PCAWin_3Cr_2_5cm;   //!
   TBranch        *b_PCAWin_dist_2_5cm;   //!
   TBranch        *b_PCAWin_npts_2_5cm;   //!
   TBranch        *b_shrStart_2_5cm;   //!
   TBranch        *b_shrStartMCS_2_5cm;   //!
   TBranch        *b_shrMCSAS_2_5cm;   //!
   TBranch        *b_shrPCA1CAS_2_5cm;   //!
   TBranch        *b_shrPCA2CAS_2_5cm;   //!
   TBranch        *b_shrPCA3CAS_2_5cm;   //!
   TBranch        *b__shrPCA1CMed_2_5cm;   //!
   TBranch        *b_DeltaMed;   //!
   TBranch        *b_DeltaMed1h;   //!
   TBranch        *b_DeltaMed2h;   //!
   TBranch        *b_DeltaRMS;   //!
   TBranch        *b_DeltaRMS1h;   //!
   TBranch        *b_DeltaRMS2h;   //!
   TBranch        *b_CylFrac_1cm;   //!
   TBranch        *b_CylFrac1h_1cm;   //!
   TBranch        *b_CylFrac2h_1cm;   //!
   TBranch        *b_CylFrac_2cm;   //!
   TBranch        *b_CylFrac1h_2cm;   //!
   TBranch        *b_CylFrac2h_2cm;   //!
   TBranch        *b_CylFrac_3cm;   //!
   TBranch        *b_CylFrac1h_3cm;   //!
   TBranch        *b_CylFrac2h_3cm;   //!
   TBranch        *b_CylFrac_4cm;   //!
   TBranch        *b_CylFrac1h_4cm;   //!
   TBranch        *b_CylFrac2h_4cm;   //!
   TBranch        *b_CylFrac_5cm;   //!
   TBranch        *b_CylFrac1h_5cm;   //!
   TBranch        *b_CylFrac2h_5cm;   //!
   TBranch        *b_NeutrinoEnergy0;   //!
   TBranch        *b_NeutrinoEnergy1;   //!
   TBranch        *b_NeutrinoEnergy2;   //!
   TBranch        *b_SliceCaloEnergy0;   //!
   TBranch        *b_SliceCaloEnergy1;   //!
   TBranch        *b_SliceCaloEnergy2;   //!
   TBranch        *b_pi0_mcgamma0_e;   //!
   TBranch        *b_pi0_mcgamma0_px;   //!
   TBranch        *b_pi0_mcgamma0_py;   //!
   TBranch        *b_pi0_mcgamma0_pz;   //!
   TBranch        *b_pi0_mcrcdot0;   //!
   TBranch        *b_pi0_mcrce0;   //!
   TBranch        *b_pi0_mcgamma1_e;   //!
   TBranch        *b_pi0_mcgamma1_px;   //!
   TBranch        *b_pi0_mcgamma1_py;   //!
   TBranch        *b_pi0_mcgamma1_pz;   //!
   TBranch        *b_pi0_mcrcdot1;   //!
   TBranch        *b_pi0_mcrce1;   //!
   TBranch        *b_pi0_nshower;   //!
   TBranch        *b_pi0_ntrack;   //!
   TBranch        *b_pi0_ngamma;   //!
   TBranch        *b_pi0_radlen1;   //!
   TBranch        *b_pi0_radlen2;   //!
   TBranch        *b_pi0_dot1;   //!
   TBranch        *b_pi0_dot2;   //!
   TBranch        *b_pi0_energy1_Y;   //!
   TBranch        *b_pi0_energy2_Y;   //!
   TBranch        *b_pi0_dir1_x;   //!
   TBranch        *b_pi0_dir1_y;   //!
   TBranch        *b_pi0_dir1_z;   //!
   TBranch        *b_pi0_dir2_x;   //!
   TBranch        *b_pi0_dir2_y;   //!
   TBranch        *b_pi0_dir2_z;   //!
   TBranch        *b_pi0_dedx1_Y;   //!
   TBranch        *b_pi0_dedx2_Y;   //!
   TBranch        *b_pi0_dedx1_fit_Y;   //!
   TBranch        *b_pi0_dedx2_fit_Y;   //!
   TBranch        *b_pi0_energy1_V;   //!
   TBranch        *b_pi0_energy2_V;   //!
   TBranch        *b_pi0_dedx1_V;   //!
   TBranch        *b_pi0_dedx2_V;   //!
   TBranch        *b_pi0_dedx1_fit_V;   //!
   TBranch        *b_pi0_dedx2_fit_V;   //!
   TBranch        *b_pi0_energy1_U;   //!
   TBranch        *b_pi0_energy2_U;   //!
   TBranch        *b_pi0_dedx1_U;   //!
   TBranch        *b_pi0_dedx2_U;   //!
   TBranch        *b_pi0_dedx1_fit_U;   //!
   TBranch        *b_pi0_dedx2_fit_U;   //!
   TBranch        *b_pi0_shrscore1;   //!
   TBranch        *b_pi0_shrscore2;   //!
   TBranch        *b_pi0_gammadot;   //!
   TBranch        *b_pi0_mass_Y;   //!
   TBranch        *b_pi0_mass_V;   //!
   TBranch        *b_pi0_mass_U;   //!
   TBranch        *b_pi0_rc_vtx_x;   //!
   TBranch        *b_pi0_rc_vtx_y;   //!
   TBranch        *b_pi0_rc_vtx_z;   //!
   TBranch        *b_pi0truth_gamma_parent;   //!
   TBranch        *b_pi0truth_elec_edep;   //!
   TBranch        *b_pi0truth_elec_etot;   //!
   TBranch        *b_pi0truth_elec_dist;   //!
   TBranch        *b_pi0truth_elec_parent;   //!
   TBranch        *b_pi0truth_gamma1_tid;   //!
   TBranch        *b_pi0truth_gamma1_edep;   //!
   TBranch        *b_pi0truth_gamma1_etot;   //!
   TBranch        *b_pi0truth_gamma1_dist;   //!
   TBranch        *b_pi0truth_gamma1_elec1;   //!
   TBranch        *b_pi0truth_gamma1_elec2;   //!
   TBranch        *b_pi0truth_gamma1_xpos;   //!
   TBranch        *b_pi0truth_gamma1_ypos;   //!
   TBranch        *b_pi0truth_gamma1_zpos;   //!
   TBranch        *b_pi0truth_gamma2_tid;   //!
   TBranch        *b_pi0truth_gamma2_edep;   //!
   TBranch        *b_pi0truth_gamma2_etot;   //!
   TBranch        *b_pi0truth_gamma2_dist;   //!
   TBranch        *b_pi0truth_gamma2_elec1;   //!
   TBranch        *b_pi0truth_gamma2_elec2;   //!
   TBranch        *b_pi0truth_gamma2_xpos;   //!
   TBranch        *b_pi0truth_gamma2_ypos;   //!
   TBranch        *b_pi0truth_gamma2_zpos;   //!
   TBranch        *b_pi0truth_gammadot;   //!
   TBranch        *b_pi0truth_run;   //!
   TBranch        *b_pi0truth_sub;   //!
   TBranch        *b_pi0truth_evt;   //!
   TBranch        *b_nflag_pl1;   //!
   TBranch        *b_nnoise_pl1;   //!
   TBranch        *b_nslhits_pl1;   //!
   TBranch        *b_nslnoise_pl1;   //!
   TBranch        *b_nhits_pl1;   //!
   TBranch        *b_frac_slnoise_pl1;   //!
   TBranch        *b_secondshower_U_charge;   //!
   TBranch        *b_secondshower_U_nhit;   //!
   TBranch        *b_secondshower_U_vtxdist;   //!
   TBranch        *b_secondshower_U_eigenratio;   //!
   TBranch        *b_secondshower_U_dot;   //!
   TBranch        *b_secondshower_U_dir;   //!
   TBranch        *b_secondshower_V_charge;   //!
   TBranch        *b_secondshower_V_nhit;   //!
   TBranch        *b_secondshower_V_vtxdist;   //!
   TBranch        *b_secondshower_V_eigenratio;   //!
   TBranch        *b_secondshower_V_dot;   //!
   TBranch        *b_secondshower_V_dir;   //!
   TBranch        *b_secondshower_Y_charge;   //!
   TBranch        *b_secondshower_Y_nhit;   //!
   TBranch        *b_secondshower_Y_vtxdist;   //!
   TBranch        *b_secondshower_Y_eigenratio;   //!
   TBranch        *b_secondshower_Y_dot;   //!
   TBranch        *b_secondshower_Y_dir;   //!
   TBranch        *b_shr_dedx_u_v;   //!
   TBranch        *b_shr_dedx_v_v;   //!
   TBranch        *b_shr_dedx_y_v;   //!
   TBranch        *b_shr_energy_u_v;   //!
   TBranch        *b_shr_energy_v_v;   //!
   TBranch        *b_shr_energy_y_v;   //!
   TBranch        *b_shr_pfp_id_v;   //!
   TBranch        *b_shr_start_x_v;   //!
   TBranch        *b_shr_start_y_v;   //!
   TBranch        *b_shr_start_z_v;   //!
   TBranch        *b_shr_dist_v;   //!
   TBranch        *b_shr_start_U_v;   //!
   TBranch        *b_shr_start_V_v;   //!
   TBranch        *b_shr_px_v;   //!
   TBranch        *b_shr_py_v;   //!
   TBranch        *b_shr_pz_v;   //!
   TBranch        *b_shr_openangle_v;   //!
   TBranch        *b_shr_theta_v;   //!
   TBranch        *b_shr_phi_v;   //!
   TBranch        *b_shr_pitch_u_v;   //!
   TBranch        *b_shr_pitch_v_v;   //!
   TBranch        *b_shr_pitch_y_v;   //!
   TBranch        *b_shr_tkfit_nhits_v;   //!
   TBranch        *b_shr_tkfit_start_x_v;   //!
   TBranch        *b_shr_tkfit_start_y_v;   //!
   TBranch        *b_shr_tkfit_start_z_v;   //!
   TBranch        *b_shr_tkfit_start_U_v;   //!
   TBranch        *b_shr_tkfit_start_V_v;   //!
   TBranch        *b_shr_tkfit_theta_v;   //!
   TBranch        *b_shr_tkfit_phi_v;   //!
   TBranch        *b_shr_tkfit_pitch_u_v;   //!
   TBranch        *b_shr_tkfit_pitch_v_v;   //!
   TBranch        *b_shr_tkfit_pitch_y_v;   //!
   TBranch        *b_shr_tkfit_dedx_u_v;   //!
   TBranch        *b_shr_tkfit_dedx_v_v;   //!
   TBranch        *b_shr_tkfit_dedx_y_v;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_u_v;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_v_v;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_y_v;   //!
   TBranch        *b_shr_tkfit_dedx_nhits_u_v;   //!
   TBranch        *b_shr_tkfit_dedx_nhits_v_v;   //!
   TBranch        *b_shr_tkfit_dedx_nhits_y_v;   //!
   TBranch        *b_shr_llr_pid_u_v;   //!
   TBranch        *b_shr_llr_pid_v_v;   //!
   TBranch        *b_shr_llr_pid_y_v;   //!
   TBranch        *b_shr_llr_pid_v;   //!
   TBranch        *b_shr_llr_pid_score_v;   //!
   TBranch        *b_shr_moliere_avg_v;   //!
   TBranch        *b_shr_moliere_rms_v;   //!
   TBranch        *b_evnunhits;   //!
   TBranch        *b_evlepnhits;   //!
   TBranch        *b_evpronhits;   //!
   TBranch        *b_evpi1nhits;   //!
   TBranch        *b_evpi0nhits;   //!
   TBranch        *b_evneunhits;   //!
   TBranch        *b_evgamnhits;   //!
   TBranch        *b_evothnhits;   //!
   TBranch        *b_slnunhits;   //!
   TBranch        *b_sllepnhits;   //!
   TBranch        *b_slpronhits;   //!
   TBranch        *b_slpi1nhits;   //!
   TBranch        *b_slpi0nhits;   //!
   TBranch        *b_slneunhits;   //!
   TBranch        *b_slgamnhits;   //!
   TBranch        *b_slothnhits;   //!
   TBranch        *b_pfnunhits;   //!
   TBranch        *b_pflepnhits;   //!
   TBranch        *b_pfpronhits;   //!
   TBranch        *b_pfpi1nhits;   //!
   TBranch        *b_pfpi0nhits;   //!
   TBranch        *b_pfneunhits;   //!
   TBranch        *b_pfgamnhits;   //!
   TBranch        *b_pfothnhits;   //!
   TBranch        *b_nu_completeness_from_pfp;   //!
   TBranch        *b_nu_purity_from_pfp;   //!
   TBranch        *b_trk_bragg_p_v;   //!
   TBranch        *b_trk_bragg_mu_v;   //!
   TBranch        *b_trk_bragg_mip_v;   //!
   TBranch        *b_trk_pida_v;   //!
   TBranch        *b_trk_pid_chipr_v;   //!
   TBranch        *b_trk_pid_chipi_v;   //!
   TBranch        *b_trk_pid_chika_v;   //!
   TBranch        *b_trk_pid_chimu_v;   //!
   TBranch        *b_trk_bragg_p_u_v;   //!
   TBranch        *b_trk_bragg_mu_u_v;   //!
   TBranch        *b_trk_bragg_mip_u_v;   //!
   TBranch        *b_trk_pida_u_v;   //!
   TBranch        *b_trk_pid_chipr_u_v;   //!
   TBranch        *b_trk_pid_chipi_u_v;   //!
   TBranch        *b_trk_pid_chika_u_v;   //!
   TBranch        *b_trk_pid_chimu_u_v;   //!
   TBranch        *b_trk_bragg_p_v_v;   //!
   TBranch        *b_trk_bragg_mu_v_v;   //!
   TBranch        *b_trk_bragg_mip_v_v;   //!
   TBranch        *b_trk_pida_v_v;   //!
   TBranch        *b_trk_pid_chipr_v_v;   //!
   TBranch        *b_trk_pid_chipi_v_v;   //!
   TBranch        *b_trk_pid_chika_v_v;   //!
   TBranch        *b_trk_pid_chimu_v_v;   //!
   TBranch        *b_trk_pfp_id_v;   //!
   TBranch        *b_trk_dir_x_v;   //!
   TBranch        *b_trk_dir_y_v;   //!
   TBranch        *b_trk_dir_z_v;   //!
   TBranch        *b_trk_start_x_v;   //!
   TBranch        *b_trk_start_y_v;   //!
   TBranch        *b_trk_start_z_v;   //!
   TBranch        *b_trk_sce_start_x_v;   //!
   TBranch        *b_trk_sce_start_y_v;   //!
   TBranch        *b_trk_sce_start_z_v;   //!
   TBranch        *b_trk_end_x_v;   //!
   TBranch        *b_trk_end_y_v;   //!
   TBranch        *b_trk_end_z_v;   //!
   TBranch        *b_trk_sce_end_x_v;   //!
   TBranch        *b_trk_sce_end_y_v;   //!
   TBranch        *b_trk_sce_end_z_v;   //!
   TBranch        *b_trk_distance_v;   //!
   TBranch        *b_trk_theta_v;   //!
   TBranch        *b_trk_phi_v;   //!
   TBranch        *b_trk_len_v;   //!
   TBranch        *b_trk_mcs_muon_mom_v;   //!
   TBranch        *b_trk_range_muon_mom_v;   //!
   TBranch        *b_trk_energy_proton_v;   //!
   TBranch        *b_trk_energy_muon_v;   //!
   TBranch        *b_trk_calo_energy_u_v;   //!
   TBranch        *b_trk_calo_energy_v_v;   //!
   TBranch        *b_trk_calo_energy_y_v;   //!
   TBranch        *b_trk_llr_pid_u_v;   //!
   TBranch        *b_trk_llr_pid_v_v;   //!
   TBranch        *b_trk_llr_pid_y_v;   //!
   TBranch        *b_trk_llr_pid_v;   //!
   TBranch        *b_trk_llr_pid_score_v;   //!
   TBranch        *b_bdt_nuNCpi0;   //!
   TBranch        *b_bdt_numuCCpi0;   //!
   TBranch        *b_bdt_numuCC;   //!
   TBranch        *b_bdt_ext;   //!
   TBranch        *b_bdt_cosmic;   //!
   TBranch        *b_bdt_global;   //!
   //TBranch        *b_bdt_global;   //!
   TBranch        *b_bdt_pi0_np;   //!
   TBranch        *b_bdt_nonpi0_np;   //!
   TBranch        *b_bdt_bkg_0p;   //!
   TBranch        *b_anglediff_Y;   //!
   TBranch        *b_anglediff_V;   //!
   TBranch        *b_anglediff_U;   //!
   TBranch        *b_trkpid;   //!

   twoproton_nuwro_overlay(TTree *tree=0);
   virtual ~twoproton_nuwro_overlay();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     Which_Run();
   virtual void     Define_Histograms(); //defines histograms. works for all samples
   virtual void     Fill_Histograms_Mine(int i, double wgt, int mc_n_threshold_muon, int mc_n_threshold_proton, int mc_n_threshold_pion0, int mc_n_threshold_pionpm, bool fv);
   virtual void     Fill_Histograms_Raquel(int i, double wgt, bool fv);
   virtual void     Fill_Track_Plots(int which_cut,float value, int pdg, bool contained_start,bool contained_end, double wgt); //fills the track variables 
   virtual void     Fill_Histograms_Particles(int mc_n_threshold_muon, int mc_n_threshold_proton, int mc_n_threshold_pion0, double mc_n_threshold_pionpm, bool fv,TVector3 vMuon, TLorentzVector muon, TVector3 vLead, TLorentzVector lead, TVector3 vRec, TLorentzVector rec, double wgt);
   virtual void     Fill_Histograms_Particles_Raquel(TVector3 vMuon, TLorentzVector muon, TVector3 vLead, TLorentzVector lead, TVector3 vRec, TLorentzVector rec, double wgt, bool fv);
   virtual void     Fill_Mine(int i, int j, double wgt);
   virtual void     Fill_Raquel(int i, int j, double wgt);
   virtual void     Fill_Particles(int j, TVector3 vMuon, TLorentzVector muon, TVector3 vLead, TLorentzVector lead, TVector3 vRec, TLorentzVector rec, double wgt);
   virtual void     Fill_Particles_Raquel(int j, TVector3 vMuon, TLorentzVector muon, TVector3 vLead, TLorentzVector lead, TVector3 vRec, TLorentzVector rec, double wgt);
   virtual void     Write_Histograms();

 private:

   //Stuff to determine which Run this is:
   ////////////////////////////////////////
   char response;
   const char* directory;
   const char* file;
   double pot_wgt;

   //Defining all the histograms becaues I wrote the classes stupidly
   //////////////////////////////////////////////////////////////////
 
   //PFP Specific Histograms                                                                                              
   static const int num = 4;
   const char * total[num] = {"npfp","vtx_npfp","ntrack","nshower"};
   TH1D * h_pfp_overlay[num]; //overlay
   TH1D* h_pfp[num]; //bnb, ext, dirt

  //Correlation Histograms for Vertex                                                                                       
  static const int num2d = 3;
  const char * total2d[num2d] = {"reco","truth","truth_sce"};
  const char * labelx[num2d] = {"reco x","true x","true x + sce"};
  const char * labely[num2d] = {"reco y","true y","true y + sce"};
  TH2D *h_correlation_overlay[num2d];
      
  //Varaibles I wanted to check after each cut
  static const int  number = 7; //number cuts                                                                        
  static const int  number2 = 11; //categories I defined                                                            
  static const int  number3 = 10; //categories raquel defined
  const char * point[number] ={"_before_selection","_after_fv","_after_three_pfps","_after_track_cut","_after_connection_cut","_after_pid","_after_reco_mom"}; //this defines histograms after each cut    
  const char * channel[number2]={"_total","_cc0p0pi","_cc1p0pi","_cc2p0pi","_ccNp0pi",
				 "_ccNp1pi","_ccNpNpi","_ccnue","_outfv","_nc","_other"}; //these are the channels I defined        
  const char * channel2[number3] = {"_total","_ccQE","_ccCOH","_ccMEC","_ccRES","_ccDIS",
				     "_ccNue","_nc","_outfv","_other"};//these are the channels that raquel defined
  
  TH1D* h_vtx_x_overlay[number][number2]; //reco x: overlay
  TH1D* h_vtx_y_overlay[number][number2]; //reco y: overlay
  TH1D* h_vtx_z_overlay[number][number2]; //reco z: overlay
  TH1D* h_vtx_x_raquel[number][number3]; //reco x: overlay                                                                          
  TH1D* h_vtx_y_raquel[number][number3]; //reco y: overlay                                                                        
  TH1D* h_vtx_z_raquel[number][number3]; //reco z: overlay                                                                       
  TH1D* h_topological_score_overlay[number][number2]; //overlay
  TH1D* h_topological_score_raquel[number][number3]; //raquel
  TH1D* h_cosmic_impact_parameter_overlay[number][number2]; //overlay
  TH1D* h_cosmic_impact_parameter_raquel[number][number3]; //raquel

  TH1D* h_vtx_x_mc[number][number2]; //mc x
  TH1D* h_vtx_y_mc[number][number2]; //mc y 
  TH1D* h_vtx_z_mc[number][number2]; //mc z
  TH1D* h_vtx_x_mc_sce[number][number2]; //mc+sce x
  TH1D* h_vtx_y_mc_sce[number][number2]; //mc+sce y
  TH1D* h_vtx_z_mc_sce[number][number2]; //mc+sce z
  TH1D* h_q2[number][number2]; //mc q2
  TH1D* h_X[number][number2]; //mc x
  TH1D* h_Y[number][number2]; //mc y
  TH1D* h_Pt[number][number2]; //mc Pt
 
  TH1D* h_vtx_x_mc_raquel[number][number3]; //mc x raquel                                                                                   
  TH1D* h_vtx_y_mc_raquel[number][number3]; //mc y                                                                                   
  TH1D* h_vtx_z_mc_raquel[number][number3]; //mc z                                                                                   
  TH1D* h_vtx_x_mc_sce_raquel[number][number3]; //mc+sce x                                                                          
  TH1D* h_vtx_y_mc_sce_raquel[number][number3]; //mc+sce y                                                                           
  TH1D* h_vtx_z_mc_sce_raquel[number][number3]; //mc+sce z                                                                           
  TH1D* h_q2_raquel[number][number3]; //mc q2                                                                                        
  TH1D* h_X_raquel[number][number3]; //mc x                                                                                          
  TH1D* h_Y_raquel[number][number3]; //mc y                                                                                          
  TH1D* h_Pt_raquel[number][number3]; //mc Pt   
  
  //Check of MEC weight
  TH1D* h_mec_wgt[number];
 
  //Track related variables after certain cuts
  static const int num_track = 4;
  const char* variable[num_track] = {"_track_score","_track_vertex_distance","_track_length","_track_pid"};
  static const int track_cut = 3;
  const char* which_track_cut[track_cut] = {"_after_3_pfps","_after_track_score","_after_distance_cut"};
  static const int num_part = 11;
  const char* particle[num_part] = {"_total","_proton_contained","_proton_uncontained","_muon_contained","_muon_uncontained","_pionpm","_pion0","_electron","_gamma","_kaon","_other"};
  TH1D* h_track_overlay[num_track][track_cut][num_part]; //overlay
  int num_bins_track[num_track] = {30,10,50,50};
  double xlim_low_track[num_track] = {0.0,0.0,0.0,-1.0};
  double xlim_high_track[num_track] = {1.0,10.0,50.0,1.0};
 
  //Effieincy Plots and Migration Matrices of XSec Variables
  ////////////////////////////////////////////////////////////
  TGraph* eff_graph = new TGraph(number); //efficiency as function of cuts                                                                                                                                                           
  TGraph* pur_graph = new TGraph(number); //efficiency as function of purity 
  TH2D* h_test_recoil;
  TH2D* h_test_leading;
  int recoil_theta = 0;
  int leading_theta = 0;

  //All the single particle plots
  static const int num_var = 4;
  const char* var[num_var] = {"_mom","_E","_theta","_phi"};
  int num_bins[num_var] = {50,50,30,10};
  double xlim_low[num_var] = {0,0,-1.5,-3.15}; //0.2 normally first -1.5                                                        
  double xlim_high_recoil[num_var] = {0.8,0.35,1.5,3.15};
  double xlim_high_leading[num_var] = {1.2,0.6,1.5,3.15}; //1.5 normally in first, 1.2                                          
  double xlim_high_muon[num_var]={2.5,1,1.5,3.15}; //2.5 first, 1.5 third     
  const char* xlabel[num_var] ={"P [GeV/c]","E [GeV]","cos(#theta)","#phi [Rad]"};
  TH1D* h_muon_overlay[num_var][number2]; //overlay
  TH1D* h_recoil_overlay[num_var][number2];
  TH1D* h_leading_overlay[num_var][number2];
  TH1D* h_muon_raquel[num_var][number3];
  TH1D* h_recoil_raquel[num_var][number3];
  TH1D* h_leading_raquel[num_var][number3];
    
  //Interesting Physics plots
  ////////////////////////////

  //mine
  TH1D* h_opening_angle_protons_lab_overlay[number2];
  TH1D* h_opening_angle_protons_com_overlay[number2];
  TH1D* h_opening_angle_mu_leading_overlay[number2];
  TH1D* h_opening_angle_mu_both_overlay[number2];
  TH1D* h_delta_PT_overlay[number2];
  TH1D* h_delta_alphaT_overlay[number2];
  TH1D* h_delta_phiT_overlay[number2];
  TH1D* h_pn_overlay[number2];
  TH1D* h_nu_E_overlay[number2];
  TH1D* h_mom_struck_nuc_overlay[number2];
  TH1D* h_tot_pz_overlay[number2];
  TH1D* h_tot_E_overlay[number2];
  TH1D* h_tot_E_minus_beam_overlay[number2];
  TH1D* h_E_resolution_overlay[number2];
  TH1D* h_PT_squared_overlay[number2];

  //raquel
  TH1D* h_opening_angle_protons_lab_raquel[number3];
  TH1D* h_opening_angle_protons_com_raquel[number3];
  TH1D* h_opening_angle_mu_leading_raquel[number3];
  TH1D* h_opening_angle_mu_both_raquel[number3];
  TH1D* h_delta_PT_raquel[number3];
  TH1D* h_delta_alphaT_raquel[number3];
  TH1D* h_delta_phiT_raquel[number3];
  TH1D* h_pn_raquel[number3];
  TH1D* h_nu_E_raquel[number3];
  TH1D* h_mom_struck_nuc_raquel[number3];
  TH1D* h_tot_pz_raquel[number3];
  TH1D* h_tot_E_raquel[number3];
  TH1D* h_tot_E_minus_beam_raquel[number3];
  TH1D* h_E_resolution_raquel[number3];
  TH1D* h_PT_squared_raquel[number3];

  vector<TH1*> h_list; //vector of all the 1D histograms
  vector<TH2*> h_list_2D; //vector of all the 2D histograms

  //defining class stuff
  helper_funcs cuts; //helper_funcs.h
  variables variables; //variables_funcs.h

  //Other parameters:                                                                                                                                                                                  
  double open_angle; //note this is the cos(opening angle)                                                                                                                                             
  double open_angle_mu; //note this is the cos(opening angle)                                                         
  double open_angle_mu_proton;
  double delta_pT; //stv delta_pT                                                                                                                                                                      
  double delta_alphaT; //stv delta_alphaT                                                                                                                                                              
  double delta_phiT; //stv delta_phiT                                                                                                                                                                  
  double cos_gamma_lab; //cos(opening angle) in lab                                                                                                                                                    
  double cos_gamma_cm; //cos(opening angle) in cm                                                                                                                                                      
  double En; //energy of struck nucleon                                                                                                                                                                
  double p_struck_nuc; //momentum of the struck nucleon                                                                                                                                                
  double pz_tot;
  double mec_wgt;

  //number of generated event/channel                                                                             
  int cc0p0pi[number+1] = {0}; //number-1 is after pid. number then would be particle specifics...should be same as number-1
  int cc1p0pi[number+1] = {0};
  int cc2p0pi[number+1] = {0};
  int ccNp0pi[number+1] = {0};
  int ccNp1pi[number+1] = {0};
  int ccNpNpi[number+1] = {0};
  int ccnue[number+1] = {0};
  int outfv[number+1] = {0};
  int nc[number+1] = {0};
  int other[number+1] = {0};
   
  //number of generated event/channel                                                                             
  int qel[number+1] = {0};
  int res[number+1] = {0};
  int mec[number+1] = {0};
  int coh[number+1] = {0};
  int dis[number+1] = {0};
  int ccnue_raquel[number+1] = {0};
  int outfv_raquel[number+1] = {0};
  int nc_raquel[number+1] = {0};
  int other_raquel[number+1] = {0};
  int res_count[4] = {0};
  
  int other_else = 0;
  int neutron = 0;
  int neutrino = 0;
  int zeros = 0;
  int total_protons = 0;
  int contain = 0;
  int uncontain = 0;

  int denom_contained =0; //checking number of events that are uncontained and contained
  int denom_uncontained = 0;
  int num_contained = 0;
  int num_uncontained = 0;

  int leading_true_angle_bad = 0;
  int recoil_true_angle_bad = 0;

};

#endif

#ifdef twoproton_nuwro_overlay_cxx
twoproton_nuwro_overlay::twoproton_nuwro_overlay(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  std::cout<<"Which Run Are we Looking at?"<<std::endl;
  std::cout<<" 1 = Run 1 \n 2 = Run 2 \n 3 = Run 3"<<std::endl;
  std::cin>>response;
  
  if(response =='1'){
    file = "/uboone/data/users/apapadop/searchingfornues/high_stat_prodgenie_bnb_nu_overlay_DetVar_Run1_NuWro_reco2_reco2.root";
  } else if(response == '2'){
    file = "/uboone/data/users/apapadop/searchingfornues/high_stat_prodgenie_bnb_nu_overlay_DetVar_Run2_NuWro_reco2_reco2.root";
  } else if(response == '3'){
    file = "/uboone/data/users/apapadop/searchingfornues/high_stat_prodgenie_bnb_nu_overlay_DetVar_Run3_NuWro_reco2_reco2.root";
  } else{
    std::cout<<"Invalid Response. Please Type 1, 2, or 3 for Run 1,Run 2, and Run 3 samples respectively."<<std::endl;
  }
  
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(Form("%s",file));
    if (!f || !f->IsOpen()) {
      f = new TFile(Form("%s",file));
    }
    TDirectory * dir = (TDirectory*)f->Get(Form("%s:/nuselection",file));
    dir->GetObject("NeutrinoSelectionFilter",tree);
  }
  Init(tree);

  //Run the program with the correct file
  Loop();
}

twoproton_nuwro_overlay::~twoproton_nuwro_overlay()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void twoproton_nuwro_overlay::Which_Run(){
  if(response =='1'){
    directory = "Run1";
    pot_wgt = 0.396;
  } else if(response == '2'){
    directory = "Run2";
    pot_wgt = 0.853;
  } else if(response == '3'){
    directory ="Run3";
    pot_wgt = 0.882;
  }  
} //end of which_run

void twoproton_nuwro_overlay::Define_Histograms(){

    //Total Histograms                                                                                                       
    for(int i=0; i < num; i++){
      h_pfp_overlay[i] = new TH1D(Form("h_%s_overlay",total[i]),Form("h_%s_overlay",total[i]),10,0,10);
      h_list.push_back(h_pfp_overlay[i]);
    }

    //Correlation Histograms                                                                                                           
    for(int i=0; i < num2d; i++){
      h_correlation_overlay[i] = new TH2D(Form("h_correlation_overlay_%s",total2d[i]),Form(";%s ;%s",labelx[i],labely[i]),40,0,275,40,-125,-125);
      h_list_2D.push_back(h_correlation_overlay[i]);
    }
  
    //Now to do the channel seperated histograms
    for(int i=0; i< number; i++){
      h_mec_wgt[i] = new TH1D(Form("h_mec_wgt%s",point[i]),Form("h_mec_wgt%s; Weight; Events",point[i]),50,-1,1);
      h_list.push_back(h_mec_wgt[i]);

      for(int j=0; j < number2; j++){ //my stuff                                                                                        
	h_vtx_x_overlay[i][j]=new TH1D(Form("h_vtx_x%s%s",point[i],channel[j]),Form("h_vtx_x%s%s",point[i],channel[j]),50,0,250);
	h_vtx_y_overlay[i][j]=new TH1D(Form("h_vtx_y%s%s",point[i],channel[j]),Form("h_vtx_y%s%s",point[i],channel[j]),50,-125,125);
	h_vtx_z_overlay[i][j]=new TH1D(Form("h_vtx_z%s%s",point[i],channel[j]),Form("h_vtx_z%s%s",point[i],channel[j]),50,0,1050);
	h_vtx_x_mc[i][j]=new TH1D(Form("h_vtx_x_mc%s%s",point[i],channel[j]),Form("h_vtx_x_mc%s%s",point[i],channel[j]),40,0,275);
	h_vtx_y_mc[i][j]=new TH1D(Form("h_vtx_y_mc%s%s",point[i],channel[j]),Form("h_vtx_y_mc%s%s",point[i],channel[j]),40,-125,125);
	h_vtx_z_mc[i][j]=new TH1D(Form("h_vtx_z_mc%s%s",point[i],channel[j]),Form("h_vtx_z_mc%s%s",point[i],channel[j]),50,0,1050);
	h_vtx_x_mc_sce[i][j]=new TH1D(Form("h_vtx_x_mc_sce%s%s",point[i],channel[j]),Form("h_vtx_x_mc_sce%s%s",point[i],channel[j]),40,0,275);
	h_vtx_y_mc_sce[i][j]=new TH1D(Form("h_vtx_y_mc_sce%s%s",point[i],channel[j]),Form("h_vtx_y_mc_sce%s%s",point[i],channel[j]),40,-125,125);
	h_vtx_z_mc_sce[i][j]=new TH1D(Form("h_vtx_z_mc_sce%s%s",point[i],channel[j]),Form("h_vtx_z_mc_sce%s%s",point[i],channel[j]),50,0,1050);
	h_q2[i][j] = new TH1D(Form("h_q2%s%s",point[i],channel[j]),Form("h_q2_x%s%s",point[i],channel[j]),20,0,2);
	h_X[i][j] = new TH1D(Form("h_X%s%s",point[i],channel[j]),Form("h_X_x%s%s",point[i],channel[j]),20,0,2);
	h_Y[i][j] = new TH1D(Form("h_Y%s%s",point[i],channel[j]),Form("h_Y_x%s%s",point[i],channel[j]),10,0,1);
	h_Pt[i][j] = new TH1D(Form("h_Pt%s%s",point[i],channel[j]),Form("h_Pt_x%s%s",point[i],channel[j]),20,0,2);
	h_cosmic_impact_parameter_overlay[i][j] = new TH1D(Form("h_cosmic_impact_parameter%s%s",point[i],channel[j]),Form("h_cosmic_impact_parameter%s%s; Cosmic Impact Distance (cm); No. Events",point[i],channel[j]),20,0,200);
	h_topological_score_overlay[i][j] = new TH1D(Form("h_topological_score%s%s",point[i],channel[j]),Form("h_topological_score%s%s; Topological Score; No. Events",point[i],channel[j]),50,0.0,1.0); //30 for wouter, 50 for steven

	h_list.push_back(h_topological_score_overlay[i][j]);
	h_list.push_back(h_cosmic_impact_parameter_overlay[i][j]);
	h_list.push_back(h_vtx_x_overlay[i][j]);
	h_list.push_back(h_vtx_y_overlay[i][j]);
	h_list.push_back(h_vtx_z_overlay[i][j]);
	h_list.push_back(h_vtx_x_mc[i][j]);
	h_list.push_back(h_vtx_y_mc[i][j]);
	h_list.push_back(h_vtx_z_mc[i][j]);
	h_list.push_back(h_vtx_x_mc_sce[i][j]);
	h_list.push_back(h_vtx_y_mc_sce[i][j]);
	h_list.push_back(h_vtx_z_mc_sce[i][j]);
	h_list.push_back(h_q2[i][j]);
	h_list.push_back(h_X[i][j]);
	h_list.push_back(h_Y[i][j]);
	h_list.push_back(h_Pt[i][j]);

      }

      for(int j=0; j < number3; j++){//raquel's stuff
	h_vtx_x_raquel[i][j]=new TH1D(Form("h_vtx_x_raquel%s%s",point[i],channel2[j]),Form("h_vtx_x_raquel%s%s",point[i],channel2[j]),50,0,250);
	h_vtx_y_raquel[i][j]=new TH1D(Form("h_vtx_y_raquel%s%s",point[i],channel2[j]),Form("h_vtx_y_raquel%s%s",point[i],channel2[j]),50,-125,125);
	h_vtx_z_raquel[i][j]=new TH1D(Form("h_vtx_z_raquel%s%s",point[i],channel2[j]),Form("h_vtx_z_raquel%s%s",point[i],channel2[j]),50,0,1050);
	h_vtx_x_mc_raquel[i][j]=new TH1D(Form("h_vtx_x_m_raquel%s%s",point[i],channel2[j]),Form("h_vtx_x_mc_raquel%s%s",point[i],channel2[j]),40,0,275);
	h_vtx_y_mc_raquel[i][j]=new TH1D(Form("h_vtx_y_mc_raquel%s%s",point[i],channel2[j]),Form("h_vtx_y_mc_raquel%s%s",point[i],channel2[j]),40,-125,125);
	h_vtx_z_mc_raquel[i][j]=new TH1D(Form("h_vtx_z_mc_raquel%s%s",point[i],channel2[j]),Form("h_vtx_z_mc_raquel%s%s",point[i],channel2[j]),50,0,1050);
	h_vtx_x_mc_sce_raquel[i][j]=new TH1D(Form("h_vtx_x_mc_sce_raquel%s%s",point[i],channel2[j]),Form("h_vtx_x_mc_sce_raquel%s%s",point[i],channel2[j]),40,0,275);
	h_vtx_y_mc_sce_raquel[i][j]=new TH1D(Form("h_vtx_y_mc_sce-raquel%s%s",point[i],channel2[j]),Form("h_vtx_y_mc_sce_raquel%s%s",point[i],channel2[j]),40,-125,125);
	h_vtx_z_mc_sce_raquel[i][j]=new TH1D(Form("h_vtx_z_mc_sce_raquel%s%s",point[i],channel2[j]),Form("h_vtx_z_mc_sce_raquel%s%s",point[i],channel2[j]),50,0,1050);
	h_q2_raquel[i][j] = new TH1D(Form("h_q2_raquel%s%s",point[i],channel2[j]),Form("h_q2_raquel%s%s",point[i],channel2[j]),20,0,2);
	h_X_raquel[i][j] = new TH1D(Form("h_X_raquel%s%s",point[i],channel2[j]),Form("h_X_raquel%s%s",point[i],channel2[j]),20,0,2);
	h_Y_raquel[i][j] = new TH1D(Form("h_Y_raquel%s%s",point[i],channel2[j]),Form("h_Y_raquel%s%s",point[i],channel2[j]),10,0,1);
	h_Pt_raquel[i][j] = new TH1D(Form("h_Pt_raquel%s%s",point[i],channel2[j]),Form("h_Pt_raquel%s%s",point[i],channel2[j]),20,0,2);
	h_cosmic_impact_parameter_raquel[i][j] = new TH1D(Form("h_cosmic_impact_parameter_raquel%s%s",point[i],channel2[j]),Form("h_cosmic_impact_parameter_raquel%s%s; Cosmic Impact Distance (cm); No. Events",point[i],channel2[j]),20,0,200);
	h_topological_score_raquel[i][j] = new TH1D(Form("h_topological_score_raquel%s%s",point[i],channel2[j]),Form("h_topological_score_raquel%s%s; Topological Score; No. Events",point[i],channel2[j]),50,0.0,1.0); //30 wouter/50 steven

	h_list.push_back(h_topological_score_raquel[i][j]);
	h_list.push_back(h_cosmic_impact_parameter_raquel[i][j]);
	h_list.push_back(h_vtx_x_raquel[i][j]);
	h_list.push_back(h_vtx_y_raquel[i][j]);
	h_list.push_back(h_vtx_z_raquel[i][j]);
	h_list.push_back(h_vtx_x_mc_raquel[i][j]);
	h_list.push_back(h_vtx_y_mc_raquel[i][j]);
	h_list.push_back(h_vtx_z_mc_raquel[i][j]);
	h_list.push_back(h_vtx_x_mc_sce_raquel[i][j]);
	h_list.push_back(h_vtx_y_mc_sce_raquel[i][j]);
	h_list.push_back(h_vtx_z_mc_sce_raquel[i][j]);
	h_list.push_back(h_q2_raquel[i][j]);
	h_list.push_back(h_X_raquel[i][j]);
	h_list.push_back(h_Y_raquel[i][j]);
	h_list.push_back(h_Pt_raquel[i][j]);

      }
    }

    //track variables
    for(int j=0; j < num_track; j++){
      for(int k=0; k < track_cut; k++){
	for(int i=0; i < num_part; i++){
	  h_track_overlay[j][k][i] = new TH1D(Form("h_track%s%s%s",variable[j],which_track_cut[k],particle[i]),Form("h_track%s%s%s",variable[j],which_track_cut[k],particle[i]),num_bins_track[j],xlim_low_track[j],xlim_high_track[j]);
	  h_list.push_back(h_track_overlay[j][k][i]);
	}
      }
    }
    
    //Efficiency plots
    ///////////////////////////////////////
    h_test_recoil = new TH2D("h_test_recoil","h_test_recoil",50,-10.0,10.0,50,0.25,0.85);
    h_test_leading = new TH2D("h_test_leading","h_test_leading",50,-10.0,10.0,50,0.25,1.2);
    h_list_2D.push_back(h_test_recoil);
    h_list_2D.push_back(h_test_leading);

    //Particle specific plots
    /////////////////////////
    for(int j = 0; j < num_var; j++){
      for(int k = 0; k < number2; k++){
	if(use_xsec_binning == true){
	  if(j == 0){//momentum
	    //muon
	    const Int_t bins_muon = 6; 
	    Double_t edges_muon[bins_muon+1] = {0.1,0.2,0.3,0.5,0.7,1.3,2.5};
	    h_muon_overlay[j][k] = new TH1D(Form("h_muon%s%s",var[j],channel[k]),Form(" h_muon%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),bins_muon,edges_muon);
	    //recoil
	    const Int_t bins_recoil = 5;
            Double_t edges_recoil[bins_recoil+1] = {0.25,0.35,0.45,0.55,0.65,1.2};
	    h_recoil_overlay[j][k] = new TH1D(Form("h_recoil%s%s",var[j],channel[k]),Form("h_recoil%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),bins_recoil,edges_recoil);
	    //leading
	    const Int_t bins_leading = 6;
            Double_t edges_leading[bins_leading+1] = {0.25,0.35,0.45,0.55,0.65,0.75,1.2};
	    h_leading_overlay[j][k] = new TH1D(Form("h_leading%s%s",var[j],channel[k]),Form("h_leading%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),bins_leading,edges_leading);

	  }else if (j == 2){//costheta
	    const Int_t bins_theta = 10;
	    Double_t edges_theta[bins_theta+1] = {-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0};
	    h_muon_overlay[j][k] = new TH1D(Form("h_muon%s%s",var[j],channel[k]),Form(" h_muon%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),bins_theta,edges_theta);
	    h_recoil_overlay[j][k] = new TH1D(Form("h_recoil%s%s",var[j],channel[k]),Form("h_recoil%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),bins_theta,edges_theta);
	    h_leading_overlay[j][k] = new TH1D(Form("h_leading%s%s",var[j],channel[k]),Form("h_leading%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),bins_theta,edges_theta);

	  }else{//phi and energy
	    h_muon_overlay[j][k] = new TH1D(Form("h_muon%s%s",var[j],channel[k]),Form(" h_muon%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_muon[j]);
            h_recoil_overlay[j][k] = new TH1D(Form("h_recoil%s%s",var[j],channel[k]),Form("h_recoil%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_recoil[j]);
            h_leading_overlay[j][k] = new TH1D(Form("h_leading%s%s",var[j],channel[k]),Form("h_leading%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_leading[j]);
	  }
	} else if (use_xsec_binning == false){
	  h_muon_overlay[j][k] = new TH1D(Form("h_muon%s%s",var[j],channel[k]),Form(" h_muon%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_muon[j]);
          h_recoil_overlay[j][k] = new TH1D(Form("h_recoil%s%s",var[j],channel[k]),Form("h_recoil%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_recoil[j]);
          h_leading_overlay[j][k] = new TH1D(Form("h_leading%s%s",var[j],channel[k]),Form("h_leading%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_leading[j]);
	} //end of false statement

	h_list.push_back(h_muon_overlay[j][k]);
	h_list.push_back(h_recoil_overlay[j][k]);
	h_list.push_back(h_leading_overlay[j][k]);

      } //end loop over channels
    } //end loop over variables

    for(int j = 0; j < num_var; j++){
      for(int k = 0; k < number3; k++){
	if(use_xsec_binning == true){
	  if(j == 0){ //momentum
	    //muon
	    const Int_t bins_muon = 6;
            Double_t edges_muon[bins_muon+1] = {0.1,0.2,0.3,0.5,0.7,1.3,2.5};
	    h_muon_raquel[j][k] = new TH1D(Form("h_muon_raquel%s%s",var[j],channel2[k]),Form(" h_muon_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),bins_muon,edges_muon);
	    //recoil
	    const Int_t bins_recoil = 5;
            Double_t edges_recoil[bins_recoil+1] = {0.25,0.35,0.45,0.55,0.65,1.2};
	    h_recoil_raquel[j][k] = new TH1D(Form("h_recoil_raquel%s%s",var[j],channel2[k]),Form("h_recoil_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),bins_recoil,edges_recoil);
	    //leading
	    const Int_t bins_leading = 6;
            Double_t edges_leading[bins_leading+1] = {0.25,0.35,0.45,0.55,0.65,0.75,1.2};
	    h_leading_raquel[j][k] = new TH1D(Form("h_leading_raquel%s%s",var[j],channel2[k]),Form("h_leading_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),bins_leading,edges_leading);

	  } else if (j == 2){ //theta
	    const Int_t bins_theta = 10;
            Double_t edges_theta[bins_theta+1] = {-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0};
	    h_muon_raquel[j][k] = new TH1D(Form("h_muon_raquel%s%s",var[j],channel2[k]),Form(" h_muon_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),bins_theta,edges_theta);
	    h_recoil_raquel[j][k] = new TH1D(Form("h_recoil_raquel%s%s",var[j],channel2[k]),Form("h_recoil_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),bins_theta,edges_theta);
	    h_leading_raquel[j][k] = new TH1D(Form("h_leading_raquel%s%s",var[j],channel2[k]),Form("h_leading_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),bins_theta,edges_theta);

	  } else{ //phi and energy
	    h_muon_raquel[j][k] = new TH1D(Form("h_muon_raquel%s%s",var[j],channel2[k]),Form(" h_muon_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_muon[j]);
	    h_recoil_raquel[j][k] = new TH1D(Form("h_recoil_raquel%s%s",var[j],channel2[k]),Form("h_recoil_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_recoil[j]);
	    h_leading_raquel[j][k] = new TH1D(Form("h_leading_raquel%s%s",var[j],channel2[k]),Form("h_leading_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_leading[j]);
	  }
	}else if(use_xsec_binning == false){
	  h_muon_raquel[j][k] = new TH1D(Form("h_muon_raquel%s%s",var[j],channel2[k]),Form(" h_muon_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_muon[j]);
	  h_recoil_raquel[j][k] = new TH1D(Form("h_recoil_raquel%s%s",var[j],channel2[k]),Form("h_recoil_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_recoil[j]);
	  h_leading_raquel[j][k] = new TH1D(Form("h_leading_raquel%s%s",var[j],channel2[k]),Form("h_leading_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_leading[j]);
	} //end of false statement

	h_list.push_back(h_muon_raquel[j][k]);
	h_list.push_back(h_recoil_raquel[j][k]);
	h_list.push_back(h_leading_raquel[j][k]);

      } //end of loop over channels
    } //end loop over variables

    //more particle specific plots
    for(int i = 0; i < number2; i++){
      if(use_xsec_binning == true){
	const Int_t bins = 10;
	Double_t edges[bins+1] = {-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0};
	h_opening_angle_protons_lab_overlay[i] = new TH1D(Form("h_opening_angle_protons_lab%s",channel[i]),Form("h_opening_angle_protons_lab%s; Opening Angle btwn Two Protons; Counts",channel[i]),bins,edges); //50, 0, 1.5                           
	h_opening_angle_protons_com_overlay[i] = new TH1D(Form("h_opening_angle_protons_com%s",channel[i]),Form("h_opening_angle_protons_com%s; cos(#gamma_{COM}); Counts",channel[i]),bins,edges);
	h_opening_angle_mu_leading_overlay[i] = new TH1D(Form("h_opening_angle_mu_leading%s",channel[i]),Form("h_opening_angle_mu_leading%s;Opening Angle btwn Muon and Leading Proton; Counts",channel[i]),bins,edges);
	h_opening_angle_mu_both_overlay[i] = new TH1D(Form("h_opening_angle_mu_both%s",channel[i]),Form("h_opening_angle_mu_both%s; Opening Angle btwn Muon and Total Proton Momentum; Counts",channel[i]),bins,edges);

	const Int_t bins_stv_mom = 4;
        Double_t edges_stv_mom[bins_stv_mom+1] = {0,0.2,0.4,0.6,1.0};
	h_delta_PT_overlay[i] = new TH1D(Form("h_delta_PT%s",channel[i]),Form("h_deltaPT%s;#delta P_{T} [GeV/c];Counts",channel[i]),bins_stv_mom,edges_stv_mom);  
	
	const Int_t bins_stv_angles = 6;
        Double_t edges_stv_angles[bins_stv_angles+1] = {0,30,60,90,120,150,180};
	h_delta_alphaT_overlay[i] = new TH1D(Form("h_delta_alphaT%s",channel[i]),Form("h_delta_alphaT%s; #delta #alpha_{T} [Deg.];Counts",channel[i]),bins_stv_angles,edges_stv_angles); //0,180                     
	h_delta_phiT_overlay[i] = new TH1D(Form("h_delta_phiT%s",channel[i]),Form("h_delta_phiT%s; #delta #phi_{T} [Deg.];Counts",channel[i]),bins_stv_angles,edges_stv_angles);                                                               

	const Int_t bins_neutron = 20;
	Double_t edges_neutron[bins_neutron+1] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
	h_pn_overlay[i] = new TH1D(Form("h_pn%s",channel[i]),Form("h_pn%s; p_{n} [GeV/c];Counts",channel[i]),bins_neutron,edges_neutron);
  
	const Int_t bins_nuE = 6;
        Double_t edges_nuE[bins_nuE+1] = {0,0.3,0.5,0.7,0.9,1.2,4.0};
	h_nu_E_overlay[i] =  new TH1D(Form("h_nu_E%s",channel[i]),Form("h_nu_E%s; Total Energy; Counts;",channel[i]),bins_nuE,edges_nuE);

	//basic bitches
	h_mom_struck_nuc_overlay[i] = new TH1D(Form("h_mom_struck_nuc%s",channel[i]),Form("h_mom_struck_nuc%s; P_{Init}; Counts", channel[i]),30, 0, 1);
	h_tot_pz_overlay[i] = new TH1D(Form("h_tot_pz%s",channel[i]),Form("h_tot_pz%s; P_{Z}^{Total}; Counts",channel[i]), 20, 0, 2);
	h_tot_E_overlay[i] = new TH1D(Form("h_tot_E%s",channel[i]),Form("h_tot_E%s; Total Energy; Counts;",channel[i]),50,0,2.5);
	h_tot_E_minus_beam_overlay[i] = new TH1D(Form("h_tot_E_minus_beam%s",channel[i]),Form("h_tot_E_minus_beam%s; Total Energy Remaining (MeV/c); Counts;",channel[i]),100,-100,0);
	h_E_resolution_overlay[i] = new TH1D(Form("h_E_resolution%s",channel[i]),Form("h_E_resolution%s; Energy Resolution (GeV/c); Counts",channel[i]),100,-1.0,1.0);
	h_PT_squared_overlay[i] = new TH1D(Form("h_PT_squared%s",channel[i]),Form("h_PT_squared%s; P_{T}^{2}; Counts", channel[i]),50,0,5);

      }else if( use_xsec_binning == false){
	h_opening_angle_protons_lab_overlay[i] = new TH1D(Form("h_opening_angle_protons_lab%s",channel[i]),Form("h_opening_angle_protons_lab%s; Opening Angle btwn Two Protons; Counts",channel[i]),30,-1.5,1.5); //50, 0, 1.5        
	h_opening_angle_protons_com_overlay[i] = new TH1D(Form("h_opening_angle_protons_com%s",channel[i]),Form("h_opening_angle_protons_com%s; cos(#gamma_{COM}); Counts",channel[i]),30,-1.5,1.5);
	h_opening_angle_mu_leading_overlay[i] = new TH1D(Form("h_opening_angle_mu_leading%s",channel[i]),Form("h_opening_angle_mu_leading%s;Opening Angle btwn Muon and Leading Proton; Counts",channel[i]),30,-1.5,1.5);
	h_opening_angle_mu_both_overlay[i] = new TH1D(Form("h_opening_angle_mu_both%s",channel[i]),Form("h_opening_angle_mu_both%s; Opening Angle btwn Muon and Total Proton Momentum; Counts",channel[i]),30,-1.5,1.5);
	h_delta_PT_overlay[i] = new TH1D(Form("h_delta_PT%s",channel[i]),Form("h_deltaPT%s;#delta P_{T} [GeV/c];Counts",channel[i]),15,0,1); //normally 10 bins
	h_delta_alphaT_overlay[i] = new TH1D(Form("h_delta_alphaT%s",channel[i]),Form("h_delta_alphaT%s; #delta #alpha_{T} [Deg.];Counts",channel[i]),10,0,180); //0,180                 
	h_delta_phiT_overlay[i] = new TH1D(Form("h_delta_phiT%s",channel[i]),Form("h_delta_phiT%s; #delta #phi_{T} [Deg.];Counts",channel[i]),10,0,180); //0,180                   
	h_pn_overlay[i] = new TH1D(Form("h_pn%s",channel[i]),Form("h_pn%s; p_{n} [GeV/c];Counts",channel[i]),50,0,1.0);
	h_nu_E_overlay[i] =  new TH1D(Form("h_nu_E%s",channel[i]),Form("h_nu_E%s; Total Energy; Counts;",channel[i]),50,0,2.5);
	h_mom_struck_nuc_overlay[i] = new TH1D(Form("h_mom_struck_nuc%s",channel[i]),Form("h_mom_struck_nuc%s; P_{Init}; Counts", channel[i]),30, 0, 1);
	h_tot_pz_overlay[i] = new TH1D(Form("h_tot_pz%s",channel[i]),Form("h_tot_pz%s; P_{Z}^{Total}; Counts",channel[i]), 20, 0, 2);
	h_tot_E_overlay[i] = new TH1D(Form("h_tot_E%s",channel[i]),Form("h_tot_E%s; Total Energy; Counts;",channel[i]),50,0,2.5);
	h_tot_E_minus_beam_overlay[i] = new TH1D(Form("h_tot_E_minus_beam%s",channel[i]),Form("h_tot_E_minus_beam%s; Total Energy Remaining (MeV/c); Counts;",channel[i]),100,-100,0);
	h_E_resolution_overlay[i] = new TH1D(Form("h_E_resolution%s",channel[i]),Form("h_E_resolution%s; Energy Resolution (GeV/c); Counts",channel[i]),100,-1.0,1.0);
	h_PT_squared_overlay[i] = new TH1D(Form("h_PT_squared%s",channel[i]),Form("h_PT_squared%s; P_{T}^{2}; Counts", channel[i]),50,0,5);
      } //end of false statement

      h_list.push_back(h_PT_squared_overlay[i]);
      h_list.push_back(h_E_resolution_overlay[i]);
      h_list.push_back(h_opening_angle_mu_both_overlay[i]);
      h_list.push_back(h_nu_E_overlay[i]);
      h_list.push_back(h_tot_E_overlay[i]);
      h_list.push_back(h_tot_E_minus_beam_overlay[i]);
      h_list.push_back(h_mom_struck_nuc_overlay[i]);
      h_list.push_back(h_tot_pz_overlay[i]);
      h_list.push_back(h_opening_angle_protons_lab_overlay[i]);
      h_list.push_back(h_opening_angle_protons_com_overlay[i]);
      h_list.push_back(h_opening_angle_mu_leading_overlay[i]);
      h_list.push_back(h_delta_PT_overlay[i]);
      h_list.push_back(h_delta_alphaT_overlay[i]);
      h_list.push_back(h_delta_phiT_overlay[i]);
      h_list.push_back(h_pn_overlay[i]);

    }

    for(int i = 0; i < number3; i++){
      if(use_xsec_binning == true){

	const Int_t bins = 10;
        Double_t edges[bins+1] = {-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0};
	h_opening_angle_protons_lab_raquel[i] = new TH1D(Form("h_opening_angle_protons_lab_raquel%s",channel2[i]),Form("h_opening_angle_protons_lab_raquel%s; Opening Angle btwn Two Protons; Counts",channel2[i]),bins,edges); 
	h_opening_angle_protons_com_raquel[i] = new TH1D(Form("h_opening_angle_protons_com_raquel%s",channel2[i]),Form("h_opening_angle_protons_com_raquel%s; cos(#gamma_{COM}); Counts",channel2[i]),bins,edges);
        h_opening_angle_mu_leading_raquel[i] = new TH1D(Form("h_opening_angle_mu_leading_raquel%s",channel2[i]),Form("h_opening_angle_mu_leading_raquel%s;Opening Angle btwn Muon and Leading Proton; Counts",channel2[i]),bins,edges);
        h_opening_angle_mu_both_raquel[i] = new TH1D(Form("h_opening_angle_mu_both_raquel%s",channel2[i]),Form("h_opening_angle_mu_both_raquel%s; Opening Angle btwn Muon and Total Proton Momentum; Counts",channel2[i]),bins,edges);
        
	const Int_t bins_stv_mom = 4;
        Double_t edges_stv_mom[bins_stv_mom+1] = {0,0.2,0.4,0.6,1.0};
	h_delta_PT_raquel[i] = new TH1D(Form("h_delta_PT_raquel%s",channel2[i]),Form("h_deltaPT_raquel%s;#delta P_{T} [GeV/c];Counts",channel2[i]),bins_stv_mom,edges_stv_mom);                
 
	const Int_t bins_stv_angles = 6;
        Double_t edges_stv_angles[bins_stv_angles+1] = {0,30,60,90,120,150,180};
        h_delta_alphaT_raquel[i] = new TH1D(Form("h_delta_alphaT_raquel%s",channel2[i]),Form("h_delta_alphaT_raquel%s; #delta #alpha_{T} [Deg.];Counts",channel2[i]),bins_stv_angles,edges_stv_angles);
	h_delta_phiT_raquel[i] = new TH1D(Form("h_delta_phiT_raquel%s",channel2[i]),Form("h_delta_phiT_raquel%s; #delta #phi_{T} [Deg.];Counts",channel2[i]),bins_stv_angles,edges_stv_angles);                                            

	const Int_t bins_neutron = 20;
	Double_t edges_neutron[bins_neutron+1] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
	h_pn_raquel[i] = new TH1D(Form("h_pn_raquel%s",channel2[i]),Form("h_pn_raquel%s; p_{n} [GeV/c];Counts",channel2[i]),bins_neutron,edges_neutron);

	const Int_t bins_nuE = 6;
        Double_t edges_nuE[bins_nuE+1] = {0,0.3,0.5,0.7,0.9,1.2,4.0};
	h_nu_E_raquel[i] =  new TH1D(Form("h_nu_E_raquel%s",channel2[i]),Form("h_nu_E_raquel%s; Total Energy; Counts;",channel2[i]),bins_nuE,edges_nuE);

	//basic bitches
        h_mom_struck_nuc_raquel[i] = new TH1D(Form("h_mom_struck_nuc_raquel%s",channel2[i]),Form("h_mom_struck_nuc_raquel%s; P_{Init}; Counts", channel2[i]),30, 0, 1);
        h_tot_pz_raquel[i] = new TH1D(Form("h_tot_pz_raquel%s",channel2[i]),Form("h_tot_pz_raquel%s; P_{Z}^{Total}; Counts",channel2[i]), 20, 0, 2);
        h_tot_E_raquel[i] = new TH1D(Form("h_tot_E_raquel%s",channel2[i]),Form("h_tot_E_raquel%s; Total Energy; Counts;",channel2[i]),50,0,2.5);
        h_tot_E_minus_beam_raquel[i] = new TH1D(Form("h_tot_E_minus_beam_raquel%s",channel2[i]),Form("h_tot_E_minus_beam_raquel%s; Total Energy; Counts;",channel2[i]),100,-100,0);
        h_E_resolution_raquel[i] = new TH1D(Form("h_E_resolution_raquel%s",channel2[i]),Form("h_E_resolution_raquel%s; Energy Resolution (GeV/c); Counts",channel2[i]),100,-1.0,1.0);
        h_PT_squared_raquel[i] = new TH1D(Form("h_PT_squared_raquel%s",channel2[i]),Form("h_PT_squared_raquel%s; P_{T}^{2}; Counts", channel2[i]),50,0,5);

      }else if( use_xsec_binning == false){
	h_opening_angle_protons_lab_raquel[i] = new TH1D(Form("h_opening_angle_protons_lab_raquel%s",channel2[i]),Form("h_opening_angle_protons_lab_raquel%s; Opening Angle btwn Two Protons; Counts",channel2[i]),30,-1.5,1.5); //50, 0, 1.5           
	h_opening_angle_protons_com_raquel[i] = new TH1D(Form("h_opening_angle_protons_com_raquel%s",channel2[i]),Form("h_oopening_angle_protons_com__raquel%s; cos(#gamma_{COM}); Counts",channel2[i]),30,-1.5,1.5);
	h_opening_angle_mu_leading_raquel[i] = new TH1D(Form("h_opening_angle_mu_leading_raquel%s",channel2[i]),Form("h_opening_angle_mu_leading_raquel%s;Opening Angle btwn Muon and Leading Proton; Counts",channel2[i]),30,-1.5,1.5);
	h_opening_angle_mu_both_raquel[i] = new TH1D(Form("h_opening_angle_mu_both_raquel%s",channel2[i]),Form("h_opening_angle_mu_both_raquel%s; Opening Angle btwn Muon and Total Proton Momentum; Counts",channel2[i]),30,-1.5,1.5);
	h_delta_PT_raquel[i] = new TH1D(Form("h_delta_PT_raquel%s",channel2[i]),Form("h_deltaPT_raquel%s;#delta P_{T} [GeV/c];Counts",channel2[i]),15,0,1); //normally 10 bins                                                                    
	h_delta_alphaT_raquel[i] = new TH1D(Form("h_delta_alphaT_raquel%s",channel2[i]),Form("h_delta_alphaT_raquel%s; #delta #alpha_{T} [Deg.];Counts",channel2[i]),10,0,180); //0,180                                                           
	h_delta_phiT_raquel[i] = new TH1D(Form("h_delta_phiT_raquel%s",channel2[i]),Form("h_delta_phiT_raquel%s; #delta #phi_{T} [Deg.];Counts",channel2[i]),10,0,180); //0,180                                                                   
	h_pn_raquel[i] = new TH1D(Form("h_pn_raquel%s",channel2[i]),Form("h_pn_raquel%s; p_{n} [GeV/c];Counts",channel2[i]),50,0,1.0);
	h_nu_E_raquel[i] =  new TH1D(Form("h_nu_E_raquel%s",channel2[i]),Form("h_nu_E_raquel%s; Total Energy; Counts;",channel2[i]),50,0,2.5);
	h_mom_struck_nuc_raquel[i] = new TH1D(Form("h_mom_struck_nuc_raquel%s",channel2[i]),Form("h_mom_struck_nuc_raquel%s; P_{Init}; Counts", channel2[i]),30, 0, 1);
	h_tot_pz_raquel[i] = new TH1D(Form("h_tot_pz_raquel%s",channel2[i]),Form("h_tot_pz_raquel%s; P_{Z}^{Total}; Counts",channel2[i]), 20, 0, 2);
	h_tot_E_raquel[i] = new TH1D(Form("h_tot_E_raquel%s",channel2[i]),Form("h_tot_E_raquel%s; Total Energy; Counts;",channel2[i]),50,0,2.5);
	h_tot_E_minus_beam_raquel[i] = new TH1D(Form("h_tot_E_minus_beam_raquel%s",channel2[i]),Form("h_tot_E_minus_beam_raquel%s; Total Energy; Counts;",channel2[i]),100,-100,0);
	h_E_resolution_raquel[i] = new TH1D(Form("h_E_resolution_raquel%s",channel2[i]),Form("h_E_resolution_raquel%s; Energy Resolution (GeV/c); Counts",channel2[i]),100,-1.0,1.0);
	h_PT_squared_raquel[i] = new TH1D(Form("h_PT_squared_raquel%s",channel2[i]),Form("h_PT_squared_raquel%s; P_{T}^{2}; Counts", channel2[i]),50,0,5);
    } //end of false statement           

      h_list.push_back(h_PT_squared_raquel[i]);
      h_list.push_back(h_E_resolution_raquel[i]);
      h_list.push_back(h_opening_angle_mu_both_raquel[i]);
      h_list.push_back(h_nu_E_raquel[i]);
      h_list.push_back(h_tot_E_raquel[i]);
      h_list.push_back(h_tot_E_minus_beam_raquel[i]);
      h_list.push_back(h_mom_struck_nuc_raquel[i]);
      h_list.push_back(h_tot_pz_raquel[i]);
      h_list.push_back(h_opening_angle_protons_lab_raquel[i]);
      h_list.push_back(h_opening_angle_protons_com_raquel[i]);
      h_list.push_back(h_opening_angle_mu_leading_raquel[i]);
      h_list.push_back(h_delta_PT_raquel[i]);
      h_list.push_back(h_delta_alphaT_raquel[i]);
      h_list.push_back(h_delta_phiT_raquel[i]);
      h_list.push_back(h_pn_raquel[i]);

    }

    //make sure to handle the weights correcly
    for (int i = 0; i < h_list.size(); i++){
      h_list[i]->Sumw2();
     }
    for(int i = 0; i < h_list_2D.size(); i++){
     h_list_2D[i]->Sumw2();
     }
}

void twoproton_nuwro_overlay::Fill_Track_Plots(int which_cut, float value, int pdg, bool contained_start, bool contained_end,double wgt){

  //which cut identifies which cut we are looking at. 0 = after 3 pfps, 1 = after track score, 2 = after vertex distance
  h_track_overlay[0][which_cut][0]->Fill(trk_score_v->at(value),wgt); //fills the total
  h_track_overlay[1][which_cut][0]->Fill(trk_distance_v->at(value),wgt);
  h_track_overlay[2][which_cut][0]->Fill(trk_len_v->at(value),wgt);
  h_track_overlay[3][which_cut][0]->Fill(trk_llr_pid_score_v->at(value),wgt);  

  if(pdg == 2212 || pdg == -2212){
    total_protons++;
    if(contained_start == true && contained_end == true){
      h_track_overlay[0][which_cut][1]->Fill(trk_score_v->at(value),wgt); //fills the contained protons
      h_track_overlay[1][which_cut][1]->Fill(trk_distance_v->at(value),wgt);
      h_track_overlay[2][which_cut][1]->Fill(trk_len_v->at(value),wgt);
      h_track_overlay[3][which_cut][1]->Fill(trk_llr_pid_score_v->at(value),wgt);  
      contain++;
    } else if (contained_start == true && contained_end == false){
      h_track_overlay[0][which_cut][2]->Fill(trk_score_v->at(value),wgt); //fills the uncontained protons                                                                                                                 
      h_track_overlay[1][which_cut][2]->Fill(trk_distance_v->at(value),wgt);
      h_track_overlay[2][which_cut][2]->Fill(trk_len_v->at(value),wgt);
      h_track_overlay[3][which_cut][2]->Fill(trk_llr_pid_score_v->at(value),wgt);
      uncontain++;
    }

  } else if(pdg == 13 || pdg == -13){
    if(contained_start == true && contained_end == true){ //contained muons
      h_track_overlay[0][which_cut][3]->Fill(trk_score_v->at(value),wgt); 
      h_track_overlay[1][which_cut][3]->Fill(trk_distance_v->at(value),wgt);
      h_track_overlay[2][which_cut][3]->Fill(trk_len_v->at(value),wgt);
      h_track_overlay[3][which_cut][3]->Fill(trk_llr_pid_score_v->at(value),wgt);  

    }else if(contained_start == true && contained_end == false){ //uncontained muons
      h_track_overlay[0][which_cut][4]->Fill(trk_score_v->at(value),wgt);
      h_track_overlay[1][which_cut][4]->Fill(trk_distance_v->at(value),wgt);
      h_track_overlay[2][which_cut][4]->Fill(trk_len_v->at(value),wgt);
      h_track_overlay[3][which_cut][4]->Fill(trk_llr_pid_score_v->at(value),wgt);
    }

  } else if(pdg == 211 || pdg == -211) {
      h_track_overlay[0][which_cut][5]->Fill(trk_score_v->at(value),wgt); //fills the pionpm
      h_track_overlay[1][which_cut][5]->Fill(trk_distance_v->at(value),wgt);
      h_track_overlay[2][which_cut][5]->Fill(trk_len_v->at(value),wgt);
      h_track_overlay[3][which_cut][5]->Fill(trk_llr_pid_score_v->at(value),wgt);  

  } else if(pdg == 111) {
    h_track_overlay[0][which_cut][6]->Fill(trk_score_v->at(value),wgt); //fills the pion0
    h_track_overlay[1][which_cut][6]->Fill(trk_distance_v->at(value),wgt);
    h_track_overlay[2][which_cut][6]->Fill(trk_len_v->at(value),wgt);
    h_track_overlay[3][which_cut][6]->Fill(trk_llr_pid_score_v->at(value),wgt);  

  } else if(pdg == 11 || pdg == -11){
    h_track_overlay[0][which_cut][7]->Fill(trk_score_v->at(value),wgt); //fills the electron
    h_track_overlay[1][which_cut][7]->Fill(trk_distance_v->at(value),wgt);
    h_track_overlay[2][which_cut][7]->Fill(trk_len_v->at(value),wgt);
    h_track_overlay[3][which_cut][7]->Fill(trk_llr_pid_score_v->at(value),wgt);  

  } else if(pdg == 22){
    h_track_overlay[0][which_cut][8]->Fill(trk_score_v->at(value),wgt); //fills the gamma
    h_track_overlay[1][which_cut][8]->Fill(trk_distance_v->at(value),wgt);
    h_track_overlay[2][which_cut][8]->Fill(trk_len_v->at(value),wgt);
    h_track_overlay[3][which_cut][8]->Fill(trk_llr_pid_score_v->at(value),wgt);  

  } else if(pdg == 321 || pdg == -321 || pdg == 311){
    h_track_overlay[0][which_cut][9]->Fill(trk_score_v->at(value),wgt); //fills the kaon
    h_track_overlay[1][which_cut][9]->Fill(trk_distance_v->at(value),wgt);
    h_track_overlay[2][which_cut][9]->Fill(trk_len_v->at(value),wgt);
    h_track_overlay[3][which_cut][9]->Fill(trk_llr_pid_score_v->at(value),wgt);  

  } else {
    if(_debug) std::cout<<"[FILL_TRACK_PLOTS] Here is the Value of the PDG in the Else Loop: "<<pdg<<std::endl;
    other_else++;
    if(pdg == 2112){
      neutron++;
    }
    if( pdg == 14 || pdg == -14){
      neutrino++;
    }
    if(pdg == 0){
      zeros++;
    }
    h_track_overlay[0][which_cut][10]->Fill(trk_score_v->at(value),wgt); //fills the else
    h_track_overlay[1][which_cut][10]->Fill(trk_distance_v->at(value),wgt);
    h_track_overlay[2][which_cut][10]->Fill(trk_len_v->at(value),wgt);
    h_track_overlay[3][which_cut][10]->Fill(trk_llr_pid_score_v->at(value),wgt);  
  }
}

void twoproton_nuwro_overlay::Fill_Mine(int i, int j, double wgt){
  //index i indicates at which point the histograms are being filled 
  //index j represents what channel we are filling                                                                
  h_vtx_x_overlay[i][j]->Fill(reco_nu_vtx_sce_x,wgt);
  h_vtx_y_overlay[i][j]->Fill(reco_nu_vtx_sce_y,wgt);
  h_vtx_z_overlay[i][j]->Fill(reco_nu_vtx_sce_z,wgt);

  h_vtx_x_mc[i][j]->Fill(true_nu_vtx_x,wgt);
  h_vtx_y_mc[i][j]->Fill(true_nu_vtx_y,wgt);
  h_vtx_z_mc[i][j]->Fill(true_nu_vtx_z,wgt);

  h_vtx_x_mc_sce[i][j]->Fill(true_nu_vtx_sce_x,wgt);
  h_vtx_y_mc_sce[i][j]->Fill(true_nu_vtx_sce_y,wgt);
  h_vtx_z_mc_sce[i][j]->Fill(true_nu_vtx_sce_z,wgt);

  h_topological_score_overlay[i][j]->Fill(topological_score,wgt); //remember to add POT weight
  h_cosmic_impact_parameter_overlay[i][j]->Fill(CosmicIP,wgt); //remember to add POT weight
  //h_q2[i][j]->Fill(mc_q2,wgt);
  //h_X[i][j]->Fill(mc_X,wgt);
  //h_Y[i][j]->Fill(mc_Y,wgt);
  //h_Pt[i][j]->Fill(mc_Pt,wgt);
 }

void twoproton_nuwro_overlay::Fill_Raquel(int i, int j, double wgt){
  //index i indicates at which point the histograms are being filled                                                 
  //index j represents what channel we are filling          
  h_vtx_x_raquel[i][j]->Fill(reco_nu_vtx_sce_x,wgt);
  h_vtx_y_raquel[i][j]->Fill(reco_nu_vtx_sce_y,wgt);
  h_vtx_z_raquel[i][j]->Fill(reco_nu_vtx_sce_z,wgt);

  h_vtx_x_mc_raquel[i][j]->Fill(true_nu_vtx_x,wgt);
  h_vtx_y_mc_raquel[i][j]->Fill(true_nu_vtx_y,wgt);
  h_vtx_z_mc_raquel[i][j]->Fill(true_nu_vtx_z,wgt);

  h_vtx_x_mc_sce_raquel[i][j]->Fill(true_nu_vtx_sce_x,wgt);
  h_vtx_y_mc_sce_raquel[i][j]->Fill(true_nu_vtx_sce_y,wgt);
  h_vtx_z_mc_sce_raquel[i][j]->Fill(true_nu_vtx_sce_z,wgt);

  h_topological_score_raquel[i][j]->Fill(topological_score,wgt); //remember to add POT weight
  h_cosmic_impact_parameter_raquel[i][j]->Fill(CosmicIP,wgt); //remember to add POT weight
  //h_q2_raquel[i][j]->Fill(mc_q2,wgt);
  ///h_X_raquel[i][j]->Fill(mc_X,wgt);
  //h_Y_raquel[i][j]->Fill(mc_Y,wgt);
  //h_Pt_raquel[i][j]->Fill(mc_Pt,wgt);
 }


void twoproton_nuwro_overlay::Fill_Histograms_Mine(int i, double wgt, int mc_n_threshold_muon, int mc_n_threshold_proton, int mc_n_threshold_pion0, int mc_n_threshold_pionpm, bool fv){
  Fill_Mine(i,0,wgt);
  //cc0p0pi                                                                                                                                  
  if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 0 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    Fill_Mine(i,1,wgt);
    cc0p0pi[i]++;
    //cc1p0pi                                                                                                                                  
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 1 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    Fill_Mine(i,2,wgt);
    cc1p0pi[i]++;
    //cc2p0pi                                                                                                           
  } else if (ccnc == 0 && nu_pdg == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 2 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    Fill_Mine(i,3,wgt);
    cc2p0pi[i]++;
    //ccNp0pi                                                                                                                                  
  } else if (ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton > 2 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    Fill_Mine(i,4,wgt);
    ccNp0pi[i]++;
    //ccNp1pi                                                                                                                                   
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton >= 0 && (mc_n_threshold_pion0 == 1 || mc_n_threshold_pionpm == 1) && fv == true){
    Fill_Mine(i,5,wgt);
    ccNp1pi[i]++;
    //ccNpNpi                                                                                                                                   
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton >= 0 && (mc_n_threshold_pion0 > 1 || mc_n_threshold_pionpm > 1) && fv == true){
    Fill_Mine(i,6,wgt);
    ccNpNpi[i]++;
    //CC NUE                                                                                                                                   
  } else if(ccnc == 0 && abs(nu_pdg) == 12 && fv == true){
    Fill_Mine(i,7,wgt);
    ccnue[i]++;
  //OUT OF FV                                                                                                                                
  } else if(fv == false){                                                                                                                  
    Fill_Mine(i,8,wgt);
    outfv[i]++;
    //NC                                                                                                                                       
  } else if(ccnc == 1 && fv == true){
    Fill_Mine(i,9,wgt);
    nc[i]++;
    //else                                                                                                                                     
  } else{
    Fill_Mine(i,10,wgt);
    other[i]++;
  }
}

void twoproton_nuwro_overlay::Fill_Histograms_Raquel(int i, double wgt, bool fv){
  Fill_Raquel(i,0, wgt);
  //CCQE
  if(ccnc == 0 && interaction == 0 && abs(nu_pdg) == 14 && fv==true){
    Fill_Raquel(i,1, wgt);
    qel[i]++;
    //CCCoh                                                                                                                                                                                             
  } else if(ccnc == 0 && interaction == 3 && abs(nu_pdg) == 14 && fv == true){
    Fill_Raquel(i,2, wgt);
    coh[i]++;
    //CCMEC                                                                                                                                                                                             
  } else if(ccnc == 0 && interaction == 10 && abs(nu_pdg) == 14 && fv==true){
    Fill_Raquel(i,3, wgt);
    mec[i]++;

    //checking the mec weight
    if(std::isfinite(weightTune) && weightTune <= 100.) {
      mec_wgt = weightSplineTimesTune;
    } else {
      mec_wgt = 1 * weightSpline;
    }
    h_mec_wgt[i]->Fill(mec_wgt,1);

    //CCRES                                                                                                                                                                                             
  } else if(ccnc == 0 && interaction == 1 && abs(nu_pdg) == 14 && fv==true){
    //if(mc_n_threshold_muon == 1 && mc_n_threshold_proton == 2){
    //  res_count[0]++;
    //}else if (mc_n_threshold_muon == 1 && mc_n_threshold_proton == 1 && (mc_n_threshold_pion0 == 1 || mc_n_threshold_pionpm == 1)){
    //  res_count[1]++;
    //}else if(mc_n_threshold_muon == 1 && mc_n_threshold_proton > 2){
    //  res_count[2]++;
    //}else{
    //  res_count[3]++;
    //}
    Fill_Raquel(i,4,wgt);
    res[i]++;
    //CCDIS                                                                                                                                                                                             
  } else if(ccnc == 0 && interaction == 2 && abs(nu_pdg) == 14 && fv==true){
  Fill_Raquel(i,5,wgt);
    dis[i]++;
    //CCNue                                                                                                                                                                                             
  } else if(ccnc == 0  && abs(nu_pdg) == 12 && fv ==true){
  Fill_Raquel(i,6, wgt);
    ccnue_raquel[i]++;
    //NC                                                                                                                                                                                                
  } else if(ccnc == 1 && fv == true){
    Fill_Raquel(i,7, wgt);
    nc_raquel[i]++;
    //OUT OF FV                                                                                                                                                                                          
  } else if(fv == false){
    Fill_Raquel(i,8, wgt);
    outfv_raquel[i]++;
    //Other                                                                                                                                                                                             
  }else{
    Fill_Raquel(i,9, wgt);
    other_raquel[i]++;
  } 
}

void twoproton_nuwro_overlay::Fill_Histograms_Particles(int mc_n_threshold_muon, int mc_n_threshold_proton, int mc_n_threshold_pion0, double mc_n_threshold_pionpm, bool fv,TVector3 vMuon, TLorentzVector muon, TVector3 vLead, TLorentzVector lead, TVector3 vRec, TLorentzVector rec, double wgt){
  Fill_Particles(0,vMuon,muon,vLead,lead,vRec,rec,wgt);

  //cc0p0pi                                                                                                                               
  if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 0 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    Fill_Particles(1, vMuon,muon,vLead,lead,vRec,rec,wgt);
    cc0p0pi[number]++;

    //cc1p0pi                                                                                                                            
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 1 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    Fill_Particles(2, vMuon,muon,vLead,lead,vRec,rec,wgt);
    cc1p0pi[number]++;

    //cc2p0pi                                                                                   
  } else if (ccnc == 0 && nu_pdg == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 2 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    Fill_Particles(3, vMuon,muon,vLead,lead,vRec,rec,wgt);
    cc2p0pi[number]++;
    //ccNp0pi                                                                                                                            
  } else if (ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton > 2 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    Fill_Particles(4, vMuon,muon,vLead,lead,vRec,rec,wgt);
    ccNp0pi[number]++;
    //ccNp1pi                                                                                                                            
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton >= 0 && (mc_n_threshold_pion0 == 1 || mc_n_threshold_pionpm == 1) && fv == true){
    Fill_Particles(5, vMuon,muon,vLead,lead,vRec,rec,wgt);
    ccNp1pi[number]++;
    //ccNpNpi                                                                                                                            
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton >= 0 && (mc_n_threshold_pion0 > 1 || mc_n_threshold_pionpm > 1) && fv == true){
    Fill_Particles(6, vMuon,muon,vLead,lead,vRec,rec,wgt);
    ccNpNpi[number]++;
    //CC NUE                                                                                                                             
  } else if(ccnc == 0 && abs(nu_pdg) == 12 && fv == true){
    Fill_Particles(7, vMuon,muon,vLead,lead,vRec,rec,wgt);
    ccnue[number]++;
    //OUT OF FV                                                                                                                            
  } else if(fv == false){                                                                                                                     
    Fill_Particles(8, vMuon,muon,vLead,lead,vRec,rec,wgt);
    outfv[number]++;
    //NC                                                                                                                                 
  } else if(ccnc == 1 && fv == true){
    Fill_Particles(9, vMuon,muon,vLead,lead,vRec,rec,wgt);
    nc[number]++;    
    //else                                                                                                                               
  } else{
    Fill_Particles(10, vMuon,muon,vLead,lead,vRec,rec,wgt);
    other[number]++;
  }

}

void twoproton_nuwro_overlay::Fill_Histograms_Particles_Raquel(TVector3 vMuon, TLorentzVector muon, TVector3 vLead, TLorentzVector lead, TVector3 vRec, TLorentzVector rec, double wgt, bool fv){
  Fill_Particles_Raquel(0, vMuon,muon,vLead,lead,vRec,rec,wgt);

  //CCQE
  if(ccnc == 0 && interaction == 0 && fv==true){
    Fill_Particles_Raquel(1, vMuon,muon,vLead,lead,vRec,rec,wgt);
    qel[number]++;

    //CCCoh                                                                                                                                                                                             
  } else if(ccnc == 0 && interaction == 3 && fv == true){
    Fill_Particles_Raquel(2, vMuon,muon,vLead,lead,vRec,rec,wgt);
    coh[number]++;

    //CCMEC                                                                                                                                                                                             
  } else if(ccnc == 0 && interaction == 10 && fv==true){
    Fill_Particles_Raquel(3, vMuon,muon,vLead,lead,vRec,rec,wgt);
    mec[number]++;

    //CCRES                                                                                                                                                                                             
  } else if(ccnc == 0 && interaction == 1 && fv==true){
    Fill_Particles_Raquel(4, vMuon,muon,vLead,lead,vRec,rec,wgt);
    res[number]++;

    //CCDIS                                                                                                                                                                                             
  } else if(ccnc == 0 && interaction == 2 && fv==true){
    Fill_Particles_Raquel(5, vMuon,muon,vLead,lead,vRec,rec,wgt);
    dis[number]++;

    //CCNue                                                                                                                                                                                             
  } else if(ccnc == 0 && abs(nu_pdg) == 12 && fv ==true){
    Fill_Particles_Raquel(6, vMuon,muon,vLead,lead,vRec,rec,wgt);
    ccnue_raquel[number]++;

    //NC                                                                                                                                                                                                
  } else if(ccnc == 1 && fv == true){
    Fill_Particles_Raquel(7, vMuon,muon,vLead,lead,vRec,rec,wgt);
    nc_raquel[number]++;

    //OUT OF FV                                                                                                                                                                                          
  } else if(fv == false){
    Fill_Particles_Raquel(8, vMuon,muon,vLead,lead,vRec,rec,wgt);
    outfv_raquel[number]++;

    //Other                                                                                                                                                                                             
  }else{
    Fill_Particles_Raquel(9, vMuon,muon,vLead,lead,vRec,rec,wgt);
    other_raquel[number]++;
  } 
} //end of fill particles raquel


void twoproton_nuwro_overlay::Fill_Particles(int j, TVector3 vMuon, TLorentzVector muon, TVector3 vLead, TLorentzVector lead, TVector3 vRec, TLorentzVector rec, double wgt){

  //Going to use Calculate Variables function inside of variables.h. Returns Following 
  // 1) vector: momenta(muon_mom,lead_mom,rec_mom);                                  
  // 2) vector: Energies(KE_muon, TotE_muon, KE_Lead, TotE_Lead, KE_Rec, TotE_Rec);
  // 3) vector: detector_angles(muon_theta,muon_phi,lead_theta,lead_phi,recoil_theta,recoil_phi);
  // 4) vector: opening_angles(opening_angle_protons_lab,opening_angle_protons_mu_leading,opening_angle_protons_mu_both);  
  // 5) double: opening_angle_protons_COM 
  // 6) vector: STVS(delta_pT,delta_alphaT,delta_phiT);  
  // 7) double: calculated_nu_E                                                                                                                                                                                                                      
  
  variables.Calculate_Variables(vMuon,vLead,vRec,add_protons);

  double muon_mom = variables.momenta[0];
  double lead_mom = variables.momenta[1];
  double recoil_mom = variables.momenta[2];

  h_muon_overlay[0][j]->Fill(variables.momenta[0],wgt); //muon mom
  h_muon_overlay[1][j]->Fill(variables.Energies[0],wgt); //muon energy
  h_muon_overlay[2][j]->Fill(variables.detector_angles[0],wgt); //muon theta
  h_muon_overlay[3][j]->Fill(variables.detector_angles[1],wgt); //muon phi

  h_leading_overlay[0][j]->Fill(variables.momenta[1],wgt); //leading momentum
  h_leading_overlay[1][j]->Fill(variables.Energies[2],wgt); //leading energy
  h_leading_overlay[2][j]->Fill(variables.detector_angles[2],wgt); //leading theta
  h_leading_overlay[3][j]->Fill(variables.detector_angles[3],wgt); //leading phi

  h_recoil_overlay[0][j]->Fill(variables.momenta[2],wgt); //recoil momentum
  h_recoil_overlay[1][j]->Fill(variables.Energies[4],wgt); //recoil energy
  h_recoil_overlay[2][j]->Fill(variables.detector_angles[4],wgt); //recoil theta
  h_recoil_overlay[3][j]->Fill(variables.detector_angles[5],wgt); //recoil phi


  //Make sure to fill those resolution plots
  



  //Beam Stuff
  double EMuon =variables.Energies[0];
  double ELead = variables.Energies[2];
  double ERec = variables.Energies[4];
  double E_tot = (EMuon + MASS_MUON) + ELead + ERec;
  TVector3 PT_miss(vMuon[0]+vLead[0]+vRec[0],vMuon[1]+vRec[1]+vLead[1],0);
  double Eneutrino = variables.calculated_nu_E;//(EMuon+MASS_MUON) + ELead + ERec +((PT_miss.Mag2())/(2*35.37)) + 0.0304;
  TVector3 vBeam(0.,0.,Eneutrino); // z-direction is defined along the neutrino direction                               
  TVector3 vq = vBeam - vMuon; // Momentum transfer                                                                     
  TVector3 vmiss = vLead - vq; // Missing momentum         
  double E_tot_minus_beam = (E_tot - Eneutrino) * 1000;

  if(_debug) std::cout<<"[FILL_PARTICLES] Value of PT_miss magnitude: "<<PT_miss.Mag()<<std::endl;
  if(_debug) std::cout<<"[FILL_PARTICLES] Value of PT_miss magnitude2: "<<PT_miss.Mag2()<<std::endl;
  if(_debug) std::cout<<"[FILL_PARTICLES] Value of PT_miss magnitude2 divided by 2*35.37: "<<(PT_miss.Mag2())/(2.0*35.37)<<std::endl;
  if(_debug) std::cout<<"[FILL_PARTICLES] Value of Eneutrino: "<<Eneutrino<<std::endl;
  if(_debug) std::cout<<"[FILL_PARTICLES] Value of E_tot_minus_beam: "<<E_tot_minus_beam<<std::endl;

  //Struck nucleon Momentum:                                                                                            
  TVector3 vector_sum(vMuon[0] + vLead[0] + vRec[0], vMuon[1] + vLead[1] + vRec[1], vMuon[2] + vLead[2] + vRec[2]);
  TVector3 p_struck_nuc_vector(vector_sum[0], vector_sum[1], 0);
  double p_struck_nuc = p_struck_nuc_vector.Mag();
  double pz_tot = vLead[2] + vRec[2];

  //Energy of Struck Nucleon:
  double En = std::sqrt(std::pow(MASS_NEUTRON,2) + vmiss.Mag2()); //energy of struck nucleon 

  h_opening_angle_protons_lab_overlay[j]->Fill(variables.opening_angles[0],wgt);
  h_opening_angle_protons_com_overlay[j]->Fill(variables.opening_angle_protons_COM,wgt);  
  h_opening_angle_mu_leading_overlay[j]->Fill(variables.opening_angles[1],wgt);
  h_opening_angle_mu_both_overlay[j]->Fill(variables.opening_angles[2],wgt);
  h_delta_PT_overlay[j]->Fill(variables.stvs[0],wgt);
  h_delta_alphaT_overlay[j]->Fill(variables.stvs[1],wgt);
  h_delta_phiT_overlay[j]->Fill(variables.stvs[2],wgt);
  h_pn_overlay[j]->Fill(variables.stvs[3],wgt);
  h_nu_E_overlay[j]->Fill(Eneutrino,wgt);  
  h_mom_struck_nuc_overlay[j]->Fill(p_struck_nuc,wgt);
  h_tot_pz_overlay[j]->Fill(pz_tot,wgt);
  h_tot_E_overlay[j]->Fill(E_tot,wgt);
  h_tot_E_minus_beam_overlay[j]->Fill(E_tot_minus_beam,wgt);
  h_E_resolution_overlay[j]->Fill(Eneutrino - double(nu_e) ,wgt);
  h_PT_squared_overlay[j]->Fill(PT_miss.Mag2(),wgt);
  
  variables.momenta.clear();
  variables.detector_angles.clear();
  variables.opening_angles.clear();
  variables.stvs.clear();
  variables.Energies.clear();

}

void twoproton_nuwro_overlay::Fill_Particles_Raquel(int j, TVector3 vMuon, TLorentzVector muon, TVector3 vLead, TLorentzVector lead, TVector3 vRec, TLorentzVector rec, double wgt){

  //Going to use Calculate Variables function inside of variables.h. Returns Following 
  // 1) vector: momenta(muon_mom,lead_mom,rec_mom); 
  // 2) vector: Energies(KE_muon, TotE_muon, KE_Lead, TotE_Lead, KE_Rec, TotE_Rec);
  // 3) vector: detector_angles(muon_theta,muon_phi,lead_theta,lead_phi,recoil_theta,recoil_phi); 
  // 4) vector: opening_angles(opening_angle_protons_lab,opening_angle_protons_mu_leading,opening_angle_protons_mu_both);   
  // 5) double: opening_angle_protons_COM  
  // 6) vector: STVS(delta_pT,delta_alphaT,delta_phiT); 
  // 7) double: calculated_nu_E                                                                                                                                                                                                                      
  
  variables.Calculate_Variables(vMuon,vLead,vRec,add_protons);
  
  h_muon_raquel[0][j]->Fill(variables.momenta[0],wgt);
  h_leading_raquel[0][j]->Fill(variables.momenta[1],wgt);
  h_recoil_raquel[0][j]->Fill(variables.momenta[2],wgt);
  
  h_muon_raquel[1][j]->Fill(variables.Energies[0],wgt);
  h_leading_raquel[1][j]->Fill(variables.Energies[2],wgt);
  h_recoil_raquel[1][j]->Fill(variables.Energies[4],wgt);

  h_muon_raquel[2][j]->Fill(variables.detector_angles[0],wgt);
  h_leading_raquel[2][j]->Fill(variables.detector_angles[2],wgt);
  h_recoil_raquel[2][j]->Fill(variables.detector_angles[4],wgt);

  h_muon_raquel[3][j]->Fill(variables.detector_angles[1],wgt);
  h_leading_raquel[3][j]->Fill(variables.detector_angles[3],wgt);
  h_recoil_raquel[3][j]->Fill(variables.detector_angles[5],wgt);
  
  //Beam Stuff
  double EMuon =variables.Energies[0];
  double ELead = variables.Energies[2];
  double ERec = variables.Energies[4];
  double E_tot = (EMuon + MASS_MUON) + ELead + ERec;
  TVector3 PT_miss(vMuon[0]+vLead[0]+vRec[0],vMuon[1]+vRec[1]+vLead[1],0);
  double Eneutrino = variables.calculated_nu_E;//(EMuon+MASS_MUON) + ELead + ERec +((PT_miss.Mag2())/(2*35.37)) + 0.0304;
  TVector3 vBeam(0.,0.,Eneutrino); // z-direction is defined along the neutrino direction                               
  TVector3 vq = vBeam - vMuon; // Momentum transfer                                                                     
  TVector3 vmiss = vLead - vq; // Missing momentum         
  double E_tot_minus_beam = (E_tot - Eneutrino) * 1000;

  if(_debug) std::cout<<"[FILL_PARTICLES_RAQUEL] Value of PT_miss magnitude: "<<PT_miss.Mag()<<std::endl;
  if(_debug) std::cout<<"[FILL_PARTICLES_RAQUEL] Value of PT_miss magnitude2: "<<PT_miss.Mag2()<<std::endl;
  if(_debug) std::cout<<"[FILL_PARTICLES_RAQUEL] Value of PT_miss magnitude2 divided by 2*35.37: "<<(PT_miss.Mag2())/(2.0*35.37)<<std::endl;
  if(_debug) std::cout<<"[FILL_PARTICLES_RAQUEL] Value of Eneutrino: "<<Eneutrino<<std::endl;
  if(_debug) std::cout<<"[FILL_PARTICLES_RAQUEL] Value of E_tot_minus_beam: "<<E_tot_minus_beam<<std::endl;

  //Struck nucleon Momentum:                                                                                            
  TVector3 vector_sum(vMuon[0] + vLead[0] + vRec[0], vMuon[1] + vLead[1] + vRec[1], vMuon[2] + vLead[2] + vRec[2]);
  TVector3 p_struck_nuc_vector(vector_sum[0], vector_sum[1], 0);
  double p_struck_nuc = p_struck_nuc_vector.Mag();
  double pz_tot = vLead[2] + vRec[2];

  //Energy of Struck Nucleon:
  double En = std::sqrt(std::pow(MASS_NEUTRON,2) + vmiss.Mag2()); //energy of struck nucleon 

  h_opening_angle_protons_lab_raquel[j]->Fill(variables.opening_angles[0],wgt);
  h_opening_angle_protons_com_raquel[j]->Fill(variables.opening_angle_protons_COM,wgt);  
  h_opening_angle_mu_leading_raquel[j]->Fill(variables.opening_angles[1],wgt);
  h_opening_angle_mu_both_raquel[j]->Fill(variables.opening_angles[2],wgt);
  h_delta_PT_raquel[j]->Fill(variables.stvs[0],wgt);
  h_delta_alphaT_raquel[j]->Fill(variables.stvs[1],wgt);
  h_delta_phiT_raquel[j]->Fill(variables.stvs[2],wgt);
  h_pn_raquel[j]->Fill(variables.stvs[3],wgt);
  h_nu_E_raquel[j]->Fill(Eneutrino,wgt);
  h_mom_struck_nuc_raquel[j]->Fill(p_struck_nuc,wgt);
  h_tot_pz_raquel[j]->Fill(pz_tot,wgt);
  h_tot_E_raquel[j]->Fill(E_tot,wgt);
  h_tot_E_minus_beam_raquel[j]->Fill(E_tot_minus_beam,wgt);
  h_E_resolution_raquel[j]->Fill(Eneutrino - double(nu_e) ,wgt);
  h_PT_squared_raquel[j]->Fill(PT_miss.Mag2(),wgt);

  variables.momenta.clear();
  variables.detector_angles.clear();
  variables.opening_angles.clear();
  variables.stvs.clear();
  variables.Energies.clear();

}

void twoproton_nuwro_overlay::Write_Histograms(){
  for(int i=0; i< h_list_2D.size(); i++){
    h_list_2D[i]->Write();
  }
  for(int i = 0; i < h_list.size(); i++){
    h_list[i]->Write();
  }
}

Int_t twoproton_nuwro_overlay::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t twoproton_nuwro_overlay::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void twoproton_nuwro_overlay::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   dtrk_x_boundary = 0;
   dtrk_y_boundary = 0;
   dtrk_z_boundary = 0;
   dshr_x_boundary = 0;
   dshr_y_boundary = 0;
   dshr_z_boundary = 0;
   dvtx_x_boundary = 0;
   dvtx_y_boundary = 0;
   dvtx_z_boundary = 0;
   dtrk_boundary = 0;
   dvtx_boundary = 0;
   dshr_boundary = 0;
   dmc_boundary = 0;
   pfp_slice_idx = 0;
   backtracked_pdg = 0;
   backtracked_e = 0;
   backtracked_tid = 0;
   backtracked_purity = 0;
   backtracked_completeness = 0;
   backtracked_overlay_purity = 0;
   backtracked_px = 0;
   backtracked_py = 0;
   backtracked_pz = 0;
   backtracked_start_x = 0;
   backtracked_start_y = 0;
   backtracked_start_z = 0;
   backtracked_start_t = 0;
   backtracked_start_U = 0;
   backtracked_start_V = 0;
   backtracked_start_Y = 0;
   backtracked_sce_start_x = 0;
   backtracked_sce_start_y = 0;
   backtracked_sce_start_z = 0;
   backtracked_sce_start_U = 0;
   backtracked_sce_start_V = 0;
   backtracked_sce_start_Y = 0;
   pfp_generation_v = 0;
   pfp_trk_daughters_v = 0;
   pfp_shr_daughters_v = 0;
   trk_score_v = 0;
   pfpdg = 0;
   pfnhits = 0;
   pfnplanehits_U = 0;
   pfnplanehits_V = 0;
   pfnplanehits_Y = 0;
   pfpplanesubclusters_U = 0;
   pfpplanesubclusters_V = 0;
   pfpplanesubclusters_Y = 0;
   pfpplanesubhitfracmax_U = 0;
   pfpplanesubhitfracmax_V = 0;
   pfpplanesubhitfracmax_Y = 0;
   mc_pdg = 0;
   mc_E = 0;
   mc_vx = 0;
   mc_vy = 0;
   mc_vz = 0;
   mc_endx = 0;
   mc_endy = 0;
   mc_endz = 0;
   mc_px = 0;
   mc_py = 0;
   mc_pz = 0;
   mc_completeness = 0;
   mc_purity = 0;
   endmuonprocess = 0;
   weights = 0;
   weightsFlux = 0;
   weightsGenie = 0;
   weightsReint = 0;
   cosmic_flashmatch_score_v = 0;
   X_SpcPts_v = 0;
   Y_SpcPts_v = 0;
   Z_SpcPts_v = 0;
   trkshrscore_v = 0;
   PCAWin_1Cr_5cm = 0;
   PCAWin_2Cr_5cm = 0;
   PCAWin_3Cr_5cm = 0;
   PCAWin_dist_5cm = 0;
   PCAWin_npts_5cm = 0;
   PCAWin_1Cr_2_5cm = 0;
   PCAWin_2Cr_2_5cm = 0;
   PCAWin_3Cr_2_5cm = 0;
   PCAWin_dist_2_5cm = 0;
   PCAWin_npts_2_5cm = 0;
   shr_dedx_u_v = 0;
   shr_dedx_v_v = 0;
   shr_dedx_y_v = 0;
   shr_energy_u_v = 0;
   shr_energy_v_v = 0;
   shr_energy_y_v = 0;
   shr_pfp_id_v = 0;
   shr_start_x_v = 0;
   shr_start_y_v = 0;
   shr_start_z_v = 0;
   shr_dist_v = 0;
   shr_start_U_v = 0;
   shr_start_V_v = 0;
   shr_px_v = 0;
   shr_py_v = 0;
   shr_pz_v = 0;
   shr_openangle_v = 0;
   shr_theta_v = 0;
   shr_phi_v = 0;
   shr_pitch_u_v = 0;
   shr_pitch_v_v = 0;
   shr_pitch_y_v = 0;
   shr_tkfit_nhits_v = 0;
   shr_tkfit_start_x_v = 0;
   shr_tkfit_start_y_v = 0;
   shr_tkfit_start_z_v = 0;
   shr_tkfit_start_U_v = 0;
   shr_tkfit_start_V_v = 0;
   shr_tkfit_theta_v = 0;
   shr_tkfit_phi_v = 0;
   shr_tkfit_pitch_u_v = 0;
   shr_tkfit_pitch_v_v = 0;
   shr_tkfit_pitch_y_v = 0;
   shr_tkfit_dedx_u_v = 0;
   shr_tkfit_dedx_v_v = 0;
   shr_tkfit_dedx_y_v = 0;
   shr_tkfit_gap10_dedx_u_v = 0;
   shr_tkfit_gap10_dedx_v_v = 0;
   shr_tkfit_gap10_dedx_y_v = 0;
   shr_tkfit_dedx_nhits_u_v = 0;
   shr_tkfit_dedx_nhits_v_v = 0;
   shr_tkfit_dedx_nhits_y_v = 0;
   shr_llr_pid_u_v = 0;
   shr_llr_pid_v_v = 0;
   shr_llr_pid_y_v = 0;
   shr_llr_pid_v = 0;
   shr_llr_pid_score_v = 0;
   shr_moliere_avg_v = 0;
   shr_moliere_rms_v = 0;
   pfnunhits = 0;
   pflepnhits = 0;
   pfpronhits = 0;
   pfpi1nhits = 0;
   pfpi0nhits = 0;
   pfneunhits = 0;
   pfgamnhits = 0;
   pfothnhits = 0;
   trk_bragg_p_v = 0;
   trk_bragg_mu_v = 0;
   trk_bragg_mip_v = 0;
   trk_pida_v = 0;
   trk_pid_chipr_v = 0;
   trk_pid_chipi_v = 0;
   trk_pid_chika_v = 0;
   trk_pid_chimu_v = 0;
   trk_bragg_p_u_v = 0;
   trk_bragg_mu_u_v = 0;
   trk_bragg_mip_u_v = 0;
   trk_pida_u_v = 0;
   trk_pid_chipr_u_v = 0;
   trk_pid_chipi_u_v = 0;
   trk_pid_chika_u_v = 0;
   trk_pid_chimu_u_v = 0;
   trk_bragg_p_v_v = 0;
   trk_bragg_mu_v_v = 0;
   trk_bragg_mip_v_v = 0;
   trk_pida_v_v = 0;
   trk_pid_chipr_v_v = 0;
   trk_pid_chipi_v_v = 0;
   trk_pid_chika_v_v = 0;
   trk_pid_chimu_v_v = 0;
   trk_pfp_id_v = 0;
   trk_dir_x_v = 0;
   trk_dir_y_v = 0;
   trk_dir_z_v = 0;
   trk_start_x_v = 0;
   trk_start_y_v = 0;
   trk_start_z_v = 0;
   trk_sce_start_x_v = 0;
   trk_sce_start_y_v = 0;
   trk_sce_start_z_v = 0;
   trk_end_x_v = 0;
   trk_end_y_v = 0;
   trk_end_z_v = 0;
   trk_sce_end_x_v = 0;
   trk_sce_end_y_v = 0;
   trk_sce_end_z_v = 0;
   trk_distance_v = 0;
   trk_theta_v = 0;
   trk_phi_v = 0;
   trk_len_v = 0;
   trk_mcs_muon_mom_v = 0;
   trk_range_muon_mom_v = 0;
   trk_energy_proton_v = 0;
   trk_energy_muon_v = 0;
   trk_calo_energy_u_v = 0;
   trk_calo_energy_v_v = 0;
   trk_calo_energy_y_v = 0;
   trk_llr_pid_u_v = 0;
   trk_llr_pid_v_v = 0;
   trk_llr_pid_y_v = 0;
   trk_llr_pid_v = 0;
   trk_llr_pid_score_v = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("selected", &selected, &b_selected);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("sub", &sub, &b_sub);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("trk_id", &trk_id, &b_trk_pfp_id);
   fChain->SetBranchAddress("shr_id", &shr_id, &b_shr_pfp_id);
   fChain->SetBranchAddress("trk2_id", &trk2_id, &b_trk2_pfp_id);
   fChain->SetBranchAddress("shr2_id", &shr2_id, &b_shr2_pfp_id);
   fChain->SetBranchAddress("shr_energy_tot", &shr_energy_tot, &b_shr_energy_tot);
   fChain->SetBranchAddress("shr_energy", &shr_energy, &b_shr_energy);
   fChain->SetBranchAddress("shr_energy_tot_cali", &shr_energy_tot_cali, &b_shr_energy_tot_cali);
   fChain->SetBranchAddress("shr_energy_cali", &shr_energy_cali, &b_shr_energy_cali);
   fChain->SetBranchAddress("shr_theta", &shr_theta, &b_shr_theta);
   fChain->SetBranchAddress("shr_phi", &shr_phi, &b_shr_phi);
   fChain->SetBranchAddress("shr_pca_0", &shr_pca_0, &b_shr_pca_0);
   fChain->SetBranchAddress("shr_pca_1", &shr_pca_1, &b_shr_pca_1);
   fChain->SetBranchAddress("shr_pca_2", &shr_pca_2, &b_shr_pca_2);
   fChain->SetBranchAddress("shr_px", &shr_px, &b_shr_px);
   fChain->SetBranchAddress("shr_py", &shr_py, &b_shr_py);
   fChain->SetBranchAddress("shr_pz", &shr_pz, &b_shr_pz);
   fChain->SetBranchAddress("shr_openangle", &shr_openangle, &b_shr_openangle);
   fChain->SetBranchAddress("shr_tkfit_start_x", &shr_tkfit_start_x, &b_shr_tkfit_start_x);
   fChain->SetBranchAddress("shr_tkfit_start_y", &shr_tkfit_start_y, &b_shr_tkfit_start_y);
   fChain->SetBranchAddress("shr_tkfit_start_z", &shr_tkfit_start_z, &b_shr_tkfit_start_z);
   fChain->SetBranchAddress("shr_tkfit_theta", &shr_tkfit_theta, &b_shr_tkfit_theta);
   fChain->SetBranchAddress("shr_tkfit_phi", &shr_tkfit_phi, &b_shr_tkfit_phi);
   fChain->SetBranchAddress("shr_start_x", &shr_start_x, &b_shr_start_x);
   fChain->SetBranchAddress("shr_start_y", &shr_start_y, &b_shr_start_y);
   fChain->SetBranchAddress("shr_start_z", &shr_start_z, &b_shr_start_z);
   fChain->SetBranchAddress("shr_dedx_Y", &shr_dedx_Y, &b_shr_dedx_Y);
   fChain->SetBranchAddress("shr_dedx_V", &shr_dedx_V, &b_shr_dedx_V);
   fChain->SetBranchAddress("shr_dedx_U", &shr_dedx_U, &b_shr_dedx_U);
   fChain->SetBranchAddress("shr_dedx_Y_cali", &shr_dedx_Y_cali, &b_shr_dedx_Y_cali);
   fChain->SetBranchAddress("shr_dedx_V_cali", &shr_dedx_V_cali, &b_shr_dedx_V_cali);
   fChain->SetBranchAddress("shr_dedx_U_cali", &shr_dedx_U_cali, &b_shr_dedx_U_cali);
   fChain->SetBranchAddress("shr_tkfit_dedx_Y", &shr_tkfit_dedx_Y, &b_shr_tkfit_dedx_Y);
   fChain->SetBranchAddress("shr_tkfit_dedx_V", &shr_tkfit_dedx_V, &b_shr_tkfit_dedx_V);
   fChain->SetBranchAddress("shr_tkfit_dedx_U", &shr_tkfit_dedx_U, &b_shr_tkfit_dedx_U);
   fChain->SetBranchAddress("shr_tkfit_dedx_max", &shr_tkfit_dedx_max, &b_shr_tkfit_dedx_max);
   fChain->SetBranchAddress("shr_tkfit_nhits_Y", &shr_tkfit_nhits_Y, &b_shr_tkfit_nhits_Y);
   fChain->SetBranchAddress("shr_tkfit_nhits_V", &shr_tkfit_nhits_V, &b_shr_tkfit_nhits_V);
   fChain->SetBranchAddress("shr_tkfit_nhits_U", &shr_tkfit_nhits_U, &b_shr_tkfit_nhits_U);
   fChain->SetBranchAddress("shr_llrpid_dedx_Y", &shr_llrpid_dedx_Y, &b_shr_llrpid_dedx_Y);
   fChain->SetBranchAddress("shr_llrpid_dedx_V", &shr_llrpid_dedx_V, &b_shr_llrpid_dedx_V);
   fChain->SetBranchAddress("shr_llrpid_dedx_U", &shr_llrpid_dedx_U, &b_shr_llrpid_dedx_U);
   fChain->SetBranchAddress("shr_llrpid_dedx", &shr_llrpid_dedx, &b_shr_llrpid_dedx);
   fChain->SetBranchAddress("shr_tkfit_dedx_Y_alt", &shr_tkfit_dedx_Y_alt, &b_shr_tkfit_dedx_Y_alt);
   fChain->SetBranchAddress("shr_tkfit_dedx_V_alt", &shr_tkfit_dedx_V_alt, &b_shr_tkfit_dedx_V_alt);
   fChain->SetBranchAddress("shr_tkfit_dedx_U_alt", &shr_tkfit_dedx_U_alt, &b_shr_tkfit_dedx_U_alt);
   fChain->SetBranchAddress("shr_tkfit_nhits_Y_alt", &shr_tkfit_nhits_Y_alt, &b_shr_tkfit_nhits_Y_alt);
   fChain->SetBranchAddress("shr_tkfit_nhits_V_alt", &shr_tkfit_nhits_V_alt, &b_shr_tkfit_nhits_V_alt);
   fChain->SetBranchAddress("shr_tkfit_nhits_U_alt", &shr_tkfit_nhits_U_alt, &b_shr_tkfit_nhits_U_alt);
   fChain->SetBranchAddress("trkfit", &trkfit, &b__trkfit);
   fChain->SetBranchAddress("shr_tkfit_npoints", &shr_tkfit_npoints, &b_shr_tkfit_npoints);
   fChain->SetBranchAddress("shr_tkfit_npointsvalid", &shr_tkfit_npointsvalid, &b_shr_tkfit_npointsvalid);
   fChain->SetBranchAddress("shr_trkfitmedangle", &shr_trkfitmedangle, &b_shr_trkfitmedangle);
   fChain->SetBranchAddress("shrmoliereavg", &shrmoliereavg, &b_shrmoliereavg);
   fChain->SetBranchAddress("shrmoliererms", &shrmoliererms, &b_shrmoliererms);
   fChain->SetBranchAddress("shr1shr2moliereavg", &shr1shr2moliereavg, &b_shr1shr2moliereavg);
   fChain->SetBranchAddress("shr1shr2moliererms", &shr1shr2moliererms, &b_shr1shr2moliererms);
   fChain->SetBranchAddress("shr1trk1moliereavg", &shr1trk1moliereavg, &b_shr1trk1moliereavg);
   fChain->SetBranchAddress("shr1trk1moliererms", &shr1trk1moliererms, &b_shr1trk1moliererms);
   fChain->SetBranchAddress("shr1trk2moliereavg", &shr1trk2moliereavg, &b_shr1trk2moliereavg);
   fChain->SetBranchAddress("shr1trk2moliererms", &shr1trk2moliererms, &b_shr1trk2moliererms);
   fChain->SetBranchAddress("ismerged", &ismerged, &b_ismerged);
   fChain->SetBranchAddress("merge_bestdot", &merge_bestdot, &b_merge_bestdot);
   fChain->SetBranchAddress("merge_bestdist", &merge_bestdist, &b_merge_bestdist);
   fChain->SetBranchAddress("merge_vtx_x", &merge_vtx_x, &b_merge_vtx_x);
   fChain->SetBranchAddress("merge_vtx_y", &merge_vtx_y, &b_merge_vtx_y);
   fChain->SetBranchAddress("merge_vtx_z", &merge_vtx_z, &b_merge_vtx_z);
   fChain->SetBranchAddress("merge_tk_ipfp", &merge_tk_ipfp, &b_merge_tk_ipfp);
   fChain->SetBranchAddress("shr_tkfit_2cm_dedx_Y", &shr_tkfit_2cm_dedx_Y, &b_shr_tkfit_2cm_dedx_Y);
   fChain->SetBranchAddress("shr_tkfit_2cm_dedx_V", &shr_tkfit_2cm_dedx_V, &b_shr_tkfit_2cm_dedx_V);
   fChain->SetBranchAddress("shr_tkfit_2cm_dedx_U", &shr_tkfit_2cm_dedx_U, &b_shr_tkfit_2cm_dedx_U);
   fChain->SetBranchAddress("shr_tkfit_2cm_nhits_Y", &shr_tkfit_2cm_nhits_Y, &b_shr_tkfit_2cm_nhits_Y);
   fChain->SetBranchAddress("shr_tkfit_2cm_nhits_V", &shr_tkfit_2cm_nhits_V, &b_shr_tkfit_2cm_nhits_V);
   fChain->SetBranchAddress("shr_tkfit_2cm_nhits_U", &shr_tkfit_2cm_nhits_U, &b_shr_tkfit_2cm_nhits_U);
   fChain->SetBranchAddress("shr_tkfit_gap05_dedx_Y", &shr_tkfit_gap05_dedx_Y, &b_shr_tkfit_gap05_dedx_Y);
   fChain->SetBranchAddress("shr_tkfit_gap05_dedx_V", &shr_tkfit_gap05_dedx_V, &b_shr_tkfit_gap05_dedx_V);
   fChain->SetBranchAddress("shr_tkfit_gap05_dedx_U", &shr_tkfit_gap05_dedx_U, &b_shr_tkfit_gap05_dedx_U);
   fChain->SetBranchAddress("shr_tkfit_gap05_nhits_Y", &shr_tkfit_gap05_nhits_Y, &b_shr_tkfit_gap05_nhits_Y);
   fChain->SetBranchAddress("shr_tkfit_gap05_nhits_V", &shr_tkfit_gap05_nhits_V, &b_shr_tkfit_gap05_nhits_V);
   fChain->SetBranchAddress("shr_tkfit_gap05_nhits_U", &shr_tkfit_gap05_nhits_U, &b_shr_tkfit_gap05_nhits_U);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_Y", &shr_tkfit_gap10_dedx_Y, &b_shr_tkfit_gap10_dedx_Y);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_V", &shr_tkfit_gap10_dedx_V, &b_shr_tkfit_gap10_dedx_V);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_U", &shr_tkfit_gap10_dedx_U, &b_shr_tkfit_gap10_dedx_U);
   fChain->SetBranchAddress("shr_tkfit_gap10_nhits_Y", &shr_tkfit_gap10_nhits_Y, &b_shr_tkfit_gap10_nhits_Y);
   fChain->SetBranchAddress("shr_tkfit_gap10_nhits_V", &shr_tkfit_gap10_nhits_V, &b_shr_tkfit_gap10_nhits_V);
   fChain->SetBranchAddress("shr_tkfit_gap10_nhits_U", &shr_tkfit_gap10_nhits_U, &b_shr_tkfit_gap10_nhits_U);
   fChain->SetBranchAddress("shr_chipr", &shr_chipr, &b_shr_chipr);
   fChain->SetBranchAddress("shr_chimu", &shr_chimu, &b_shr_chimu);
   fChain->SetBranchAddress("shr_bragg_p", &shr_bragg_p, &b_shr_bragg_p);
   fChain->SetBranchAddress("shr_bragg_mu", &shr_bragg_mu, &b_shr_bragg_mu);
   fChain->SetBranchAddress("shr_bragg_mip", &shr_bragg_mip, &b_shr_bragg_mip);
   fChain->SetBranchAddress("shr_bragg_kaon", &shr_bragg_kaon, &b_shr_bragg_kaon);
   fChain->SetBranchAddress("shr_bragg_pion", &shr_bragg_pion, &b_shr_bragg_pion);
   fChain->SetBranchAddress("tksh_distance", &tksh_distance, &b_tksh_distance);
   fChain->SetBranchAddress("tksh_angle", &tksh_angle, &b_tksh_angle);
   fChain->SetBranchAddress("shr_distance", &shr_distance, &b_shr_distance);
   fChain->SetBranchAddress("shr_score", &shr_score, &b_shr_score);
   fChain->SetBranchAddress("shr_bkt_pdg", &shr_bkt_pdg, &b_shr_bkt_pdg);
   fChain->SetBranchAddress("shr_bkt_purity", &shr_bkt_purity, &b_shr_bkt_purity);
   fChain->SetBranchAddress("shr_bkt_completeness", &shr_bkt_completeness, &b_shr_bkt_completeness);
   fChain->SetBranchAddress("shr_bkt_E", &shr_bkt_E, &b_shr_bkt_E);
   fChain->SetBranchAddress("trk_len", &trk_len, &b_trk_len);
   fChain->SetBranchAddress("trk_theta", &trk_theta, &b_trk_theta);
   fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
   fChain->SetBranchAddress("trk_energy", &trk_energy, &b_trk_energy);
   fChain->SetBranchAddress("trk_energy_muon", &trk_energy_muon, &b_trk_energy_muon);
   fChain->SetBranchAddress("trk_energy_muon_mcs", &trk_energy_muon_mcs, &b_trk_energy_muon_mcs);
   fChain->SetBranchAddress("trk_energy_tot", &trk_energy_tot, &b_trk_energy_tot);
   fChain->SetBranchAddress("trk_energy_muon_tot", &trk_energy_muon_tot, &b_trk_energy_muon_tot);
   fChain->SetBranchAddress("trk_distance", &trk_distance, &b_trk_distance);
   fChain->SetBranchAddress("trk_score", &trk_score, &b_trk_score);
   fChain->SetBranchAddress("trk_bkt_pdg", &trk_bkt_pdg, &b_trk_bkt_pdg);
   fChain->SetBranchAddress("trk_bkt_purity", &trk_bkt_purity, &b_trk_bkt_purity);
   fChain->SetBranchAddress("trk_bkt_completeness", &trk_bkt_completeness, &b_trk_bkt_completeness);
   fChain->SetBranchAddress("trk_bkt_E", &trk_bkt_E, &b_trk_bkt_E);
   fChain->SetBranchAddress("trk_chipr_best", &trk_chipr_best, &b_trk_chipr_best);
   fChain->SetBranchAddress("trk_chipr_worst", &trk_chipr_worst, &b_trk_chipr_worst);
   fChain->SetBranchAddress("trk_chimu_best", &trk_chimu_best, &b_trk_chimu_best);
   fChain->SetBranchAddress("trk_chimu_worst", &trk_chimu_worst, &b_trk_chimu_worst);
   fChain->SetBranchAddress("trk_chipr", &trk_chipr, &b_trk_chipr);
   fChain->SetBranchAddress("trk_chimu", &trk_chimu, &b_trk_chimu);
   fChain->SetBranchAddress("trk_pida", &trk_pida, &b_trk_pida);
   fChain->SetBranchAddress("trk_bragg_p", &trk_bragg_p, &b_trk_bragg_p);
   fChain->SetBranchAddress("trk_bragg_mu", &trk_bragg_mu, &b_trk_bragg_mu);
   fChain->SetBranchAddress("trk_bragg_mip", &trk_bragg_mip, &b_trk_bragg_mip);
   fChain->SetBranchAddress("trk_bragg_kaon", &trk_bragg_kaon, &b_trk_bragg_kaon);
   fChain->SetBranchAddress("trk_bragg_pion", &trk_bragg_pion, &b_trk_bragg_pion);
   fChain->SetBranchAddress("trk_hits_max", &trk_hits_max, &b_trk_hits_max);
   fChain->SetBranchAddress("shr_hits_max", &shr_hits_max, &b_shr_hits_max);
   fChain->SetBranchAddress("trk_hits_2nd", &trk_hits_2nd, &b_trk_hits_2nd);
   fChain->SetBranchAddress("shr_hits_2nd", &shr_hits_2nd, &b_shr_hits_2nd);
   fChain->SetBranchAddress("trkshrhitdist0", &trkshrhitdist0, &b_trkshrhitdist0);
   fChain->SetBranchAddress("trkshrhitdist1", &trkshrhitdist1, &b_trkshrhitdist1);
   fChain->SetBranchAddress("trkshrhitdist2", &trkshrhitdist2, &b_trkshrhitdist2);
   fChain->SetBranchAddress("trk2shrhitdist0", &trk2shrhitdist0, &b_trk2shrhitdist0);
   fChain->SetBranchAddress("trk2shrhitdist1", &trk2shrhitdist1, &b_trk2shrhitdist1);
   fChain->SetBranchAddress("trk2shrhitdist2", &trk2shrhitdist2, &b_trk2shrhitdist2);
   fChain->SetBranchAddress("trk1trk2hitdist0", &trk1trk2hitdist0, &b_trk1trk2hitdist0);
   fChain->SetBranchAddress("trk1trk2hitdist1", &trk1trk2hitdist1, &b_trk1trk2hitdist1);
   fChain->SetBranchAddress("trk1trk2hitdist2", &trk1trk2hitdist2, &b_trk1trk2hitdist2);
   fChain->SetBranchAddress("total_hits_y", &total_hits_y, &b_total_hits_y);
   fChain->SetBranchAddress("extra_energy_y", &extra_energy_y, &b_extra_energy_y);
   fChain->SetBranchAddress("trk_energy_hits_tot", &trk_energy_hits_tot, &b_trk_energy_hits_tot);
   fChain->SetBranchAddress("subcluster", &subcluster, &b_subcluster);
   fChain->SetBranchAddress("shrsubclusters0", &shrsubclusters0, &b_shrsubclusters0);
   fChain->SetBranchAddress("shrsubclusters1", &shrsubclusters1, &b_shrsubclusters1);
   fChain->SetBranchAddress("shrsubclusters2", &shrsubclusters2, &b_shrsubclusters2);
   fChain->SetBranchAddress("shrclusfrac0", &shrclusfrac0, &b_shrclusfrac0);
   fChain->SetBranchAddress("shrclusfrac1", &shrclusfrac1, &b_shrclusfrac1);
   fChain->SetBranchAddress("shrclusfrac2", &shrclusfrac2, &b_shrclusfrac2);
   fChain->SetBranchAddress("shrclusdir0", &shrclusdir0, &b_shrclusdir0);
   fChain->SetBranchAddress("shrclusdir1", &shrclusdir1, &b_shrclusdir1);
   fChain->SetBranchAddress("shrclusdir2", &shrclusdir2, &b_shrclusdir2);
   fChain->SetBranchAddress("shr_hits_tot", &shr_hits_tot, &b_shr_hits_tot);
   fChain->SetBranchAddress("shr_hits_y_tot", &shr_hits_y_tot, &b_shr_hits_y_tot);
   fChain->SetBranchAddress("shr_hits_u_tot", &shr_hits_u_tot, &b_shr_hits_u_tot);
   fChain->SetBranchAddress("shr_hits_v_tot", &shr_hits_v_tot, &b_shr_hits_v_tot);
   fChain->SetBranchAddress("trk_hits_tot", &trk_hits_tot, &b_trk_hits_tot);
   fChain->SetBranchAddress("trk_hits_y_tot", &trk_hits_y_tot, &b_trk_hits_y_tot);
   fChain->SetBranchAddress("trk_hits_u_tot", &trk_hits_u_tot, &b_trk_hits_u_tot);
   fChain->SetBranchAddress("trk_hits_v_tot", &trk_hits_v_tot, &b_trk_hits_v_tot);
   fChain->SetBranchAddress("_elecclusters_U_charge", &_elecclusters_U_charge, &b_elecclusters_U_charge);
   fChain->SetBranchAddress("_elecclusters_V_charge", &_elecclusters_V_charge, &b_elecclusters_V_charge);
   fChain->SetBranchAddress("_elecclusters_Y_charge", &_elecclusters_Y_charge, &b_elecclusters_Y_charge);
   fChain->SetBranchAddress("_elecclusters_U_N", &_elecclusters_U_N, &b_elecclusters_U_N);
   fChain->SetBranchAddress("_elecclusters_V_N", &_elecclusters_V_N, &b_elecclusters_V_N);
   fChain->SetBranchAddress("_elecclusters_Y_N", &_elecclusters_Y_N, &b_elecclusters_Y_N);
   fChain->SetBranchAddress("n_tracks_contained", &n_tracks_contained, &b_n_tracks_contained);
   fChain->SetBranchAddress("n_showers_contained", &n_showers_contained, &b_n_showers_contained);
   fChain->SetBranchAddress("matched_E", &matched_E, &b_matched_E);
   fChain->SetBranchAddress("hits_ratio", &hits_ratio, &b_hits_ratio);
   fChain->SetBranchAddress("contained_fraction", &contained_fraction, &b_contained_fraction);
   fChain->SetBranchAddress("sps_contained_fraction", &sps_contained_fraction, &b_sps_contained_fraction);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("p", &p, &b_p);
   fChain->SetBranchAddress("pt_assume_muon", &pt_assume_muon, &b_pt_assume_muon);
   fChain->SetBranchAddress("p_assume_muon", &p_assume_muon, &b_p_assume_muon);
   fChain->SetBranchAddress("reco_e", &reco_e, &b_reco_e);
   fChain->SetBranchAddress("dvtx", &dvtx, &b_dvtx);
   fChain->SetBranchAddress("dtrk", &dtrk, &b_dtrk);
   fChain->SetBranchAddress("contained_sps_ratio", &contained_sps_ratio, &b_contained_sps_ratio);
   fChain->SetBranchAddress("dtrk_x_boundary", &dtrk_x_boundary, &b_dtrk_x_boundary);
   fChain->SetBranchAddress("dtrk_y_boundary", &dtrk_y_boundary, &b_dtrk_y_boundary);
   fChain->SetBranchAddress("dtrk_z_boundary", &dtrk_z_boundary, &b_dtrk_z_boundary);
   fChain->SetBranchAddress("dshr_x_boundary", &dshr_x_boundary, &b_dshr_x_boundary);
   fChain->SetBranchAddress("dshr_y_boundary", &dshr_y_boundary, &b_dshr_y_boundary);
   fChain->SetBranchAddress("dshr_z_boundary", &dshr_z_boundary, &b_dshr_z_boundary);
   fChain->SetBranchAddress("dvtx_x_boundary", &dvtx_x_boundary, &b_dvtx_x_boundary);
   fChain->SetBranchAddress("dvtx_y_boundary", &dvtx_y_boundary, &b_dvtx_y_boundary);
   fChain->SetBranchAddress("dvtx_z_boundary", &dvtx_z_boundary, &b_dvtx_z_boundary);
   fChain->SetBranchAddress("dtrk_boundary", &dtrk_boundary, &b_dtrk_boundary);
   fChain->SetBranchAddress("dvtx_boundary", &dvtx_boundary, &b_dvtx_boundary);
   fChain->SetBranchAddress("dshr_boundary", &dshr_boundary, &b_dshr_boundary);
   fChain->SetBranchAddress("dmc_boundary", &dmc_boundary, &b_dmc_boundary);
   fChain->SetBranchAddress("CosmicIP", &CosmicIP, &b_CosmicIP);
   fChain->SetBranchAddress("CosmicIPAll3D", &CosmicIPAll3D, &b_CosmicIPAll3D);
   fChain->SetBranchAddress("CosmicDirAll3D", &CosmicDirAll3D, &b_CosmicDirAll3D);
   fChain->SetBranchAddress("CosmicIPAll2DEnds", &CosmicIPAll2DEnds, &b_CosmicIPAll2DEnds);
   fChain->SetBranchAddress("CosmicDirAll2DEnds", &CosmicDirAll2DEnds, &b_CosmicDirAll2DEnds);
   fChain->SetBranchAddress("CosmicIPAll2DOvlp", &CosmicIPAll2DOvlp, &b_CosmicIPAll2DOvlp);
   fChain->SetBranchAddress("CosmicDirAll2DOvlp", &CosmicDirAll2DOvlp, &b_CosmicDirAll2DOvlp);
   fChain->SetBranchAddress("leeweight", &leeweight, &b_leeweight);
   fChain->SetBranchAddress("true_pt", &true_pt, &b_true_pt);
   fChain->SetBranchAddress("true_pt_visible", &true_pt_visible, &b_true_pt_visible);
   fChain->SetBranchAddress("true_p", &true_p, &b_true_p);
   fChain->SetBranchAddress("true_p_visible", &true_p_visible, &b_true_p_visible);
   fChain->SetBranchAddress("true_e_visible", &true_e_visible, &b_true_e_visible);
   fChain->SetBranchAddress("_opfilter_pe_beam", &_opfilter_pe_beam, &b_opfilter_pe_beam);
   fChain->SetBranchAddress("_opfilter_pe_veto", &_opfilter_pe_veto, &b_opfilter_pe_veto);
   fChain->SetBranchAddress("nu_pdg", &nu_pdg, &b_nu_pdg);
   fChain->SetBranchAddress("ccnc", &ccnc, &b_ccnc);
   fChain->SetBranchAddress("nu_parent_pdg", &nu_parent_pdg, &b_nu_parent_pdg);
   fChain->SetBranchAddress("nu_hadron_pdg", &nu_hadron_pdg, &b_nu_hadron_pdg);
   fChain->SetBranchAddress("nu_decay_mode", &nu_decay_mode, &b_nu_decay_mode);
   fChain->SetBranchAddress("interaction", &interaction, &b_interaction);
   fChain->SetBranchAddress("nu_e", &nu_e, &b_nu_e);
   fChain->SetBranchAddress("nu_l", &nu_l, &b_nu_l);
   fChain->SetBranchAddress("nu_pt", &nu_pt, &b_nu_pt);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("isVtxInFiducial", &isVtxInFiducial, &b_isVtxInFiducial);
   fChain->SetBranchAddress("truthFiducial", &truthFiducial, &b_truthFiducial);
   fChain->SetBranchAddress("true_nu_vtx_t", &true_nu_vtx_t, &b_true_nu_vtx_t);
   fChain->SetBranchAddress("true_nu_vtx_x", &true_nu_vtx_x, &b_true_nu_vtx_x);
   fChain->SetBranchAddress("true_nu_vtx_y", &true_nu_vtx_y, &b_true_nu_vtx_y);
   fChain->SetBranchAddress("true_nu_vtx_z", &true_nu_vtx_z, &b_true_nu_vtx_z);
   fChain->SetBranchAddress("true_nu_vtx_sce_x", &true_nu_vtx_sce_x, &b_true_nu_vtx_sce_x);
   fChain->SetBranchAddress("true_nu_vtx_sce_y", &true_nu_vtx_sce_y, &b_true_nu_vtx_sce_y);
   fChain->SetBranchAddress("true_nu_vtx_sce_z", &true_nu_vtx_sce_z, &b_true_nu_vtx_sce_z);
   fChain->SetBranchAddress("reco_nu_vtx_x", &reco_nu_vtx_x, &b_reco_nu_vtx_x);
   fChain->SetBranchAddress("reco_nu_vtx_y", &reco_nu_vtx_y, &b_reco_nu_vtx_y);
   fChain->SetBranchAddress("reco_nu_vtx_z", &reco_nu_vtx_z, &b_reco_nu_vtx_z);
   fChain->SetBranchAddress("reco_nu_vtx_sce_x", &reco_nu_vtx_sce_x, &b_reco_nu_vtx_sce_x);
   fChain->SetBranchAddress("reco_nu_vtx_sce_y", &reco_nu_vtx_sce_y, &b_reco_nu_vtx_sce_y);
   fChain->SetBranchAddress("reco_nu_vtx_sce_z", &reco_nu_vtx_sce_z, &b_reco_nu_vtx_sce_z);
   fChain->SetBranchAddress("nmuon", &nmuon, &b_nmuon);
   fChain->SetBranchAddress("muon_e", &muon_e, &b_muon_e);
   fChain->SetBranchAddress("muon_c", &muon_c, &b_muon_c);
   fChain->SetBranchAddress("muon_p", &muon_p, &b_muon_p);
   fChain->SetBranchAddress("nelec", &nelec, &b_nelec);
   fChain->SetBranchAddress("elec_e", &elec_e, &b_elec_e);
   fChain->SetBranchAddress("elec_c", &elec_c, &b_elec_c);
   fChain->SetBranchAddress("elec_p", &elec_p, &b_elec_p);
   fChain->SetBranchAddress("elec_vx", &elec_vx, &b_elec_vx);
   fChain->SetBranchAddress("elec_vy", &elec_vy, &b_elec_vy);
   fChain->SetBranchAddress("elec_vz", &elec_vz, &b_elec_vz);
   fChain->SetBranchAddress("elec_px", &elec_px, &b_elec_px);
   fChain->SetBranchAddress("elec_py", &elec_py, &b_elec_py);
   fChain->SetBranchAddress("elec_pz", &elec_pz, &b_elec_pz);
   fChain->SetBranchAddress("npi0", &npi0, &b_npi0);
   fChain->SetBranchAddress("pi0_e", &pi0_e, &b_pi0_e);
   fChain->SetBranchAddress("pi0_c", &pi0_c, &b_pi0_c);
   fChain->SetBranchAddress("pi0_p", &pi0_p, &b_pi0_p);
   fChain->SetBranchAddress("nneutron", &nneutron, &b_nneutron);
   fChain->SetBranchAddress("nproton", &nproton, &b_nproton);
   fChain->SetBranchAddress("proton_e", &proton_e, &b_proton_e);
   fChain->SetBranchAddress("proton_c", &proton_c, &b_proton_c);
   fChain->SetBranchAddress("proton_p", &proton_p, &b_proton_p);
   fChain->SetBranchAddress("npion", &npion, &b_npion);
   fChain->SetBranchAddress("pion_e", &pion_e, &b_pion_e);
   fChain->SetBranchAddress("pion_c", &pion_c, &b_pion_c);
   fChain->SetBranchAddress("pion_p", &pion_p, &b_pion_p);
   fChain->SetBranchAddress("neta", &neta, &b_neta);
   fChain->SetBranchAddress("eta_e", &eta_e, &b_eta_e);
   fChain->SetBranchAddress("nslice", &nslice, &b_nslice);
   fChain->SetBranchAddress("crtveto", &crtveto, &b_crtveto);
   fChain->SetBranchAddress("crthitpe", &crthitpe, &b_crthitpe);
   fChain->SetBranchAddress("pfp_slice_idx", &pfp_slice_idx, &b_pfp_slice_idx);
   fChain->SetBranchAddress("category", &category, &b_category);
   fChain->SetBranchAddress("backtracked_pdg", &backtracked_pdg, &b_backtracked_pdg);
   fChain->SetBranchAddress("backtracked_e", &backtracked_e, &b_backtracked_e);
   fChain->SetBranchAddress("backtracked_tid", &backtracked_tid, &b_backtracked_tid);
   fChain->SetBranchAddress("backtracked_purity", &backtracked_purity, &b_backtracked_purity);
   fChain->SetBranchAddress("backtracked_completeness", &backtracked_completeness, &b_backtracked_completeness);
   fChain->SetBranchAddress("backtracked_overlay_purity", &backtracked_overlay_purity, &b_backtracked_overlay_purity);
   fChain->SetBranchAddress("backtracked_px", &backtracked_px, &b_backtracked_px);
   fChain->SetBranchAddress("backtracked_py", &backtracked_py, &b_backtracked_py);
   fChain->SetBranchAddress("backtracked_pz", &backtracked_pz, &b_backtracked_pz);
   fChain->SetBranchAddress("backtracked_start_x", &backtracked_start_x, &b_backtracked_start_x);
   fChain->SetBranchAddress("backtracked_start_y", &backtracked_start_y, &b_backtracked_start_y);
   fChain->SetBranchAddress("backtracked_start_z", &backtracked_start_z, &b_backtracked_start_z);
   fChain->SetBranchAddress("backtracked_start_t", &backtracked_start_t, &b_backtracked_start_t);
   fChain->SetBranchAddress("backtracked_start_U", &backtracked_start_U, &b_backtracked_start_U);
   fChain->SetBranchAddress("backtracked_start_V", &backtracked_start_V, &b_backtracked_start_V);
   fChain->SetBranchAddress("backtracked_start_Y", &backtracked_start_Y, &b_backtracked_start_Y);
   fChain->SetBranchAddress("backtracked_sce_start_x", &backtracked_sce_start_x, &b_backtracked_sce_start_x);
   fChain->SetBranchAddress("backtracked_sce_start_y", &backtracked_sce_start_y, &b_backtracked_sce_start_y);
   fChain->SetBranchAddress("backtracked_sce_start_z", &backtracked_sce_start_z, &b_backtracked_sce_start_z);
   fChain->SetBranchAddress("backtracked_sce_start_U", &backtracked_sce_start_U, &b_backtracked_sce_start_U);
   fChain->SetBranchAddress("backtracked_sce_start_V", &backtracked_sce_start_V, &b_backtracked_sce_start_V);
   fChain->SetBranchAddress("backtracked_sce_start_Y", &backtracked_sce_start_Y, &b_backtracked_sce_start_Y);
   fChain->SetBranchAddress("lep_e", &lep_e, &b_lep_e);
   fChain->SetBranchAddress("pass", &pass, &b_pass);
   fChain->SetBranchAddress("swtrig", &swtrig, &b_swtrig);
   fChain->SetBranchAddress("evnhits", &evnhits, &b_evnhits);
   fChain->SetBranchAddress("slpdg", &slpdg, &b_slpdg);
   fChain->SetBranchAddress("slnhits", &slnhits, &b_slnhits);
   fChain->SetBranchAddress("n_pfps", &n_pfps, &b_n_pfps);
   fChain->SetBranchAddress("n_tracks", &n_tracks, &b_n_tracks);
   fChain->SetBranchAddress("n_showers", &n_showers, &b_n_showers);
   fChain->SetBranchAddress("pfp_generation_v", &pfp_generation_v, &b_pfp_generation_v);
   fChain->SetBranchAddress("pfp_trk_daughters_v", &pfp_trk_daughters_v, &b_pfp_trk_daughters_v);
   fChain->SetBranchAddress("pfp_shr_daughters_v", &pfp_shr_daughters_v, &b_pfp_shr_daughters_v);
   fChain->SetBranchAddress("trk_score_v", &trk_score_v, &b_trk_score_v);
   fChain->SetBranchAddress("pfpdg", &pfpdg, &b_pfpdg);
   fChain->SetBranchAddress("pfnhits", &pfnhits, &b_pfnhits);
   fChain->SetBranchAddress("pfnplanehits_U", &pfnplanehits_U, &b_pfnplanehits_U);
   fChain->SetBranchAddress("pfnplanehits_V", &pfnplanehits_V, &b_pfnplanehits_V);
   fChain->SetBranchAddress("pfnplanehits_Y", &pfnplanehits_Y, &b_pfnplanehits_Y);
   fChain->SetBranchAddress("pfpplanesubclusters_U", &pfpplanesubclusters_U, &b_pfpplanesubclusters_U);
   fChain->SetBranchAddress("pfpplanesubclusters_V", &pfpplanesubclusters_V, &b_pfpplanesubclusters_V);
   fChain->SetBranchAddress("pfpplanesubclusters_Y", &pfpplanesubclusters_Y, &b_pfpplanesubclusters_Y);
   fChain->SetBranchAddress("pfpplanesubhitfracmax_U", &pfpplanesubhitfracmax_U, &b_pfpplanesubhitfracmax_U);
   fChain->SetBranchAddress("pfpplanesubhitfracmax_V", &pfpplanesubhitfracmax_V, &b_pfpplanesubhitfracmax_V);
   fChain->SetBranchAddress("pfpplanesubhitfracmax_Y", &pfpplanesubhitfracmax_Y, &b_pfpplanesubhitfracmax_Y);
   fChain->SetBranchAddress("hits_u", &hits_u, &b_hits_u);
   fChain->SetBranchAddress("hits_v", &hits_v, &b_hits_v);
   fChain->SetBranchAddress("hits_y", &hits_y, &b_hits_y);
   fChain->SetBranchAddress("topological_score", &topological_score, &b_topological_score);
   fChain->SetBranchAddress("slclustfrac", &slclustfrac, &b_slclustfrac);
   fChain->SetBranchAddress("mc_pdg", &mc_pdg, &b_mc_pdg);
   fChain->SetBranchAddress("mc_E", &mc_E, &b_mc_E);
   fChain->SetBranchAddress("mc_vx", &mc_vx, &b_mc_vx);
   fChain->SetBranchAddress("mc_vy", &mc_vy, &b_mc_vy);
   fChain->SetBranchAddress("mc_vz", &mc_vz, &b_mc_vz);
   fChain->SetBranchAddress("mc_endx", &mc_endx, &b_mc_endx);
   fChain->SetBranchAddress("mc_endy", &mc_endy, &b_mc_endy);
   fChain->SetBranchAddress("mc_endz", &mc_endz, &b_mc_endz);
   fChain->SetBranchAddress("mc_px", &mc_px, &b_mc_px);
   fChain->SetBranchAddress("mc_py", &mc_py, &b_mc_py);
   fChain->SetBranchAddress("mc_pz", &mc_pz, &b_mc_pz);
   fChain->SetBranchAddress("mc_completeness", &mc_completeness, &b_mc_completeness);
   fChain->SetBranchAddress("mc_purity", &mc_purity, &b_mc_purity);
   fChain->SetBranchAddress("endmuonprocess", &endmuonprocess, &b_endmuonprocess);
   fChain->SetBranchAddress("endmuonmichel", &endmuonmichel, &b_endmuonmichel);
   fChain->SetBranchAddress("filter_antibdt", &filter_antibdt, &b_filter_antibdt);
   fChain->SetBranchAddress("filter_ncpi0", &filter_ncpi0, &b_filter_ncpi0);
   fChain->SetBranchAddress("filter_pi0", &filter_pi0, &b_filter_pi0);
   fChain->SetBranchAddress("filter_ccinclusive", &filter_ccinclusive, &b_filter_ccinclusive);
   fChain->SetBranchAddress("weights", &weights, &b_weights);
   fChain->SetBranchAddress("weightsFlux", &weightsFlux, &b_weightsFlux);
   fChain->SetBranchAddress("weightsGenie", &weightsGenie, &b_weightsGenie);
   fChain->SetBranchAddress("weightsReint", &weightsReint, &b_weightsReint);
   fChain->SetBranchAddress("weightSpline", &weightSpline, &b_weightSpline);
   fChain->SetBranchAddress("weightTune", &weightTune, &b_weightTune);
   fChain->SetBranchAddress("weightSplineTimesTune", &weightSplineTimesTune, &b_weightSplineTimesTune);
   fChain->SetBranchAddress("knobRPAup", &knobRPAup, &b_knobRPAup);
   fChain->SetBranchAddress("knobRPAdn", &knobRPAdn, &b_knobRPAdn);
   fChain->SetBranchAddress("knobCCMECup", &knobCCMECup, &b_knobCCMECup);
   fChain->SetBranchAddress("knobCCMECdn", &knobCCMECdn, &b_knobCCMECdn);
   fChain->SetBranchAddress("knobAxFFCCQEup", &knobAxFFCCQEup, &b_knobAxFFCCQEup);
   fChain->SetBranchAddress("knobAxFFCCQEdn", &knobAxFFCCQEdn, &b_knobAxFFCCQEdn);
   fChain->SetBranchAddress("knobVecFFCCQEup", &knobVecFFCCQEup, &b_knobVecFFCCQEup);
   fChain->SetBranchAddress("knobVecFFCCQEdn", &knobVecFFCCQEdn, &b_knobVecFFCCQEdn);
   fChain->SetBranchAddress("knobDecayAngMECup", &knobDecayAngMECup, &b_knobDecayAngMECup);
   fChain->SetBranchAddress("knobDecayAngMECdn", &knobDecayAngMECdn, &b_knobDecayAngMECdn);
   fChain->SetBranchAddress("knobThetaDelta2Npiup", &knobThetaDelta2Npiup, &b_knobThetaDelta2Npiup);
   fChain->SetBranchAddress("knobThetaDelta2Npidn", &knobThetaDelta2Npidn, &b_knobThetaDelta2Npidn);
   fChain->SetBranchAddress("flash_pe", &flash_pe, &b_flash_pe);
   fChain->SetBranchAddress("flash_time", &flash_time, &b_flash_time);
   fChain->SetBranchAddress("nu_flashmatch_score", &nu_flashmatch_score, &b_nu_flashmatch_score);
   fChain->SetBranchAddress("best_cosmic_flashmatch_score", &best_cosmic_flashmatch_score, &b_best_cosmic_flashmatch_score);
   fChain->SetBranchAddress("best_obviouscosmic_flashmatch_score", &best_obviouscosmic_flashmatch_score, &b_best_obviouscosmic_flashmatch_score);
   fChain->SetBranchAddress("cosmic_flashmatch_score_v", &cosmic_flashmatch_score_v, &b_cosmic_flashmatch_score_v);
   fChain->SetBranchAddress("mcf_nu_e", &mcf_nu_e, &b_mcf_nu_e);
   fChain->SetBranchAddress("mcf_lep_e", &mcf_lep_e, &b_mcf_lep_e);
   fChain->SetBranchAddress("mcf_actvol", &mcf_actvol, &b_mcf_actvol);
   fChain->SetBranchAddress("mcf_nmm", &mcf_nmm, &b_mcf_nmm);
   fChain->SetBranchAddress("mcf_nmp", &mcf_nmp, &b_mcf_nmp);
   fChain->SetBranchAddress("mcf_nem", &mcf_nem, &b_mcf_nem);
   fChain->SetBranchAddress("mcf_nep", &mcf_nep, &b_mcf_nep);
   fChain->SetBranchAddress("mcf_np0", &mcf_np0, &b_mcf_np0);
   fChain->SetBranchAddress("mcf_npp", &mcf_npp, &b_mcf_npp);
   fChain->SetBranchAddress("mcf_npm", &mcf_npm, &b_mcf_npm);
   fChain->SetBranchAddress("mcf_mcshr_elec_etot", &mcf_mcshr_elec_etot, &b_mcf_mcshr_elec_etot);
   fChain->SetBranchAddress("mcf_pass_ccpi0", &mcf_pass_ccpi0, &b_mcf_pass_ccpi0);
   fChain->SetBranchAddress("mcf_pass_ncpi0", &mcf_pass_ncpi0, &b_mcf_pass_ncpi0);
   fChain->SetBranchAddress("mcf_pass_ccnopi", &mcf_pass_ccnopi, &b_mcf_pass_ccnopi);
   fChain->SetBranchAddress("mcf_pass_ncnopi", &mcf_pass_ncnopi, &b_mcf_pass_ncnopi);
   fChain->SetBranchAddress("mcf_pass_cccpi", &mcf_pass_cccpi, &b_mcf_pass_cccpi);
   fChain->SetBranchAddress("mcf_pass_nccpi", &mcf_pass_nccpi, &b_mcf_pass_nccpi);
   fChain->SetBranchAddress("X_SpcPts_v", &X_SpcPts_v, &b_X_SpcPts_v);
   fChain->SetBranchAddress("Y_SpcPts_v", &Y_SpcPts_v, &b_Y_SpcPts_v);
   fChain->SetBranchAddress("Z_SpcPts_v", &Z_SpcPts_v, &b_Z_SpcPts_v);
   fChain->SetBranchAddress("shr_id_MCStool", &shr_id_MCStool, &b_shr_pfp_id);
   fChain->SetBranchAddress("shr_hits_max_MCStool", &shr_hits_max_MCStool, &b_shr_hits_max_MCStool);
   fChain->SetBranchAddress("n_showers_contained_MCStool", &n_showers_contained_MCStool, &b_n_showers_contained_MCStool);
   fChain->SetBranchAddress("trkshrscore_v", &trkshrscore_v, &b_trkshrscore_v);
   fChain->SetBranchAddress("shrPCA_1Cr", &shrPCA_1Cr, &b_shrPCA_1Cr);
   fChain->SetBranchAddress("shrPCA_2Cr", &shrPCA_2Cr, &b_shrPCA_2Cr);
   fChain->SetBranchAddress("shrPCA_3Cr", &shrPCA_3Cr, &b_shrPCA_3Cr);
   fChain->SetBranchAddress("shrPCA_1Ce", &shrPCA_1Ce, &b_shrPCA_1Ce);
   fChain->SetBranchAddress("shrPCA_2Ce", &shrPCA_2Ce, &b_shrPCA_2Ce);
   fChain->SetBranchAddress("shrPCA_3Ce", &shrPCA_3Ce, &b_shrPCA_3Ce);
   fChain->SetBranchAddress("shrPCA1CAS", &shrPCA1CAS, &b_shrPCA1CAS);
   fChain->SetBranchAddress("shrPCA2CAS", &shrPCA2CAS, &b_shrPCA2CAS);
   fChain->SetBranchAddress("shrPCA3CAS", &shrPCA3CAS, &b_shrPCA3CAS);
   fChain->SetBranchAddress("shrPCA_1Cr2h", &shrPCA_1Cr2h, &b_shrPCA_1Cr2h);
   fChain->SetBranchAddress("shrPCA_2Cr2h", &shrPCA_2Cr2h, &b_shrPCA_2Cr2h);
   fChain->SetBranchAddress("shrPCA_3Cr2h", &shrPCA_3Cr2h, &b_shrPCA_3Cr2h);
   fChain->SetBranchAddress("shrPCA_1Cr1h", &shrPCA_1Cr1h, &b_shrPCA_1Cr1h);
   fChain->SetBranchAddress("shrPCA_2Cr1h", &shrPCA_2Cr1h, &b_shrPCA_2Cr1h);
   fChain->SetBranchAddress("shrPCA_3Cr1h", &shrPCA_3Cr1h, &b_shrPCA_3Cr1h);
   fChain->SetBranchAddress("shrMCSMom", &shrMCSMom, &b_shrMCSMom);
   fChain->SetBranchAddress("shrMCSMom1h", &shrMCSMom1h, &b_shrMCSMom1h);
   fChain->SetBranchAddress("shrMCSMom2h", &shrMCSMom2h, &b_shrMCSMom2h);
   fChain->SetBranchAddress("shrPCALen", &shrPCALen, &b_shrPCALen);
   fChain->SetBranchAddress("n_shrSpcPts", &n_shrSpcPts, &b_n_shrSpcPts);
   fChain->SetBranchAddress("PCAWin_1Cr_5cm", &PCAWin_1Cr_5cm, &b_PCAWin_1Cr_5cm);
   fChain->SetBranchAddress("PCAWin_2Cr_5cm", &PCAWin_2Cr_5cm, &b_PCAWin_2Cr_5cm);
   fChain->SetBranchAddress("PCAWin_3Cr_5cm", &PCAWin_3Cr_5cm, &b_PCAWin_3Cr_5cm);
   fChain->SetBranchAddress("PCAWin_dist_5cm", &PCAWin_dist_5cm, &b_PCAWin_dist_5cm);
   fChain->SetBranchAddress("PCAWin_npts_5cm", &PCAWin_npts_5cm, &b_PCAWin_npts_5cm);
   fChain->SetBranchAddress("shrStart_5cm", &shrStart_5cm, &b_shrStart_5cm);
   fChain->SetBranchAddress("shrStartMCS_5cm", &shrStartMCS_5cm, &b_shrStartMCS_5cm);
   fChain->SetBranchAddress("shrMCSAS_5cm", &shrMCSAS_5cm, &b_shrMCSAS_5cm);
   fChain->SetBranchAddress("shrPCA1CAS_5cm", &shrPCA1CAS_5cm, &b_shrPCA1CAS_5cm);
   fChain->SetBranchAddress("shrPCA2CAS_5cm", &shrPCA2CAS_5cm, &b_shrPCA2CAS_5cm);
   fChain->SetBranchAddress("shrPCA3CAS_5cm", &shrPCA3CAS_5cm, &b_shrPCA3CAS_5cm);
   fChain->SetBranchAddress("shrPCA1CMed_5cm", &shrPCA1CMed_5cm, &b__shrPCA1CMed_5cm);
   fChain->SetBranchAddress("PCAWin_1Cr_2_5cm", &PCAWin_1Cr_2_5cm, &b_PCAWin_1Cr_2_5cm);
   fChain->SetBranchAddress("PCAWin_2Cr_2_5cm", &PCAWin_2Cr_2_5cm, &b_PCAWin_2Cr_2_5cm);
   fChain->SetBranchAddress("PCAWin_3Cr_2_5cm", &PCAWin_3Cr_2_5cm, &b_PCAWin_3Cr_2_5cm);
   fChain->SetBranchAddress("PCAWin_dist_2_5cm", &PCAWin_dist_2_5cm, &b_PCAWin_dist_2_5cm);
   fChain->SetBranchAddress("PCAWin_npts_2_5cm", &PCAWin_npts_2_5cm, &b_PCAWin_npts_2_5cm);
   fChain->SetBranchAddress("shrStart_2_5cm", &shrStart_2_5cm, &b_shrStart_2_5cm);
   fChain->SetBranchAddress("shrStartMCS_2_5cm", &shrStartMCS_2_5cm, &b_shrStartMCS_2_5cm);
   fChain->SetBranchAddress("shrMCSAS_2_5cm", &shrMCSAS_2_5cm, &b_shrMCSAS_2_5cm);
   fChain->SetBranchAddress("shrPCA1CAS_2_5cm", &shrPCA1CAS_2_5cm, &b_shrPCA1CAS_2_5cm);
   fChain->SetBranchAddress("shrPCA2CAS_2_5cm", &shrPCA2CAS_2_5cm, &b_shrPCA2CAS_2_5cm);
   fChain->SetBranchAddress("shrPCA3CAS_2_5cm", &shrPCA3CAS_2_5cm, &b_shrPCA3CAS_2_5cm);
   fChain->SetBranchAddress("shrPCA1CMed_2_5cm", &shrPCA1CMed_2_5cm, &b__shrPCA1CMed_2_5cm);
   fChain->SetBranchAddress("DeltaMed", &DeltaMed, &b_DeltaMed);
   fChain->SetBranchAddress("DeltaMed1h", &DeltaMed1h, &b_DeltaMed1h);
   fChain->SetBranchAddress("DeltaMed2h", &DeltaMed2h, &b_DeltaMed2h);
   fChain->SetBranchAddress("DeltaRMS", &DeltaRMS, &b_DeltaRMS);
   fChain->SetBranchAddress("DeltaRMS1h", &DeltaRMS1h, &b_DeltaRMS1h);
   fChain->SetBranchAddress("DeltaRMS2h", &DeltaRMS2h, &b_DeltaRMS2h);
   fChain->SetBranchAddress("CylFrac_1cm", &CylFrac_1cm, &b_CylFrac_1cm);
   fChain->SetBranchAddress("CylFrac1h_1cm", &CylFrac1h_1cm, &b_CylFrac1h_1cm);
   fChain->SetBranchAddress("CylFrac2h_1cm", &CylFrac2h_1cm, &b_CylFrac2h_1cm);
   fChain->SetBranchAddress("CylFrac_2cm", &CylFrac_2cm, &b_CylFrac_2cm);
   fChain->SetBranchAddress("CylFrac1h_2cm", &CylFrac1h_2cm, &b_CylFrac1h_2cm);
   fChain->SetBranchAddress("CylFrac2h_2cm", &CylFrac2h_2cm, &b_CylFrac2h_2cm);
   fChain->SetBranchAddress("CylFrac_3cm", &CylFrac_3cm, &b_CylFrac_3cm);
   fChain->SetBranchAddress("CylFrac1h_3cm", &CylFrac1h_3cm, &b_CylFrac1h_3cm);
   fChain->SetBranchAddress("CylFrac2h_3cm", &CylFrac2h_3cm, &b_CylFrac2h_3cm);
   fChain->SetBranchAddress("CylFrac_4cm", &CylFrac_4cm, &b_CylFrac_4cm);
   fChain->SetBranchAddress("CylFrac1h_4cm", &CylFrac1h_4cm, &b_CylFrac1h_4cm);
   fChain->SetBranchAddress("CylFrac2h_4cm", &CylFrac2h_4cm, &b_CylFrac2h_4cm);
   fChain->SetBranchAddress("CylFrac_5cm", &CylFrac_5cm, &b_CylFrac_5cm);
   fChain->SetBranchAddress("CylFrac1h_5cm", &CylFrac1h_5cm, &b_CylFrac1h_5cm);
   fChain->SetBranchAddress("CylFrac2h_5cm", &CylFrac2h_5cm, &b_CylFrac2h_5cm);
   fChain->SetBranchAddress("NeutrinoEnergy0", &NeutrinoEnergy0, &b_NeutrinoEnergy0);
   fChain->SetBranchAddress("NeutrinoEnergy1", &NeutrinoEnergy1, &b_NeutrinoEnergy1);
   fChain->SetBranchAddress("NeutrinoEnergy2", &NeutrinoEnergy2, &b_NeutrinoEnergy2);
   fChain->SetBranchAddress("SliceCaloEnergy0", &SliceCaloEnergy0, &b_SliceCaloEnergy0);
   fChain->SetBranchAddress("SliceCaloEnergy1", &SliceCaloEnergy1, &b_SliceCaloEnergy1);
   fChain->SetBranchAddress("SliceCaloEnergy2", &SliceCaloEnergy2, &b_SliceCaloEnergy2);
   fChain->SetBranchAddress("pi0_mcgamma0_e", &pi0_mcgamma0_e, &b_pi0_mcgamma0_e);
   fChain->SetBranchAddress("pi0_mcgamma0_px", &pi0_mcgamma0_px, &b_pi0_mcgamma0_px);
   fChain->SetBranchAddress("pi0_mcgamma0_py", &pi0_mcgamma0_py, &b_pi0_mcgamma0_py);
   fChain->SetBranchAddress("pi0_mcgamma0_pz", &pi0_mcgamma0_pz, &b_pi0_mcgamma0_pz);
   fChain->SetBranchAddress("pi0_mcrcdot0", &pi0_mcrcdot0, &b_pi0_mcrcdot0);
   fChain->SetBranchAddress("pi0_mcrce0", &pi0_mcrce0, &b_pi0_mcrce0);
   fChain->SetBranchAddress("pi0_mcgamma1_e", &pi0_mcgamma1_e, &b_pi0_mcgamma1_e);
   fChain->SetBranchAddress("pi0_mcgamma1_px", &pi0_mcgamma1_px, &b_pi0_mcgamma1_px);
   fChain->SetBranchAddress("pi0_mcgamma1_py", &pi0_mcgamma1_py, &b_pi0_mcgamma1_py);
   fChain->SetBranchAddress("pi0_mcgamma1_pz", &pi0_mcgamma1_pz, &b_pi0_mcgamma1_pz);
   fChain->SetBranchAddress("pi0_mcrcdot1", &pi0_mcrcdot1, &b_pi0_mcrcdot1);
   fChain->SetBranchAddress("pi0_mcrce1", &pi0_mcrce1, &b_pi0_mcrce1);
   fChain->SetBranchAddress("pi0_nshower", &pi0_nshower, &b_pi0_nshower);
   fChain->SetBranchAddress("pi0_ntrack", &pi0_ntrack, &b_pi0_ntrack);
   fChain->SetBranchAddress("pi0_ngamma", &pi0_ngamma, &b_pi0_ngamma);
   fChain->SetBranchAddress("pi0_radlen1", &pi0_radlen1, &b_pi0_radlen1);
   fChain->SetBranchAddress("pi0_radlen2", &pi0_radlen2, &b_pi0_radlen2);
   fChain->SetBranchAddress("pi0_dot1", &pi0_dot1, &b_pi0_dot1);
   fChain->SetBranchAddress("pi0_dot2", &pi0_dot2, &b_pi0_dot2);
   fChain->SetBranchAddress("pi0_energy1_Y", &pi0_energy1_Y, &b_pi0_energy1_Y);
   fChain->SetBranchAddress("pi0_energy2_Y", &pi0_energy2_Y, &b_pi0_energy2_Y);
   fChain->SetBranchAddress("pi0_dir1_x", &pi0_dir1_x, &b_pi0_dir1_x);
   fChain->SetBranchAddress("pi0_dir1_y", &pi0_dir1_y, &b_pi0_dir1_y);
   fChain->SetBranchAddress("pi0_dir1_z", &pi0_dir1_z, &b_pi0_dir1_z);
   fChain->SetBranchAddress("pi0_dir2_x", &pi0_dir2_x, &b_pi0_dir2_x);
   fChain->SetBranchAddress("pi0_dir2_y", &pi0_dir2_y, &b_pi0_dir2_y);
   fChain->SetBranchAddress("pi0_dir2_z", &pi0_dir2_z, &b_pi0_dir2_z);
   fChain->SetBranchAddress("pi0_dedx1_Y", &pi0_dedx1_Y, &b_pi0_dedx1_Y);
   fChain->SetBranchAddress("pi0_dedx2_Y", &pi0_dedx2_Y, &b_pi0_dedx2_Y);
   fChain->SetBranchAddress("pi0_dedx1_fit_Y", &pi0_dedx1_fit_Y, &b_pi0_dedx1_fit_Y);
   fChain->SetBranchAddress("pi0_dedx2_fit_Y", &pi0_dedx2_fit_Y, &b_pi0_dedx2_fit_Y);
   fChain->SetBranchAddress("pi0_energy1_V", &pi0_energy1_V, &b_pi0_energy1_V);
   fChain->SetBranchAddress("pi0_energy2_V", &pi0_energy2_V, &b_pi0_energy2_V);
   fChain->SetBranchAddress("pi0_dedx1_V", &pi0_dedx1_V, &b_pi0_dedx1_V);
   fChain->SetBranchAddress("pi0_dedx2_V", &pi0_dedx2_V, &b_pi0_dedx2_V);
   fChain->SetBranchAddress("pi0_dedx1_fit_V", &pi0_dedx1_fit_V, &b_pi0_dedx1_fit_V);
   fChain->SetBranchAddress("pi0_dedx2_fit_V", &pi0_dedx2_fit_V, &b_pi0_dedx2_fit_V);
   fChain->SetBranchAddress("pi0_energy1_U", &pi0_energy1_U, &b_pi0_energy1_U);
   fChain->SetBranchAddress("pi0_energy2_U", &pi0_energy2_U, &b_pi0_energy2_U);
   fChain->SetBranchAddress("pi0_dedx1_U", &pi0_dedx1_U, &b_pi0_dedx1_U);
   fChain->SetBranchAddress("pi0_dedx2_U", &pi0_dedx2_U, &b_pi0_dedx2_U);
   fChain->SetBranchAddress("pi0_dedx1_fit_U", &pi0_dedx1_fit_U, &b_pi0_dedx1_fit_U);
   fChain->SetBranchAddress("pi0_dedx2_fit_U", &pi0_dedx2_fit_U, &b_pi0_dedx2_fit_U);
   fChain->SetBranchAddress("pi0_shrscore1", &pi0_shrscore1, &b_pi0_shrscore1);
   fChain->SetBranchAddress("pi0_shrscore2", &pi0_shrscore2, &b_pi0_shrscore2);
   fChain->SetBranchAddress("pi0_gammadot", &pi0_gammadot, &b_pi0_gammadot);
   fChain->SetBranchAddress("pi0_mass_Y", &pi0_mass_Y, &b_pi0_mass_Y);
   fChain->SetBranchAddress("pi0_mass_V", &pi0_mass_V, &b_pi0_mass_V);
   fChain->SetBranchAddress("pi0_mass_U", &pi0_mass_U, &b_pi0_mass_U);
   fChain->SetBranchAddress("pi0_rc_vtx_x", &pi0_rc_vtx_x, &b_pi0_rc_vtx_x);
   fChain->SetBranchAddress("pi0_rc_vtx_y", &pi0_rc_vtx_y, &b_pi0_rc_vtx_y);
   fChain->SetBranchAddress("pi0_rc_vtx_z", &pi0_rc_vtx_z, &b_pi0_rc_vtx_z);
   fChain->SetBranchAddress("pi0truth_gamma_parent", &pi0truth_gamma_parent, &b_pi0truth_gamma_parent);
   fChain->SetBranchAddress("pi0truth_elec_edep", &pi0truth_elec_edep, &b_pi0truth_elec_edep);
   fChain->SetBranchAddress("pi0truth_elec_etot", &pi0truth_elec_etot, &b_pi0truth_elec_etot);
   fChain->SetBranchAddress("pi0truth_elec_dist", &pi0truth_elec_dist, &b_pi0truth_elec_dist);
   fChain->SetBranchAddress("pi0truth_elec_parent", &pi0truth_elec_parent, &b_pi0truth_elec_parent);
   fChain->SetBranchAddress("pi0truth_gamma1_tid", &pi0truth_gamma1_tid, &b_pi0truth_gamma1_tid);
   fChain->SetBranchAddress("pi0truth_gamma1_edep", &pi0truth_gamma1_edep, &b_pi0truth_gamma1_edep);
   fChain->SetBranchAddress("pi0truth_gamma1_etot", &pi0truth_gamma1_etot, &b_pi0truth_gamma1_etot);
   fChain->SetBranchAddress("pi0truth_gamma1_dist", &pi0truth_gamma1_dist, &b_pi0truth_gamma1_dist);
   fChain->SetBranchAddress("pi0truth_gamma1_elec1", &pi0truth_gamma1_elec1, &b_pi0truth_gamma1_elec1);
   fChain->SetBranchAddress("pi0truth_gamma1_elec2", &pi0truth_gamma1_elec2, &b_pi0truth_gamma1_elec2);
   fChain->SetBranchAddress("pi0truth_gamma1_xpos", &pi0truth_gamma1_xpos, &b_pi0truth_gamma1_xpos);
   fChain->SetBranchAddress("pi0truth_gamma1_ypos", &pi0truth_gamma1_ypos, &b_pi0truth_gamma1_ypos);
   fChain->SetBranchAddress("pi0truth_gamma1_zpos", &pi0truth_gamma1_zpos, &b_pi0truth_gamma1_zpos);
   fChain->SetBranchAddress("pi0truth_gamma2_tid", &pi0truth_gamma2_tid, &b_pi0truth_gamma2_tid);
   fChain->SetBranchAddress("pi0truth_gamma2_edep", &pi0truth_gamma2_edep, &b_pi0truth_gamma2_edep);
   fChain->SetBranchAddress("pi0truth_gamma2_etot", &pi0truth_gamma2_etot, &b_pi0truth_gamma2_etot);
   fChain->SetBranchAddress("pi0truth_gamma2_dist", &pi0truth_gamma2_dist, &b_pi0truth_gamma2_dist);
   fChain->SetBranchAddress("pi0truth_gamma2_elec1", &pi0truth_gamma2_elec1, &b_pi0truth_gamma2_elec1);
   fChain->SetBranchAddress("pi0truth_gamma2_elec2", &pi0truth_gamma2_elec2, &b_pi0truth_gamma2_elec2);
   fChain->SetBranchAddress("pi0truth_gamma2_xpos", &pi0truth_gamma2_xpos, &b_pi0truth_gamma2_xpos);
   fChain->SetBranchAddress("pi0truth_gamma2_ypos", &pi0truth_gamma2_ypos, &b_pi0truth_gamma2_ypos);
   fChain->SetBranchAddress("pi0truth_gamma2_zpos", &pi0truth_gamma2_zpos, &b_pi0truth_gamma2_zpos);
   fChain->SetBranchAddress("pi0truth_gammadot", &pi0truth_gammadot, &b_pi0truth_gammadot);
   fChain->SetBranchAddress("pi0truth_run", &pi0truth_run, &b_pi0truth_run);
   fChain->SetBranchAddress("pi0truth_sub", &pi0truth_sub, &b_pi0truth_sub);
   fChain->SetBranchAddress("pi0truth_evt", &pi0truth_evt, &b_pi0truth_evt);
   fChain->SetBranchAddress("nflag_pl1", &nflag_pl1, &b_nflag_pl1);
   fChain->SetBranchAddress("nnoise_pl1", &nnoise_pl1, &b_nnoise_pl1);
   fChain->SetBranchAddress("nslhits_pl1", &nslhits_pl1, &b_nslhits_pl1);
   fChain->SetBranchAddress("nslnoise_pl1", &nslnoise_pl1, &b_nslnoise_pl1);
   fChain->SetBranchAddress("nhits_pl1", &nhits_pl1, &b_nhits_pl1);
   fChain->SetBranchAddress("frac_slnoise_pl1", &frac_slnoise_pl1, &b_frac_slnoise_pl1);
   fChain->SetBranchAddress("secondshower_U_charge", &secondshower_U_charge, &b_secondshower_U_charge);
   fChain->SetBranchAddress("secondshower_U_nhit", &secondshower_U_nhit, &b_secondshower_U_nhit);
   fChain->SetBranchAddress("secondshower_U_vtxdist", &secondshower_U_vtxdist, &b_secondshower_U_vtxdist);
   fChain->SetBranchAddress("secondshower_U_eigenratio", &secondshower_U_eigenratio, &b_secondshower_U_eigenratio);
   fChain->SetBranchAddress("secondshower_U_dot", &secondshower_U_dot, &b_secondshower_U_dot);
   fChain->SetBranchAddress("secondshower_U_dir", &secondshower_U_dir, &b_secondshower_U_dir);
   fChain->SetBranchAddress("secondshower_V_charge", &secondshower_V_charge, &b_secondshower_V_charge);
   fChain->SetBranchAddress("secondshower_V_nhit", &secondshower_V_nhit, &b_secondshower_V_nhit);
   fChain->SetBranchAddress("secondshower_V_vtxdist", &secondshower_V_vtxdist, &b_secondshower_V_vtxdist);
   fChain->SetBranchAddress("secondshower_V_eigenratio", &secondshower_V_eigenratio, &b_secondshower_V_eigenratio);
   fChain->SetBranchAddress("secondshower_V_dot", &secondshower_V_dot, &b_secondshower_V_dot);
   fChain->SetBranchAddress("secondshower_V_dir", &secondshower_V_dir, &b_secondshower_V_dir);
   fChain->SetBranchAddress("secondshower_Y_charge", &secondshower_Y_charge, &b_secondshower_Y_charge);
   fChain->SetBranchAddress("secondshower_Y_nhit", &secondshower_Y_nhit, &b_secondshower_Y_nhit);
   fChain->SetBranchAddress("secondshower_Y_vtxdist", &secondshower_Y_vtxdist, &b_secondshower_Y_vtxdist);
   fChain->SetBranchAddress("secondshower_Y_eigenratio", &secondshower_Y_eigenratio, &b_secondshower_Y_eigenratio);
   fChain->SetBranchAddress("secondshower_Y_dot", &secondshower_Y_dot, &b_secondshower_Y_dot);
   fChain->SetBranchAddress("secondshower_Y_dir", &secondshower_Y_dir, &b_secondshower_Y_dir);
   fChain->SetBranchAddress("shr_dedx_u_v", &shr_dedx_u_v, &b_shr_dedx_u_v);
   fChain->SetBranchAddress("shr_dedx_v_v", &shr_dedx_v_v, &b_shr_dedx_v_v);
   fChain->SetBranchAddress("shr_dedx_y_v", &shr_dedx_y_v, &b_shr_dedx_y_v);
   fChain->SetBranchAddress("shr_energy_u_v", &shr_energy_u_v, &b_shr_energy_u_v);
   fChain->SetBranchAddress("shr_energy_v_v", &shr_energy_v_v, &b_shr_energy_v_v);
   fChain->SetBranchAddress("shr_energy_y_v", &shr_energy_y_v, &b_shr_energy_y_v);
   fChain->SetBranchAddress("shr_pfp_id_v", &shr_pfp_id_v, &b_shr_pfp_id_v);
   fChain->SetBranchAddress("shr_start_x_v", &shr_start_x_v, &b_shr_start_x_v);
   fChain->SetBranchAddress("shr_start_y_v", &shr_start_y_v, &b_shr_start_y_v);
   fChain->SetBranchAddress("shr_start_z_v", &shr_start_z_v, &b_shr_start_z_v);
   fChain->SetBranchAddress("shr_dist_v", &shr_dist_v, &b_shr_dist_v);
   fChain->SetBranchAddress("shr_start_U_v", &shr_start_U_v, &b_shr_start_U_v);
   fChain->SetBranchAddress("shr_start_V_v", &shr_start_V_v, &b_shr_start_V_v);
   fChain->SetBranchAddress("shr_px_v", &shr_px_v, &b_shr_px_v);
   fChain->SetBranchAddress("shr_py_v", &shr_py_v, &b_shr_py_v);
   fChain->SetBranchAddress("shr_pz_v", &shr_pz_v, &b_shr_pz_v);
   fChain->SetBranchAddress("shr_openangle_v", &shr_openangle_v, &b_shr_openangle_v);
   fChain->SetBranchAddress("shr_theta_v", &shr_theta_v, &b_shr_theta_v);
   fChain->SetBranchAddress("shr_phi_v", &shr_phi_v, &b_shr_phi_v);
   fChain->SetBranchAddress("shr_pitch_u_v", &shr_pitch_u_v, &b_shr_pitch_u_v);
   fChain->SetBranchAddress("shr_pitch_v_v", &shr_pitch_v_v, &b_shr_pitch_v_v);
   fChain->SetBranchAddress("shr_pitch_y_v", &shr_pitch_y_v, &b_shr_pitch_y_v);
   fChain->SetBranchAddress("shr_tkfit_nhits_v", &shr_tkfit_nhits_v, &b_shr_tkfit_nhits_v);
   fChain->SetBranchAddress("shr_tkfit_start_x_v", &shr_tkfit_start_x_v, &b_shr_tkfit_start_x_v);
   fChain->SetBranchAddress("shr_tkfit_start_y_v", &shr_tkfit_start_y_v, &b_shr_tkfit_start_y_v);
   fChain->SetBranchAddress("shr_tkfit_start_z_v", &shr_tkfit_start_z_v, &b_shr_tkfit_start_z_v);
   fChain->SetBranchAddress("shr_tkfit_start_U_v", &shr_tkfit_start_U_v, &b_shr_tkfit_start_U_v);
   fChain->SetBranchAddress("shr_tkfit_start_V_v", &shr_tkfit_start_V_v, &b_shr_tkfit_start_V_v);
   fChain->SetBranchAddress("shr_tkfit_theta_v", &shr_tkfit_theta_v, &b_shr_tkfit_theta_v);
   fChain->SetBranchAddress("shr_tkfit_phi_v", &shr_tkfit_phi_v, &b_shr_tkfit_phi_v);
   fChain->SetBranchAddress("shr_tkfit_pitch_u_v", &shr_tkfit_pitch_u_v, &b_shr_tkfit_pitch_u_v);
   fChain->SetBranchAddress("shr_tkfit_pitch_v_v", &shr_tkfit_pitch_v_v, &b_shr_tkfit_pitch_v_v);
   fChain->SetBranchAddress("shr_tkfit_pitch_y_v", &shr_tkfit_pitch_y_v, &b_shr_tkfit_pitch_y_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_u_v", &shr_tkfit_dedx_u_v, &b_shr_tkfit_dedx_u_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_v_v", &shr_tkfit_dedx_v_v, &b_shr_tkfit_dedx_v_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_y_v", &shr_tkfit_dedx_y_v, &b_shr_tkfit_dedx_y_v);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_u_v", &shr_tkfit_gap10_dedx_u_v, &b_shr_tkfit_gap10_dedx_u_v);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_v_v", &shr_tkfit_gap10_dedx_v_v, &b_shr_tkfit_gap10_dedx_v_v);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_y_v", &shr_tkfit_gap10_dedx_y_v, &b_shr_tkfit_gap10_dedx_y_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_nhits_u_v", &shr_tkfit_dedx_nhits_u_v, &b_shr_tkfit_dedx_nhits_u_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_nhits_v_v", &shr_tkfit_dedx_nhits_v_v, &b_shr_tkfit_dedx_nhits_v_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_nhits_y_v", &shr_tkfit_dedx_nhits_y_v, &b_shr_tkfit_dedx_nhits_y_v);
   fChain->SetBranchAddress("shr_llr_pid_u_v", &shr_llr_pid_u_v, &b_shr_llr_pid_u_v);
   fChain->SetBranchAddress("shr_llr_pid_v_v", &shr_llr_pid_v_v, &b_shr_llr_pid_v_v);
   fChain->SetBranchAddress("shr_llr_pid_y_v", &shr_llr_pid_y_v, &b_shr_llr_pid_y_v);
   fChain->SetBranchAddress("shr_llr_pid_v", &shr_llr_pid_v, &b_shr_llr_pid_v);
   fChain->SetBranchAddress("shr_llr_pid_score_v", &shr_llr_pid_score_v, &b_shr_llr_pid_score_v);
   fChain->SetBranchAddress("shr_moliere_avg_v", &shr_moliere_avg_v, &b_shr_moliere_avg_v);
   fChain->SetBranchAddress("shr_moliere_rms_v", &shr_moliere_rms_v, &b_shr_moliere_rms_v);
   fChain->SetBranchAddress("evnunhits", &evnunhits, &b_evnunhits);
   fChain->SetBranchAddress("evlepnhits", &evlepnhits, &b_evlepnhits);
   fChain->SetBranchAddress("evpronhits", &evpronhits, &b_evpronhits);
   fChain->SetBranchAddress("evpi1nhits", &evpi1nhits, &b_evpi1nhits);
   fChain->SetBranchAddress("evpi0nhits", &evpi0nhits, &b_evpi0nhits);
   fChain->SetBranchAddress("evneunhits", &evneunhits, &b_evneunhits);
   fChain->SetBranchAddress("evgamnhits", &evgamnhits, &b_evgamnhits);
   fChain->SetBranchAddress("evothnhits", &evothnhits, &b_evothnhits);
   fChain->SetBranchAddress("slnunhits", &slnunhits, &b_slnunhits);
   fChain->SetBranchAddress("sllepnhits", &sllepnhits, &b_sllepnhits);
   fChain->SetBranchAddress("slpronhits", &slpronhits, &b_slpronhits);
   fChain->SetBranchAddress("slpi1nhits", &slpi1nhits, &b_slpi1nhits);
   fChain->SetBranchAddress("slpi0nhits", &slpi0nhits, &b_slpi0nhits);
   fChain->SetBranchAddress("slneunhits", &slneunhits, &b_slneunhits);
   fChain->SetBranchAddress("slgamnhits", &slgamnhits, &b_slgamnhits);
   fChain->SetBranchAddress("slothnhits", &slothnhits, &b_slothnhits);
   fChain->SetBranchAddress("pfnunhits", &pfnunhits, &b_pfnunhits);
   fChain->SetBranchAddress("pflepnhits", &pflepnhits, &b_pflepnhits);
   fChain->SetBranchAddress("pfpronhits", &pfpronhits, &b_pfpronhits);
   fChain->SetBranchAddress("pfpi1nhits", &pfpi1nhits, &b_pfpi1nhits);
   fChain->SetBranchAddress("pfpi0nhits", &pfpi0nhits, &b_pfpi0nhits);
   fChain->SetBranchAddress("pfneunhits", &pfneunhits, &b_pfneunhits);
   fChain->SetBranchAddress("pfgamnhits", &pfgamnhits, &b_pfgamnhits);
   fChain->SetBranchAddress("pfothnhits", &pfothnhits, &b_pfothnhits);
   fChain->SetBranchAddress("nu_completeness_from_pfp", &nu_completeness_from_pfp, &b_nu_completeness_from_pfp);
   fChain->SetBranchAddress("nu_purity_from_pfp", &nu_purity_from_pfp, &b_nu_purity_from_pfp);
   fChain->SetBranchAddress("trk_bragg_p_v", &trk_bragg_p_v, &b_trk_bragg_p_v);
   fChain->SetBranchAddress("trk_bragg_mu_v", &trk_bragg_mu_v, &b_trk_bragg_mu_v);
   fChain->SetBranchAddress("trk_bragg_mip_v", &trk_bragg_mip_v, &b_trk_bragg_mip_v);
   fChain->SetBranchAddress("trk_pida_v", &trk_pida_v, &b_trk_pida_v);
   fChain->SetBranchAddress("trk_pid_chipr_v", &trk_pid_chipr_v, &b_trk_pid_chipr_v);
   fChain->SetBranchAddress("trk_pid_chipi_v", &trk_pid_chipi_v, &b_trk_pid_chipi_v);
   fChain->SetBranchAddress("trk_pid_chika_v", &trk_pid_chika_v, &b_trk_pid_chika_v);
   fChain->SetBranchAddress("trk_pid_chimu_v", &trk_pid_chimu_v, &b_trk_pid_chimu_v);
   fChain->SetBranchAddress("trk_bragg_p_u_v", &trk_bragg_p_u_v, &b_trk_bragg_p_u_v);
   fChain->SetBranchAddress("trk_bragg_mu_u_v", &trk_bragg_mu_u_v, &b_trk_bragg_mu_u_v);
   fChain->SetBranchAddress("trk_bragg_mip_u_v", &trk_bragg_mip_u_v, &b_trk_bragg_mip_u_v);
   fChain->SetBranchAddress("trk_pida_u_v", &trk_pida_u_v, &b_trk_pida_u_v);
   fChain->SetBranchAddress("trk_pid_chipr_u_v", &trk_pid_chipr_u_v, &b_trk_pid_chipr_u_v);
   fChain->SetBranchAddress("trk_pid_chipi_u_v", &trk_pid_chipi_u_v, &b_trk_pid_chipi_u_v);
   fChain->SetBranchAddress("trk_pid_chika_u_v", &trk_pid_chika_u_v, &b_trk_pid_chika_u_v);
   fChain->SetBranchAddress("trk_pid_chimu_u_v", &trk_pid_chimu_u_v, &b_trk_pid_chimu_u_v);
   fChain->SetBranchAddress("trk_bragg_p_v_v", &trk_bragg_p_v_v, &b_trk_bragg_p_v_v);
   fChain->SetBranchAddress("trk_bragg_mu_v_v", &trk_bragg_mu_v_v, &b_trk_bragg_mu_v_v);
   fChain->SetBranchAddress("trk_bragg_mip_v_v", &trk_bragg_mip_v_v, &b_trk_bragg_mip_v_v);
   fChain->SetBranchAddress("trk_pida_v_v", &trk_pida_v_v, &b_trk_pida_v_v);
   fChain->SetBranchAddress("trk_pid_chipr_v_v", &trk_pid_chipr_v_v, &b_trk_pid_chipr_v_v);
   fChain->SetBranchAddress("trk_pid_chipi_v_v", &trk_pid_chipi_v_v, &b_trk_pid_chipi_v_v);
   fChain->SetBranchAddress("trk_pid_chika_v_v", &trk_pid_chika_v_v, &b_trk_pid_chika_v_v);
   fChain->SetBranchAddress("trk_pid_chimu_v_v", &trk_pid_chimu_v_v, &b_trk_pid_chimu_v_v);
   fChain->SetBranchAddress("trk_pfp_id_v", &trk_pfp_id_v, &b_trk_pfp_id_v);
   fChain->SetBranchAddress("trk_dir_x_v", &trk_dir_x_v, &b_trk_dir_x_v);
   fChain->SetBranchAddress("trk_dir_y_v", &trk_dir_y_v, &b_trk_dir_y_v);
   fChain->SetBranchAddress("trk_dir_z_v", &trk_dir_z_v, &b_trk_dir_z_v);
   fChain->SetBranchAddress("trk_start_x_v", &trk_start_x_v, &b_trk_start_x_v);
   fChain->SetBranchAddress("trk_start_y_v", &trk_start_y_v, &b_trk_start_y_v);
   fChain->SetBranchAddress("trk_start_z_v", &trk_start_z_v, &b_trk_start_z_v);
   fChain->SetBranchAddress("trk_sce_start_x_v", &trk_sce_start_x_v, &b_trk_sce_start_x_v);
   fChain->SetBranchAddress("trk_sce_start_y_v", &trk_sce_start_y_v, &b_trk_sce_start_y_v);
   fChain->SetBranchAddress("trk_sce_start_z_v", &trk_sce_start_z_v, &b_trk_sce_start_z_v);
   fChain->SetBranchAddress("trk_end_x_v", &trk_end_x_v, &b_trk_end_x_v);
   fChain->SetBranchAddress("trk_end_y_v", &trk_end_y_v, &b_trk_end_y_v);
   fChain->SetBranchAddress("trk_end_z_v", &trk_end_z_v, &b_trk_end_z_v);
   fChain->SetBranchAddress("trk_sce_end_x_v", &trk_sce_end_x_v, &b_trk_sce_end_x_v);
   fChain->SetBranchAddress("trk_sce_end_y_v", &trk_sce_end_y_v, &b_trk_sce_end_y_v);
   fChain->SetBranchAddress("trk_sce_end_z_v", &trk_sce_end_z_v, &b_trk_sce_end_z_v);
   fChain->SetBranchAddress("trk_distance_v", &trk_distance_v, &b_trk_distance_v);
   fChain->SetBranchAddress("trk_theta_v", &trk_theta_v, &b_trk_theta_v);
   fChain->SetBranchAddress("trk_phi_v", &trk_phi_v, &b_trk_phi_v);
   fChain->SetBranchAddress("trk_len_v", &trk_len_v, &b_trk_len_v);
   fChain->SetBranchAddress("trk_mcs_muon_mom_v", &trk_mcs_muon_mom_v, &b_trk_mcs_muon_mom_v);
   fChain->SetBranchAddress("trk_range_muon_mom_v", &trk_range_muon_mom_v, &b_trk_range_muon_mom_v);
   fChain->SetBranchAddress("trk_energy_proton_v", &trk_energy_proton_v, &b_trk_energy_proton_v);
   fChain->SetBranchAddress("trk_energy_muon_v", &trk_energy_muon_v, &b_trk_energy_muon_v);
   fChain->SetBranchAddress("trk_calo_energy_u_v", &trk_calo_energy_u_v, &b_trk_calo_energy_u_v);
   fChain->SetBranchAddress("trk_calo_energy_v_v", &trk_calo_energy_v_v, &b_trk_calo_energy_v_v);
   fChain->SetBranchAddress("trk_calo_energy_y_v", &trk_calo_energy_y_v, &b_trk_calo_energy_y_v);
   fChain->SetBranchAddress("trk_llr_pid_u_v", &trk_llr_pid_u_v, &b_trk_llr_pid_u_v);
   fChain->SetBranchAddress("trk_llr_pid_v_v", &trk_llr_pid_v_v, &b_trk_llr_pid_v_v);
   fChain->SetBranchAddress("trk_llr_pid_y_v", &trk_llr_pid_y_v, &b_trk_llr_pid_y_v);
   fChain->SetBranchAddress("trk_llr_pid_v", &trk_llr_pid_v, &b_trk_llr_pid_v);
   fChain->SetBranchAddress("trk_llr_pid_score_v", &trk_llr_pid_score_v, &b_trk_llr_pid_score_v);
   fChain->SetBranchAddress("bdt_nuNCpi0", &bdt_nuNCpi0, &b_bdt_nuNCpi0);
   fChain->SetBranchAddress("bdt_numuCCpi0", &bdt_numuCCpi0, &b_bdt_numuCCpi0);
   fChain->SetBranchAddress("bdt_numuCC", &bdt_numuCC, &b_bdt_numuCC);
   fChain->SetBranchAddress("bdt_ext", &bdt_ext, &b_bdt_ext);
   fChain->SetBranchAddress("bdt_cosmic", &bdt_cosmic, &b_bdt_cosmic);
   fChain->SetBranchAddress("bdt_global", &bdt_global, &b_bdt_global);
   fChain->SetBranchAddress("pass_antibdt_filter", &pass_antibdt_filter, &b_bdt_global);
   fChain->SetBranchAddress("bdt_pi0_np", &bdt_pi0_np, &b_bdt_pi0_np);
   fChain->SetBranchAddress("bdt_nonpi0_np", &bdt_nonpi0_np, &b_bdt_nonpi0_np);
   fChain->SetBranchAddress("bdt_bkg_0p", &bdt_bkg_0p, &b_bdt_bkg_0p);
   fChain->SetBranchAddress("anglediff_Y", &anglediff_Y, &b_anglediff_Y);
   fChain->SetBranchAddress("anglediff_V", &anglediff_V, &b_anglediff_V);
   fChain->SetBranchAddress("anglediff_U", &anglediff_U, &b_anglediff_U);
   fChain->SetBranchAddress("trkpid", &trkpid, &b_trkpid);
   Notify();
}

Bool_t twoproton_nuwro_overlay::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void twoproton_nuwro_overlay::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t twoproton_nuwro_overlay::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef twoproton_nuwro_overlay_cxx
