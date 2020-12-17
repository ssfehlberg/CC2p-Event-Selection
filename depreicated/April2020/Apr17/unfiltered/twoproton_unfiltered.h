//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr  9 11:36:39 2020 by ROOT version 6.12/06
// from TTree tree/
// found on file: /pnfs/uboone/scratch/users/sfehlber/2020/April/Apr08/CC2p/overlay_analyzer/final.root
//////////////////////////////////////////////////////////

#ifndef twoproton_unfiltered_h
#define twoproton_unfiltered_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class twoproton_unfiltered {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           subrun;
   Int_t           event;
   Int_t           n_pfp_per_event;
   Int_t           n_trk_per_event;
   Int_t           n_shower_per_event;
   Int_t           n_neutrinos_per_event;
   Int_t           n_nu_per_event;
   Int_t           n_nu_pfp_per_event;
   Int_t           nu_PDG_per_event;
   Int_t           mc_ccnc;
   Int_t           mc_mode;
   Int_t           mc_interactiontype;
   Int_t           mc_hitnuc;
   Float_t         mc_q2;
   Float_t         mc_X;
   Float_t         mc_Y;
   Float_t         mc_Pt;
   Float_t         mc_nu_vtxx;
   Float_t         mc_nu_vtxy;
   Float_t         mc_nu_vtxz;
   Float_t         mc_nu_vtxx_sce;
   Float_t         mc_nu_vtxy_sce;
   Float_t         mc_nu_vtxz_sce;
   Float_t         mc_enu;
   Float_t         mc_wgt;
   Float_t         mc_wgt_cv;
   Float_t         mc_wgt_0_EtaNCEL;
   Float_t         mc_wgt_0_FrAbs_N;
   Float_t         mc_wgt_0_FrAbs_pi;
   Float_t         mc_wgt_0_FrCEx_N;
   Float_t         mc_wgt_0_FrCEx_pi;
   Float_t         mc_wgt_0_FrInel_N;
   Float_t         mc_wgt_0_FrInel_pi;
   Float_t         mc_wgt_0_FrPiProd_N;
   Float_t         mc_wgt_0_FrPiProd_pi;
   Float_t         mc_wgt_0_MFP_N;
   Float_t         mc_wgt_0_MFP_pi;
   Float_t         mc_wgt_0_MaCCQE;
   Float_t         mc_wgt_0_MaCCRES;
   Float_t         mc_wgt_0_MaNCEL;
   Float_t         mc_wgt_0_MaNCRES;
   Float_t         mc_wgt_0_MvCCRES;
   Float_t         mc_wgt_0_MvNCRES;
   Float_t         mc_wgt_0_NonRESBGvnCC1pi;
   Float_t         mc_wgt_0_NonRESBGvnCC2pi;
   Float_t         mc_wgt_0_NonRESBGvnNC1pi;
   Float_t         mc_wgt_0_NonRESBGvnNC2pi;
   Float_t         mc_wgt_0_NonRESBGvpCC1pi;
   Float_t         mc_wgt_0_NonRESBGvpCC2pi;
   Float_t         mc_wgt_0_NonRESBGvpNC1pi;
   Float_t         mc_wgt_0_NonRESBGvpNC2pi;
   Float_t         mc_wgt_0_NormCCMEC;
   Float_t         mc_wgt_0_NormNCMEC;
   Float_t         mc_wgt_1_EtaNCEL;
   Float_t         mc_wgt_1_FrAbs_N;
   Float_t         mc_wgt_1_FrAbs_pi;
   Float_t         mc_wgt_1_FrCEx_N;
   Float_t         mc_wgt_1_FrCEx_pi;
   Float_t         mc_wgt_1_FrInel_N;
   Float_t         mc_wgt_1_FrInel_pi;
   Float_t         mc_wgt_1_FrPiProd_N;
   Float_t         mc_wgt_1_FrPiProd_pi;
   Float_t         mc_wgt_1_MFP_N;
   Float_t         mc_wgt_1_MFP_pi;
   Float_t         mc_wgt_1_MaCCQE;
   Float_t         mc_wgt_1_MaCCRES;
   Float_t         mc_wgt_1_MaNCEL;
   Float_t         mc_wgt_1_MaNCRES;
   Float_t         mc_wgt_1_MvCCRES;
   Float_t         mc_wgt_1_MvNCRES;
   Float_t         mc_wgt_1_NonRESBGvnCC1pi;
   Float_t         mc_wgt_1_NonRESBGvnCC2pi;
   Float_t         mc_wgt_1_NonRESBGvnNC1pi;
   Float_t         mc_wgt_1_NonRESBGvnNC2pi;
   Float_t         mc_wgt_1_NonRESBGvpCC1pi;
   Float_t         mc_wgt_1_NonRESBGvpCC2pi;
   Float_t         mc_wgt_1_NonRESBGvpNC1pi;
   Float_t         mc_wgt_1_NonRESBGvpNC2pi;
   Float_t         mc_wgt_1_NormCCMEC;
   Float_t         mc_wgt_1_NormNCMEC;
   Int_t           evtwgt_genie_pm1_nfunc;
   vector<string>  *evtwgt_genie_pm1_funcname;
   vector<int>     *evtwgt_genie_pm1_nweight;
   vector<vector<double> > *evtwgt_genie_pm1_weight;
   Int_t           evtwgt_genie_multisim_nfunc;
   vector<string>  *evtwgt_genie_multisim_funcname;
   vector<int>     *evtwgt_genie_multisim_nweight;
   vector<vector<double> > *evtwgt_genie_multisim_weight;
   Int_t           evtwgt_flux_multisim_nfunc;
   vector<string>  *evtwgt_flux_multisim_funcname;
   vector<int>     *evtwgt_flux_multisim_nweight;
   vector<vector<double> > *evtwgt_flux_multisim_weight;
   Int_t           mc_nupdg;
   Int_t           mc_n_muon;
   Int_t           mc_n_proton;
   Int_t           mc_n_pionpm;
   Int_t           mc_n_pion0;
   Int_t           mc_n_electron;
   Int_t           mc_n_neutron;
   Int_t           mc_n_threshold_muon;
   Int_t           mc_n_threshold_proton;
   Int_t           mc_n_threshold_pionpm;
   Int_t           mc_n_threshold_pion0;
   Int_t           mc_n_threshold_electron;
   Int_t           mc_n_threshold_neutron;
   vector<float>   *mc_g4_mom_all;
   vector<float>   *mc_g4_mom_muon;
   vector<float>   *mc_g4_mom_proton;
   vector<float>   *mc_g4_mom_electron;
   vector<float>   *mc_g4_mom_pionpm;
   vector<float>   *mc_g4_mom_pion0;
   vector<float>   *mc_g4_mom_neutron;
   vector<float>   *mc_g4_E;
   vector<float>   *mc_g4_p;
   vector<float>   *mc_g4_mass;
   vector<float>   *mc_g4_phi;
   vector<float>   *mc_g4_theta;
   vector<int>     *mc_g4_pdg;
   vector<float>   *mc_g4_start_x;
   vector<float>   *mc_g4_start_y;
   vector<float>   *mc_g4_start_z;
   vector<float>   *mc_g4_end_x;
   vector<float>   *mc_g4_end_y;
   vector<float>   *mc_g4_end_z;
   vector<float>   *mc_g4_start_x_sce;
   vector<float>   *mc_g4_start_y_sce;
   vector<float>   *mc_g4_start_z_sce;
   vector<float>   *mc_g4_end_x_sce;
   vector<float>   *mc_g4_end_y_sce;
   vector<float>   *mc_g4_end_z_sce;
   vector<bool>    *is_from_nu_slice;
   vector<bool>    *is_primary;
   vector<bool>    *is_contained;
   vector<int>     *mc_pdg;
   vector<int>     *mc_primary;
   vector<int>     *mc_origin;
   vector<float>   *mc_length;
   vector<float>   *mc_start_x;
   vector<float>   *mc_start_y;
   vector<float>   *mc_start_z;
   vector<float>   *mc_end_x;
   vector<float>   *mc_end_y;
   vector<float>   *mc_end_z;
   vector<float>   *mc_start_x_sce;
   vector<float>   *mc_start_y_sce;
   vector<float>   *mc_start_z_sce;
   vector<float>   *mc_end_x_sce;
   vector<float>   *mc_end_y_sce;
   vector<float>   *mc_end_z_sce;
   vector<float>   *mc_theta;
   vector<float>   *mc_phi;
   vector<float>   *mc_ke;
   vector<float>   *mc_mom;
   vector<int>     *n_pfp;
   vector<int>     *n_trk;
   vector<int>     *id_pfp;
   vector<int>     *n_shower;
   vector<int>     *parentPDG;
   vector<float>   *trk_score;
   vector<float>   *KE_len;
   vector<float>   *dislen_ratio;
   vector<float>   *reco_q2;
   vector<float>   *top_score;
   vector<int>     *n_daughters;
   Int_t           vtx_n_pfp;
   vector<bool>    *has_shower;
   Float_t         reco_nu_vtxx;
   Float_t         reco_nu_vtxy;
   Float_t         reco_nu_vtxz;
   vector<float>   *reco_length;
   vector<float>   *reco_start_x;
   vector<float>   *reco_start_y;
   vector<float>   *reco_start_z;
   vector<float>   *reco_end_x;
   vector<float>   *reco_end_y;
   vector<float>   *reco_end_z;
   vector<float>   *reco_theta;
   vector<float>   *reco_phi;
   vector<float>   *reco_ke;
   vector<float>   *reco_mom;
   vector<float>   *reco_mom_muon;
   vector<float>   *reco_mom_proton;
   vector<float>   *reco_mom_pion;
   vector<int>     *nhits_0;
   vector<int>     *nhits_1;
   vector<int>     *nhits_2;
   vector<float>   *chi2_p_0;
   vector<float>   *chi2_p_1;
   vector<float>   *chi2_p_2;
   vector<float>   *chi2_mu_0;
   vector<float>   *chi2_mu_1;
   vector<float>   *chi2_mu_2;
   vector<float>   *LL3;
   vector<float>   *LL_p_0;
   vector<float>   *LL_p_1;
   vector<float>   *LL_p_2;
   vector<float>   *LL_mip_0;
   vector<float>   *LL_mip_1;
   vector<float>   *LL_mip_2;
   vector<float>   *LL_mu_0;
   vector<float>   *LL_mu_1;
   vector<float>   *LL_mu_2;
   vector<float>   *LL_back_p_0;
   vector<float>   *LL_back_p_1;
   vector<float>   *LL_back_p_2;
   vector<float>   *LL_back_mu_0;
   vector<float>   *LL_back_mu_1;
   vector<float>   *LL_back_mu_2;
   vector<float>   *TM_dedx_0;
   vector<float>   *TM_dedx_1;
   vector<float>   *TM_dedx_2;
   vector<float>   *PIDA_0;
   vector<float>   *PIDA_1;
   vector<float>   *PIDA_2;
   vector<float>   *start_dedx_0;
   vector<float>   *start_dedx_1;
   vector<float>   *start_dedx_2;
   vector<float>   *end_dedx_0;
   vector<float>   *end_dedx_1;
   vector<float>   *end_dedx_2;
   vector<float>   *ratio_dedx_0;
   vector<float>   *ratio_dedx_1;
   vector<float>   *ratio_dedx_2;
   vector<float>   *avg_dedx_0;
   vector<float>   *avg_dedx_1;
   vector<float>   *avg_dedx_2;
   vector<float>   *total_dedx_0;
   vector<float>   *total_dedx_1;
   vector<float>   *total_dedx_2;
   vector<float>   *KE_calo_0;
   vector<float>   *KE_calo_1;
   vector<float>   *KE_calo_2;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_event;   //!
   TBranch        *b_n_pfp_per_event;   //!
   TBranch        *b_n_trk_per_event;   //!
   TBranch        *b_n_shower_per_event;   //!
   TBranch        *b_n_neutrinos_per_event;   //!
   TBranch        *b_n_nu_per_event;   //!
   TBranch        *b_n_nu_pfp_per_event;   //!
   TBranch        *b_nu_PDG_per_event;   //!
   TBranch        *b_mc_ccnc;   //!
   TBranch        *b_mc_mode;   //!
   TBranch        *b_mc_interactiontype;   //!
   TBranch        *b_mc_hitnuc;   //!
   TBranch        *b_mc_q2;   //!
   TBranch        *b_mc_X;   //!
   TBranch        *b_mc_Y;   //!
   TBranch        *b_mc_Pt;   //!
   TBranch        *b_mc_nu_vtxx;   //!
   TBranch        *b_mc_nu_vtxy;   //!
   TBranch        *b_mc_nu_vtxz;   //!
   TBranch        *b_mc_nu_vtxx_sce;   //!
   TBranch        *b_mc_nu_vtxy_sce;   //!
   TBranch        *b_mc_nu_vtxz_sce;   //!
   TBranch        *b_mc_enu;   //!
   TBranch        *b_mc_wgt;   //!
   TBranch        *b_mc_wgt_cv;   //!
   TBranch        *b_mc_wgt_0_EtaNCEL;   //!
   TBranch        *b_mc_wgt_0_FrAbs_N;   //!
   TBranch        *b_mc_wgt_0_FrAbs_pi;   //!
   TBranch        *b_mc_wgt_0_FrCEx_N;   //!
   TBranch        *b_mc_wgt_0_FrCEx_pi;   //!
   TBranch        *b_mc_wgt_0_FrInel_N;   //!
   TBranch        *b_mc_wgt_0_FrInel_pi;   //!
   TBranch        *b_mc_wgt_0_FrPiProd_N;   //!
   TBranch        *b_mc_wgt_0_FrPiProd_pi;   //!
   TBranch        *b_mc_wgt_0_MFP_N;   //!
   TBranch        *b_mc_wgt_0_MFP_pi;   //!
   TBranch        *b_mc_wgt_0_MaCCQE;   //!
   TBranch        *b_mc_wgt_0_MaCCRES;   //!
   TBranch        *b_mc_wgt_0_MaNCEL;   //!
   TBranch        *b_mc_wgt_0_MaNCRES;   //!
   TBranch        *b_mc_wgt_0_MvCCRES;   //!
   TBranch        *b_mc_wgt_0_MvNCRES;   //!
   TBranch        *b_mc_wgt_0_NonRESBGvnCC1pi;   //!
   TBranch        *b_mc_wgt_0_NonRESBGvnCC2pi;   //!
   TBranch        *b_mc_wgt_0_NonRESBGvnNC1pi;   //!
   TBranch        *b_mc_wgt_0_NonRESBGvnNC2pi;   //!
   TBranch        *b_mc_wgt_0_NonRESBGvpCC1pi;   //!
   TBranch        *b_mc_wgt_0_NonRESBGvpCC2pi;   //!
   TBranch        *b_mc_wgt_0_NonRESBGvpNC1pi;   //!
   TBranch        *b_mc_wgt_0_NonRESBGvpNC2pi;   //!
   TBranch        *b_mc_wgt_0_NormCCMEC;   //!
   TBranch        *b_mc_wgt_0_NormNCMEC;   //!
   TBranch        *b_mc_wgt_1_EtaNCEL;   //!
   TBranch        *b_mc_wgt_1_FrAbs_N;   //!
   TBranch        *b_mc_wgt_1_FrAbs_pi;   //!
   TBranch        *b_mc_wgt_1_FrCEx_N;   //!
   TBranch        *b_mc_wgt_1_FrCEx_pi;   //!
   TBranch        *b_mc_wgt_1_FrInel_N;   //!
   TBranch        *b_mc_wgt_1_FrInel_pi;   //!
   TBranch        *b_mc_wgt_1_FrPiProd_N;   //!
   TBranch        *b_mc_wgt_1_FrPiProd_pi;   //!
   TBranch        *b_mc_wgt_1_MFP_N;   //!
   TBranch        *b_mc_wgt_1_MFP_pi;   //!
   TBranch        *b_mc_wgt_1_MaCCQE;   //!
   TBranch        *b_mc_wgt_1_MaCCRES;   //!
   TBranch        *b_mc_wgt_1_MaNCEL;   //!
   TBranch        *b_mc_wgt_1_MaNCRES;   //!
   TBranch        *b_mc_wgt_1_MvCCRES;   //!
   TBranch        *b_mc_wgt_1_MvNCRES;   //!
   TBranch        *b_mc_wgt_1_NonRESBGvnCC1pi;   //!
   TBranch        *b_mc_wgt_1_NonRESBGvnCC2pi;   //!
   TBranch        *b_mc_wgt_1_NonRESBGvnNC1pi;   //!
   TBranch        *b_mc_wgt_1_NonRESBGvnNC2pi;   //!
   TBranch        *b_mc_wgt_1_NonRESBGvpCC1pi;   //!
   TBranch        *b_mc_wgt_1_NonRESBGvpCC2pi;   //!
   TBranch        *b_mc_wgt_1_NonRESBGvpNC1pi;   //!
   TBranch        *b_mc_wgt_1_NonRESBGvpNC2pi;   //!
   TBranch        *b_mc_wgt_1_NormCCMEC;   //!
   TBranch        *b_mc_wgt_1_NormNCMEC;   //!
   TBranch        *b_evtwgt_genie_pm1_nfunc;   //!
   TBranch        *b_evtwgt_genie_pm1_funcname;   //!
   TBranch        *b_evtwgt_genie_pm1_nweight;   //!
   TBranch        *b_evtwgt_genie_pm1_weight;   //!
   TBranch        *b_evtwgt_genie_multisim_nfunc;   //!
   TBranch        *b_evtwgt_genie_multisim_funcname;   //!
   TBranch        *b_evtwgt_genie_multisim_nweight;   //!
   TBranch        *b_evtwgt_genie_multisim_weight;   //!
   TBranch        *b_evtwgt_flux_multisim_nfunc;   //!
   TBranch        *b_evtwgt_flux_multisim_funcname;   //!
   TBranch        *b_evtwgt_flux_multisim_nweight;   //!
   TBranch        *b_evtwgt_flux_multisim_weight;   //!
   TBranch        *b_mc_nupdg;   //!
   TBranch        *b_mc_n_muon;   //!
   TBranch        *b_mc_n_proton;   //!
   TBranch        *b_mc_n_pionpm;   //!
   TBranch        *b_mc_n_pion0;   //!
   TBranch        *b_mc_n_electron;   //!
   TBranch        *b_mc_n_neutron;   //!
   TBranch        *b_mc_n_threshold_muon;   //!
   TBranch        *b_mc_n_threshold_proton;   //!
   TBranch        *b_mc_n_threshold_pionpm;   //!
   TBranch        *b_mc_n_threshold_pion0;   //!
   TBranch        *b_mc_n_threshold_electron;   //!
   TBranch        *b_mc_n_threshold_neutron;   //!
   TBranch        *b_mc_g4_mom_all;   //!
   TBranch        *b_mc_g4_mom_muon;   //!
   TBranch        *b_mc_g4_mom_proton;   //!
   TBranch        *b_mc_g4_mom_electron;   //!
   TBranch        *b_mc_g4_mom_pionpm;   //!
   TBranch        *b_mc_g4_mom_pion0;   //!
   TBranch        *b_mc_g4_mom_neutron;   //!
   TBranch        *b_mc_g4_E;   //!
   TBranch        *b_mc_g4_p;   //!
   TBranch        *b_mc_g4_mass;   //!
   TBranch        *b_mc_g4_phi;   //!
   TBranch        *b_mc_g4_theta;   //!
   TBranch        *b_mc_g4_pdg;   //!
   TBranch        *b_mc_g4_start_x;   //!
   TBranch        *b_mc_g4_start_y;   //!
   TBranch        *b_mc_g4_start_z;   //!
   TBranch        *b_mc_g4_end_x;   //!
   TBranch        *b_mc_g4_end_y;   //!
   TBranch        *b_mc_g4_end_z;   //!
   TBranch        *b_mc_g4_start_x_sce;   //!
   TBranch        *b_mc_g4_start_y_sce;   //!
   TBranch        *b_mc_g4_start_z_sce;   //!
   TBranch        *b_mc_g4_end_x_sce;   //!
   TBranch        *b_mc_g4_end_y_sce;   //!
   TBranch        *b_mc_g4_end_z_sce;   //!
   TBranch        *b_is_from_nu_slice;   //!
   TBranch        *b_is_primary;   //!
   TBranch        *b_is_contained;   //!
   TBranch        *b_mc_pdg;   //!
   TBranch        *b_mc_primary;   //!
   TBranch        *b_mc_origin;   //!
   TBranch        *b_mc_length;   //!
   TBranch        *b_mc_start_x;   //!
   TBranch        *b_mc_start_y;   //!
   TBranch        *b_mc_start_z;   //!
   TBranch        *b_mc_end_x;   //!
   TBranch        *b_mc_end_y;   //!
   TBranch        *b_mc_end_z;   //!
   TBranch        *b_mc_start_x_sce;   //!
   TBranch        *b_mc_start_y_sce;   //!
   TBranch        *b_mc_start_z_sce;   //!
   TBranch        *b_mc_end_x_sce;   //!
   TBranch        *b_mc_end_y_sce;   //!
   TBranch        *b_mc_end_z_sce;   //!
   TBranch        *b_mc_theta;   //!
   TBranch        *b_mc_phi;   //!
   TBranch        *b_mc_ke;   //!
   TBranch        *b_mc_mom;   //!
   TBranch        *b_n_pfp;   //!
   TBranch        *b_n_trk;   //!
   TBranch        *b_id_pfp;   //!
   TBranch        *b_n_shower;   //!
   TBranch        *b_parentPDG;   //!
   TBranch        *b_trk_score;   //!
   TBranch        *b_KE_len;   //!
   TBranch        *b_dislen_ratio;   //!
   TBranch        *b_reco_q2;   //!
   TBranch        *b_top_score;   //!
   TBranch        *b_n_daughters;   //!
   TBranch        *b_vtx_n_pfp;   //!
   TBranch        *b_has_shower;   //!
   TBranch        *b_reco_nu_vtxx;   //!
   TBranch        *b_reco_nu_vtxy;   //!
   TBranch        *b_reco_nu_vtxz;   //!
   TBranch        *b_reco_length;   //!
   TBranch        *b_reco_start_x;   //!
   TBranch        *b_reco_start_y;   //!
   TBranch        *b_reco_start_z;   //!
   TBranch        *b_reco_end_x;   //!
   TBranch        *b_reco_end_y;   //!
   TBranch        *b_reco_end_z;   //!
   TBranch        *b_reco_theta;   //!
   TBranch        *b_reco_phi;   //!
   TBranch        *b_reco_ke;   //!
   TBranch        *b_reco_mom;   //!
   TBranch        *b_reco_mom_muon;   //!
   TBranch        *b_reco_mom_proton;   //!
   TBranch        *b_reco_mom_pion;   //!
   TBranch        *b_nhits_0;   //!
   TBranch        *b_nhits_1;   //!
   TBranch        *b_nhits_2;   //!
   TBranch        *b_chi2_p_0;   //!
   TBranch        *b_chi2_p_1;   //!
   TBranch        *b_chi2_p_2;   //!
   TBranch        *b_chi2_mu_0;   //!
   TBranch        *b_chi2_mu_1;   //!
   TBranch        *b_chi2_mu_2;   //!
   TBranch        *b_LL3;   //!
   TBranch        *b_LL_p_0;   //!
   TBranch        *b_LL_p_1;   //!
   TBranch        *b_LL_p_2;   //!
   TBranch        *b_LL_mip_0;   //!
   TBranch        *b_LL_mip_1;   //!
   TBranch        *b_LL_mip_2;   //!
   TBranch        *b_LL_mu_0;   //!
   TBranch        *b_LL_mu_1;   //!
   TBranch        *b_LL_mu_2;   //!
   TBranch        *b_LL_back_p_0;   //!
   TBranch        *b_LL_back_p_1;   //!
   TBranch        *b_LL_back_p_2;   //!
   TBranch        *b_LL_back_mu_0;   //!
   TBranch        *b_LL_back_mu_1;   //!
   TBranch        *b_LL_back_mu_2;   //!
   TBranch        *b_TM_dedx_0;   //!
   TBranch        *b_TM_dedx_1;   //!
   TBranch        *b_TM_dedx_2;   //!
   TBranch        *b_PIDA_0;   //!
   TBranch        *b_PIDA_1;   //!
   TBranch        *b_PIDA_2;   //!
   TBranch        *b_start_dedx_0;   //!
   TBranch        *b_start_dedx_1;   //!
   TBranch        *b_start_dedx_2;   //!
   TBranch        *b_end_dedx_0;   //!
   TBranch        *b_end_dedx_1;   //!
   TBranch        *b_end_dedx_2;   //!
   TBranch        *b_ratio_dedx_0;   //!
   TBranch        *b_ratio_dedx_1;   //!
   TBranch        *b_ratio_dedx_2;   //!
   TBranch        *b_avg_dedx_0;   //!
   TBranch        *b_avg_dedx_1;   //!
   TBranch        *b_avg_dedx_2;   //!
   TBranch        *b_total_dedx_0;   //!
   TBranch        *b_total_dedx_1;   //!
   TBranch        *b_total_dedx_2;   //!
   TBranch        *b_KE_calo_0;   //!
   TBranch        *b_KE_calo_1;   //!
   TBranch        *b_KE_calo_2;   //!

   twoproton_unfiltered(TTree *tree=0);
   virtual ~twoproton_unfiltered();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     In_FV(float FV_edge, float x, float y, float z);
   virtual void     Define_Histograms();
   virtual void     Fill_Histograms_Mine(int i);
   virtual void     Fill_Histograms_Raquel(int i);
   virtual void     Fill_Mine(int i, int j);
   virtual void     Fill_Raquel(int i, int j);
   virtual void     Write_Histograms();

   //this is where you can define all your constants that you want
 private:

   //Boolean for inFV
   bool fv;

   //POT Weight:
   static const int wgt = 1;

   //Total Histograms                                                                                                                                                                                     
   static const int num = 4;
   const char * total[num] = {"npfp","vtx_npfp","ntrack","nshower"};
   TH1D * h_overlay[num];

   //Correlation Histograms                                                                                                                                                                                 
   static const int num2d = 3;
   const char * total2d[num2d] = {"reco","truth","truth_sce"};
   const char * labelx[num2d] = {"reco x","true x","true x + sce"};
   const char * labely[num2d] = {"reco y","true y","true y + sce"};
   TH2D *h_correlation_overlay[num2d];

   //Now to define all the specific channels and their histograms                                                                                        
   //Since I am basically a plot factory now, I am going to try and do this the smart way                                                                
   ///////////////////////////////////////////////////////////////////////////////////////                                                               
   static const int  number=2; //before and after selection                                                                                                      
   static const int  number2 = 11; //categories I defined                                                                                                        
   static const int  number3 = 10; //categories raquel defined
   const char * point[number] ={"_before_selection","_after_selection"}; //this defines histograms before and after the selection                        
   const char * channel[number2]={"_total","_cc0p0pi","_cc1p0pi","_cc2p0pi","_ccNp0pi",
                                "_ccNp1pi","_ccNpNpi","_ccnue","_outfv","_nc","_other"}; //these are the channels I defined                              

   const char * channel2[number3] = {"_total","_ccQE","_ccCOH","_ccMEC","_ccRES","_ccDIS",
				     "_ccNue","_nc","_outfv","_other"};//these are the channels that raquel defined

   TH1D* h_vtx_x[number][number2]; //reco x
   TH1D* h_vtx_y[number][number2]; //reco y 
   TH1D* h_vtx_z[number][number2]; //reco z
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

   TH1D* h_vtx_x_raquel[number][number3]; //reco x                                                                                                                                                       
   TH1D* h_vtx_y_raquel[number][number3]; //reco y                                                                                                                                                      
   TH1D* h_vtx_z_raquel[number][number3]; //reco z                                                                                                                                                       
   TH1D* h_vtx_x_mc_raquel[number][number3]; //mc x                                                                                                                                                      
   TH1D* h_vtx_y_mc_raquel[number][number3]; //mc y                                                                                                                                                      
   TH1D* h_vtx_z_mc_raquel[number][number3]; //mc z                                                                                                                                                      
   TH1D* h_vtx_x_mc_sce_raquel[number][number3]; //mc+sce x                                                                                                                                              
   TH1D* h_vtx_y_mc_sce_raquel[number][number3]; //mc+sce y                                                                                                                                              
   TH1D* h_vtx_z_mc_sce_raquel[number][number3]; //mc+sce z                                                                                                                                              
   TH1D* h_q2_raquel[number][number3]; //mc q2                                                                                                                                                           
   TH1D* h_X_raquel[number][number3]; //mc x                                                                                                                                                             
   TH1D* h_Y_raquel[number][number3]; //mc y                                                                                                                                                            
   TH1D* h_Pt_raquel[number][number3]; //mc Pt   

   //number of generated event/channel                                                                                                                                                               
   int cc0p0pi[number] = {0};
   int cc1p0pi[number] = {0};
   int cc2p0pi[number] = {0};
   int ccNp0pi[number] = {0};
   int ccNp1pi[number] = {0};
   int ccNpNpi[number] = {0};
   int ccnue[number] = {0};
   int outfv[number] = {0};
   int nc[number] = {0};
   int other[number] = {0};
   
   //number of generated event/channel                                                                                                                                                              
   int qel[number] = {0};
   int res[number] = {0};
   int mec[number] = {0};
   int coh[number] = {0};
   int dis[number] = {0};
   int ccnue_raquel[number] = {0};
   int outfv_raquel[number] = {0};
   int nc_raquel[number] = {0};
   int other_raquel[number] = {0};

};

#endif

#ifdef twoproton_unfiltered_cxx
twoproton_unfiltered::twoproton_unfiltered(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/pnfs/uboone/scratch/users/sfehlber/2020/April/Apr08/CC2p/overlay_analyzer/final.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/pnfs/uboone/scratch/users/sfehlber/2020/April/Apr08/CC2p/overlay_analyzer/final.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/pnfs/uboone/scratch/users/sfehlber/2020/April/Apr08/CC2p/overlay_analyzer/final.root:/TwoProtonAna");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

twoproton_unfiltered::~twoproton_unfiltered()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void twoproton_unfiltered::Define_Histograms(){

  //Total Histograms                                                                                                                                                                                    
  for(int i=0; i < num; i++){
    h_overlay[i] = new TH1D(Form("h_%s_overlay",total[i]),Form("h_%s_overlay",total[i]),10,0,10);
  }

  //Correlation Histograms                                                                                                                                                                             
  for(int i=0; i < num2d; i++){
    h_correlation_overlay[i] = new TH2D(Form("h_correlation_overlay_%s",total2d[i]),Form(";%s ;%s",labelx[i],labely[i]),40,0,275,40,-125,-125);
  }
  
  //Now to do the channel seperated histograms
  for(int i=0; i< number; i++){

    for(int j=0; j < number2; j++){ //my stuff                                                                                                          

      h_vtx_x[i][j]=new TH1D(Form("h_vtx_x%s%s",point[i],channel[j]),Form("h_vtx_x%s%s",point[i],channel[j]),40,0,275);
      h_vtx_y[i][j]=new TH1D(Form("h_vtx_y%s%s",point[i],channel[j]),Form("h_vtx_y%s%s",point[i],channel[j]),40,-125,125);
      h_vtx_z[i][j]=new TH1D(Form("h_vtx_z%s%s",point[i],channel[j]),Form("h_vtx_z%s%s",point[i],channel[j]),50,0,1050);
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

    }

    for(int j=0; j < number3; j++){//raquel's stuff
    
      h_vtx_x_raquel[i][j]=new TH1D(Form("h_vtx_x_raquel%s%s",point[i],channel2[j]),Form("h_vtx_x_raquel%s%s",point[i],channel2[j]),40,0,275);
      h_vtx_y_raquel[i][j]=new TH1D(Form("h_vtx_y_raquel%s%s",point[i],channel2[j]),Form("h_vtx_y_raquel%s%s",point[i],channel2[j]),40,-125,125);
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

    }
  }
}

void twoproton_unfiltered::Fill_Mine(int i, int j){

  //index i indicates at which point the histograms are being filled                                                                                                                                     
  //index j represents what channel we are filling                                                                                                                                                       
  h_vtx_x[i][j]->Fill(reco_nu_vtxx,wgt);
  h_vtx_y[i][j]->Fill(reco_nu_vtxy,wgt);
  h_vtx_z[i][j]->Fill(reco_nu_vtxz,wgt);
  h_vtx_x_mc[i][j]->Fill(mc_nu_vtxx,wgt);
  h_vtx_y_mc[i][j]->Fill(mc_nu_vtxy_sce,wgt);
  h_vtx_z_mc[i][j]->Fill(mc_nu_vtxz_sce,wgt);
  h_vtx_x_mc_sce[i][j]->Fill(mc_nu_vtxx_sce,wgt);
  h_vtx_y_mc_sce[i][j]->Fill(mc_nu_vtxy_sce,wgt);
  h_vtx_z_mc_sce[i][j]->Fill(mc_nu_vtxz_sce,wgt);
  h_q2[i][j]->Fill(mc_q2,wgt);
  h_X[i][j]->Fill(mc_X,wgt);
  h_Y[i][j]->Fill(mc_Y,wgt);
  h_Pt[i][j]->Fill(mc_Pt,wgt);

}

void twoproton_unfiltered::Fill_Raquel(int i, int j){

  //index i indicates at which point the histograms are being filled                                                                                                                                     
  //index j represents what channel we are filling                                                                                                                                                       
  h_vtx_x_raquel[i][j]->Fill(reco_nu_vtxx,wgt);
  h_vtx_y_raquel[i][j]->Fill(reco_nu_vtxy,wgt);
  h_vtx_z_raquel[i][j]->Fill(reco_nu_vtxz,wgt);
  h_vtx_x_mc_raquel[i][j]->Fill(mc_nu_vtxx,wgt);
  h_vtx_y_mc_raquel[i][j]->Fill(mc_nu_vtxy_sce,wgt);
  h_vtx_z_mc_raquel[i][j]->Fill(mc_nu_vtxz_sce,wgt);
  h_vtx_x_mc_sce_raquel[i][j]->Fill(mc_nu_vtxx_sce,wgt);
  h_vtx_y_mc_sce_raquel[i][j]->Fill(mc_nu_vtxy_sce,wgt);
  h_vtx_z_mc_sce_raquel[i][j]->Fill(mc_nu_vtxz_sce,wgt);
  h_q2_raquel[i][j]->Fill(mc_q2,wgt);
  h_X_raquel[i][j]->Fill(mc_X,wgt);
  h_Y_raquel[i][j]->Fill(mc_Y,wgt);
  h_Pt_raquel[i][j]->Fill(mc_Pt,wgt);

}

void twoproton_unfiltered::Fill_Histograms_Mine(int i){

  Fill_Mine(i,0);

  //cc0p0pi                                                                                                                                  
  if(mc_ccnc == 0 && mc_nupdg == 14 && mc_n_threshold_proton == 0 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    Fill_Mine(i,1);
    cc0p0pi[i]++;
 
    //cc1p0pi                                                                                                                                  
  } else if(mc_ccnc == 0 && mc_nupdg == 14 && mc_n_threshold_proton == 1 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    Fill_Mine(i,2);
    cc1p0pi[i]++;

    //cc2p0pi                                                                                                                                    
  } else if (mc_ccnc == 0 && mc_nupdg == 14 && mc_n_threshold_proton == 2 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    Fill_Mine(i,3);
    cc2p0pi[i]++;
           
    //ccNp0pi                                                                                                                                  
  } else if (mc_ccnc == 0 && mc_nupdg == 14 && mc_n_threshold_proton > 2 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    Fill_Mine(i,4);
    ccNp0pi[i]++;

    //ccNp1pi                                                                                                                                    
  } else if(mc_ccnc == 0 && mc_nupdg == 14 && mc_n_threshold_proton >= 0 && (mc_n_threshold_pion0 == 1 || mc_n_threshold_pionpm == 1) && fv == true){
    Fill_Mine(i,5);
    ccNp1pi[i]++;
           
    //ccNpNpi                                                                                                                                   
  } else if(mc_ccnc == 0 && mc_nupdg == 14 && mc_n_threshold_proton >= 0 && (mc_n_threshold_pion0 > 1 || mc_n_threshold_pionpm > 1) && fv == true){
    Fill_Mine(i,6);
    ccNpNpi[i]++;

    //CC NUE                                                                                                                                   
  } else if(mc_ccnc == 0 && mc_nupdg == 13 && fv == true){
    Fill_Mine(i,7);
    ccnue[i]++;
           
  //OUT OF FV                                                                                                                                
  } else if(fv == false){                                                                                                                  
    Fill_Mine(i,8);
    outfv[i]++;
  
    //NC                                                                                                                                       
  } else if(mc_ccnc == 1 && fv == true){
    Fill_Mine(i,9);
    nc[i]++;
       
    //else                                                                                                                                     
  } else{
    Fill_Mine(i,10);
    other[i]++;
  }

}

void twoproton_unfiltered::Fill_Histograms_Raquel(int i){

  Fill_Raquel(i,0);
  
  //CCQE
  if(mc_ccnc == 0 && mc_mode == 0 && fv==true){
    Fill_Raquel(i,1);
    qel[i]++;
  
    //CCCoh                                                                                                                                                                                             
  } else if(mc_ccnc == 0 && mc_mode == 3 && fv == true){
    Fill_Raquel(i,2);
    coh[i]++;

    //CCMEC                                                                                                                                                                                             
  } else if(mc_ccnc == 0 && mc_mode == 10 && fv==true){
    Fill_Raquel(i,3);
    mec[i]++;

    //CCRES                                                                                                                                                                                             
  } else if(mc_ccnc == 0 && mc_mode == 1 && fv==true){
    Fill_Raquel(i,4);
    res[i]++;
    
    //CCDIS                                                                                                                                                                                             
  } else if(mc_ccnc == 0 && mc_mode == 2 && fv==true){
    Fill_Raquel(i,5);
    dis[i]++;

    //CCNue                                                                                                                                                                                             
  } else if(mc_ccnc == 0 && mc_ccnc == 0 && mc_nupdg == abs(13) && fv ==true){
    Fill_Raquel(i,6);
    ccnue_raquel[i]++;

    //NC                                                                                                                                                                                                
  } else if(mc_ccnc == 1 && fv == true){
    Fill_Raquel(i,7);
    nc_raquel[i]++;

    //OUT OF FV                                                                                                                                                                                          
  } else if(fv == false){
    Fill_Raquel(i,8);
    outfv_raquel[i]++;

    //Other                                                                                                                                                                                             
  }else{
    Fill_Raquel(i,9);
    other_raquel[i]++;
  }
}

void twoproton_unfiltered::Write_Histograms(){
  for(int i=0; i<num; i++){
    h_overlay[i]->Write();
  }

  for(int i=0; i< num2d; i++){
    h_correlation_overlay[i]->Write();
  }

  for(int i=0; i< number; i++){
    for(int j=0; j < number2; j++){ //my stuff	
       h_vtx_x[i][j]->Write();
       h_vtx_y[i][j]->Write();
       h_vtx_z[i][j]->Write();
       h_vtx_x_mc[i][j]->Write();
       h_vtx_y_mc[i][j]->Write();
       h_vtx_z_mc[i][j]->Write();
       h_vtx_x_mc_sce[i][j]->Write();
       h_vtx_y_mc_sce[i][j]->Write();
       h_vtx_z_mc_sce[i][j]->Write();
       h_q2[i][j]->Write();
       h_X[i][j]->Write();
       h_Y[i][j]->Write();
       h_Pt[i][j]->Write();

     }

    for(int j=0; j< number3; j++){//raquel's stuff
      h_vtx_x_raquel[i][j]->Write();
      h_vtx_y_raquel[i][j]->Write();
      h_vtx_z_raquel[i][j]->Write();
      h_vtx_x_mc_raquel[i][j]->Write();
      h_vtx_y_mc_raquel[i][j]->Write();
      h_vtx_z_mc_raquel[i][j]->Write();
      h_vtx_x_mc_sce_raquel[i][j]->Write();
      h_vtx_y_mc_sce_raquel[i][j]->Write();
      h_vtx_z_mc_sce_raquel[i][j]->Write();
      h_q2_raquel[i][j]->Write();
      h_X_raquel[i][j]->Write();
      h_Y_raquel[i][j]->Write();
      h_Pt_raquel[i][j]->Write();
    }
  }
}

void twoproton_unfiltered::In_FV(float FV_edge, float x, float y, float z){

  //Just defining the x,y, and z locations of the FV
  float_t xmin = 0.0 + FV_edge;
  float_t xmax = 256.35 - FV_edge;
  float_t ymin = -116.5 + FV_edge;
  float_t ymax = 116.5 - FV_edge;
  float_t zmin = 0.0 + FV_edge;
  float_t zmax = 1036.8 - FV_edge;

  if((x <= xmin || x >= xmax) || (y <= ymin || y >= ymax) || (z <= zmin || z >= zmax)) {
    fv = false;
  } else{
    fv = true;
  }

}

Int_t twoproton_unfiltered::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t twoproton_unfiltered::LoadTree(Long64_t entry)
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

void twoproton_unfiltered::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   evtwgt_genie_pm1_funcname = 0;
   evtwgt_genie_pm1_nweight = 0;
   evtwgt_genie_pm1_weight = 0;
   evtwgt_genie_multisim_funcname = 0;
   evtwgt_genie_multisim_nweight = 0;
   evtwgt_genie_multisim_weight = 0;
   evtwgt_flux_multisim_funcname = 0;
   evtwgt_flux_multisim_nweight = 0;
   evtwgt_flux_multisim_weight = 0;
   mc_g4_mom_all = 0;
   mc_g4_mom_muon = 0;
   mc_g4_mom_proton = 0;
   mc_g4_mom_electron = 0;
   mc_g4_mom_pionpm = 0;
   mc_g4_mom_pion0 = 0;
   mc_g4_mom_neutron = 0;
   mc_g4_E = 0;
   mc_g4_p = 0;
   mc_g4_mass = 0;
   mc_g4_phi = 0;
   mc_g4_theta = 0;
   mc_g4_pdg = 0;
   mc_g4_start_x = 0;
   mc_g4_start_y = 0;
   mc_g4_start_z = 0;
   mc_g4_end_x = 0;
   mc_g4_end_y = 0;
   mc_g4_end_z = 0;
   mc_g4_start_x_sce = 0;
   mc_g4_start_y_sce = 0;
   mc_g4_start_z_sce = 0;
   mc_g4_end_x_sce = 0;
   mc_g4_end_y_sce = 0;
   mc_g4_end_z_sce = 0;
   is_from_nu_slice = 0;
   is_primary = 0;
   is_contained = 0;
   mc_pdg = 0;
   mc_primary = 0;
   mc_origin = 0;
   mc_length = 0;
   mc_start_x = 0;
   mc_start_y = 0;
   mc_start_z = 0;
   mc_end_x = 0;
   mc_end_y = 0;
   mc_end_z = 0;
   mc_start_x_sce = 0;
   mc_start_y_sce = 0;
   mc_start_z_sce = 0;
   mc_end_x_sce = 0;
   mc_end_y_sce = 0;
   mc_end_z_sce = 0;
   mc_theta = 0;
   mc_phi = 0;
   mc_ke = 0;
   mc_mom = 0;
   n_pfp = 0;
   n_trk = 0;
   id_pfp = 0;
   n_shower = 0;
   parentPDG = 0;
   trk_score = 0;
   KE_len = 0;
   dislen_ratio = 0;
   reco_q2 = 0;
   top_score = 0;
   n_daughters = 0;
   has_shower = 0;
   reco_length = 0;
   reco_start_x = 0;
   reco_start_y = 0;
   reco_start_z = 0;
   reco_end_x = 0;
   reco_end_y = 0;
   reco_end_z = 0;
   reco_theta = 0;
   reco_phi = 0;
   reco_ke = 0;
   reco_mom = 0;
   reco_mom_muon = 0;
   reco_mom_proton = 0;
   reco_mom_pion = 0;
   nhits_0 = 0;
   nhits_1 = 0;
   nhits_2 = 0;
   chi2_p_0 = 0;
   chi2_p_1 = 0;
   chi2_p_2 = 0;
   chi2_mu_0 = 0;
   chi2_mu_1 = 0;
   chi2_mu_2 = 0;
   LL3 = 0;
   LL_p_0 = 0;
   LL_p_1 = 0;
   LL_p_2 = 0;
   LL_mip_0 = 0;
   LL_mip_1 = 0;
   LL_mip_2 = 0;
   LL_mu_0 = 0;
   LL_mu_1 = 0;
   LL_mu_2 = 0;
   LL_back_p_0 = 0;
   LL_back_p_1 = 0;
   LL_back_p_2 = 0;
   LL_back_mu_0 = 0;
   LL_back_mu_1 = 0;
   LL_back_mu_2 = 0;
   TM_dedx_0 = 0;
   TM_dedx_1 = 0;
   TM_dedx_2 = 0;
   PIDA_0 = 0;
   PIDA_1 = 0;
   PIDA_2 = 0;
   start_dedx_0 = 0;
   start_dedx_1 = 0;
   start_dedx_2 = 0;
   end_dedx_0 = 0;
   end_dedx_1 = 0;
   end_dedx_2 = 0;
   ratio_dedx_0 = 0;
   ratio_dedx_1 = 0;
   ratio_dedx_2 = 0;
   avg_dedx_0 = 0;
   avg_dedx_1 = 0;
   avg_dedx_2 = 0;
   total_dedx_0 = 0;
   total_dedx_1 = 0;
   total_dedx_2 = 0;
   KE_calo_0 = 0;
   KE_calo_1 = 0;
   KE_calo_2 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("n_pfp_per_event", &n_pfp_per_event, &b_n_pfp_per_event);
   fChain->SetBranchAddress("n_trk_per_event", &n_trk_per_event, &b_n_trk_per_event);
   fChain->SetBranchAddress("n_shower_per_event", &n_shower_per_event, &b_n_shower_per_event);
   fChain->SetBranchAddress("n_neutrinos_per_event", &n_neutrinos_per_event, &b_n_neutrinos_per_event);
   fChain->SetBranchAddress("n_nu_per_event", &n_nu_per_event, &b_n_nu_per_event);
   fChain->SetBranchAddress("n_nu_pfp_per_event", &n_nu_pfp_per_event, &b_n_nu_pfp_per_event);
   fChain->SetBranchAddress("nu_PDG_per_event", &nu_PDG_per_event, &b_nu_PDG_per_event);
   fChain->SetBranchAddress("mc_ccnc", &mc_ccnc, &b_mc_ccnc);
   fChain->SetBranchAddress("mc_mode", &mc_mode, &b_mc_mode);
   fChain->SetBranchAddress("mc_interactiontype", &mc_interactiontype, &b_mc_interactiontype);
   fChain->SetBranchAddress("mc_hitnuc", &mc_hitnuc, &b_mc_hitnuc);
   fChain->SetBranchAddress("mc_q2", &mc_q2, &b_mc_q2);
   fChain->SetBranchAddress("mc_X", &mc_X, &b_mc_X);
   fChain->SetBranchAddress("mc_Y", &mc_Y, &b_mc_Y);
   fChain->SetBranchAddress("mc_Pt", &mc_Pt, &b_mc_Pt);
   fChain->SetBranchAddress("mc_nu_vtxx", &mc_nu_vtxx, &b_mc_nu_vtxx);
   fChain->SetBranchAddress("mc_nu_vtxy", &mc_nu_vtxy, &b_mc_nu_vtxy);
   fChain->SetBranchAddress("mc_nu_vtxz", &mc_nu_vtxz, &b_mc_nu_vtxz);
   fChain->SetBranchAddress("mc_nu_vtxx_sce", &mc_nu_vtxx_sce, &b_mc_nu_vtxx_sce);
   fChain->SetBranchAddress("mc_nu_vtxy_sce", &mc_nu_vtxy_sce, &b_mc_nu_vtxy_sce);
   fChain->SetBranchAddress("mc_nu_vtxz_sce", &mc_nu_vtxz_sce, &b_mc_nu_vtxz_sce);
   fChain->SetBranchAddress("mc_enu", &mc_enu, &b_mc_enu);
   fChain->SetBranchAddress("mc_wgt", &mc_wgt, &b_mc_wgt);
   fChain->SetBranchAddress("mc_wgt_cv", &mc_wgt_cv, &b_mc_wgt_cv);
   fChain->SetBranchAddress("mc_wgt_0_EtaNCEL", &mc_wgt_0_EtaNCEL, &b_mc_wgt_0_EtaNCEL);
   fChain->SetBranchAddress("mc_wgt_0_FrAbs_N", &mc_wgt_0_FrAbs_N, &b_mc_wgt_0_FrAbs_N);
   fChain->SetBranchAddress("mc_wgt_0_FrAbs_pi", &mc_wgt_0_FrAbs_pi, &b_mc_wgt_0_FrAbs_pi);
   fChain->SetBranchAddress("mc_wgt_0_FrCEx_N", &mc_wgt_0_FrCEx_N, &b_mc_wgt_0_FrCEx_N);
   fChain->SetBranchAddress("mc_wgt_0_FrCEx_pi", &mc_wgt_0_FrCEx_pi, &b_mc_wgt_0_FrCEx_pi);
   fChain->SetBranchAddress("mc_wgt_0_FrInel_N", &mc_wgt_0_FrInel_N, &b_mc_wgt_0_FrInel_N);
   fChain->SetBranchAddress("mc_wgt_0_FrInel_pi", &mc_wgt_0_FrInel_pi, &b_mc_wgt_0_FrInel_pi);
   fChain->SetBranchAddress("mc_wgt_0_FrPiProd_N", &mc_wgt_0_FrPiProd_N, &b_mc_wgt_0_FrPiProd_N);
   fChain->SetBranchAddress("mc_wgt_0_FrPiProd_pi", &mc_wgt_0_FrPiProd_pi, &b_mc_wgt_0_FrPiProd_pi);
   fChain->SetBranchAddress("mc_wgt_0_MFP_N", &mc_wgt_0_MFP_N, &b_mc_wgt_0_MFP_N);
   fChain->SetBranchAddress("mc_wgt_0_MFP_pi", &mc_wgt_0_MFP_pi, &b_mc_wgt_0_MFP_pi);
   fChain->SetBranchAddress("mc_wgt_0_MaCCQE", &mc_wgt_0_MaCCQE, &b_mc_wgt_0_MaCCQE);
   fChain->SetBranchAddress("mc_wgt_0_MaCCRES", &mc_wgt_0_MaCCRES, &b_mc_wgt_0_MaCCRES);
   fChain->SetBranchAddress("mc_wgt_0_MaNCEL", &mc_wgt_0_MaNCEL, &b_mc_wgt_0_MaNCEL);
   fChain->SetBranchAddress("mc_wgt_0_MaNCRES", &mc_wgt_0_MaNCRES, &b_mc_wgt_0_MaNCRES);
   fChain->SetBranchAddress("mc_wgt_0_MvCCRES", &mc_wgt_0_MvCCRES, &b_mc_wgt_0_MvCCRES);
   fChain->SetBranchAddress("mc_wgt_0_MvNCRES", &mc_wgt_0_MvNCRES, &b_mc_wgt_0_MvNCRES);
   fChain->SetBranchAddress("mc_wgt_0_NonRESBGvnCC1pi", &mc_wgt_0_NonRESBGvnCC1pi, &b_mc_wgt_0_NonRESBGvnCC1pi);
   fChain->SetBranchAddress("mc_wgt_0_NonRESBGvnCC2pi", &mc_wgt_0_NonRESBGvnCC2pi, &b_mc_wgt_0_NonRESBGvnCC2pi);
   fChain->SetBranchAddress("mc_wgt_0_NonRESBGvnNC1pi", &mc_wgt_0_NonRESBGvnNC1pi, &b_mc_wgt_0_NonRESBGvnNC1pi);
   fChain->SetBranchAddress("mc_wgt_0_NonRESBGvnNC2pi", &mc_wgt_0_NonRESBGvnNC2pi, &b_mc_wgt_0_NonRESBGvnNC2pi);
   fChain->SetBranchAddress("mc_wgt_0_NonRESBGvpCC1pi", &mc_wgt_0_NonRESBGvpCC1pi, &b_mc_wgt_0_NonRESBGvpCC1pi);
   fChain->SetBranchAddress("mc_wgt_0_NonRESBGvpCC2pi", &mc_wgt_0_NonRESBGvpCC2pi, &b_mc_wgt_0_NonRESBGvpCC2pi);
   fChain->SetBranchAddress("mc_wgt_0_NonRESBGvpNC1pi", &mc_wgt_0_NonRESBGvpNC1pi, &b_mc_wgt_0_NonRESBGvpNC1pi);
   fChain->SetBranchAddress("mc_wgt_0_NonRESBGvpNC2pi", &mc_wgt_0_NonRESBGvpNC2pi, &b_mc_wgt_0_NonRESBGvpNC2pi);
   fChain->SetBranchAddress("mc_wgt_0_NormCCMEC", &mc_wgt_0_NormCCMEC, &b_mc_wgt_0_NormCCMEC);
   fChain->SetBranchAddress("mc_wgt_0_NormNCMEC", &mc_wgt_0_NormNCMEC, &b_mc_wgt_0_NormNCMEC);
   fChain->SetBranchAddress("mc_wgt_1_EtaNCEL", &mc_wgt_1_EtaNCEL, &b_mc_wgt_1_EtaNCEL);
   fChain->SetBranchAddress("mc_wgt_1_FrAbs_N", &mc_wgt_1_FrAbs_N, &b_mc_wgt_1_FrAbs_N);
   fChain->SetBranchAddress("mc_wgt_1_FrAbs_pi", &mc_wgt_1_FrAbs_pi, &b_mc_wgt_1_FrAbs_pi);
   fChain->SetBranchAddress("mc_wgt_1_FrCEx_N", &mc_wgt_1_FrCEx_N, &b_mc_wgt_1_FrCEx_N);
   fChain->SetBranchAddress("mc_wgt_1_FrCEx_pi", &mc_wgt_1_FrCEx_pi, &b_mc_wgt_1_FrCEx_pi);
   fChain->SetBranchAddress("mc_wgt_1_FrInel_N", &mc_wgt_1_FrInel_N, &b_mc_wgt_1_FrInel_N);
   fChain->SetBranchAddress("mc_wgt_1_FrInel_pi", &mc_wgt_1_FrInel_pi, &b_mc_wgt_1_FrInel_pi);
   fChain->SetBranchAddress("mc_wgt_1_FrPiProd_N", &mc_wgt_1_FrPiProd_N, &b_mc_wgt_1_FrPiProd_N);
   fChain->SetBranchAddress("mc_wgt_1_FrPiProd_pi", &mc_wgt_1_FrPiProd_pi, &b_mc_wgt_1_FrPiProd_pi);
   fChain->SetBranchAddress("mc_wgt_1_MFP_N", &mc_wgt_1_MFP_N, &b_mc_wgt_1_MFP_N);
   fChain->SetBranchAddress("mc_wgt_1_MFP_pi", &mc_wgt_1_MFP_pi, &b_mc_wgt_1_MFP_pi);
   fChain->SetBranchAddress("mc_wgt_1_MaCCQE", &mc_wgt_1_MaCCQE, &b_mc_wgt_1_MaCCQE);
   fChain->SetBranchAddress("mc_wgt_1_MaCCRES", &mc_wgt_1_MaCCRES, &b_mc_wgt_1_MaCCRES);
   fChain->SetBranchAddress("mc_wgt_1_MaNCEL", &mc_wgt_1_MaNCEL, &b_mc_wgt_1_MaNCEL);
   fChain->SetBranchAddress("mc_wgt_1_MaNCRES", &mc_wgt_1_MaNCRES, &b_mc_wgt_1_MaNCRES);
   fChain->SetBranchAddress("mc_wgt_1_MvCCRES", &mc_wgt_1_MvCCRES, &b_mc_wgt_1_MvCCRES);
   fChain->SetBranchAddress("mc_wgt_1_MvNCRES", &mc_wgt_1_MvNCRES, &b_mc_wgt_1_MvNCRES);
   fChain->SetBranchAddress("mc_wgt_1_NonRESBGvnCC1pi", &mc_wgt_1_NonRESBGvnCC1pi, &b_mc_wgt_1_NonRESBGvnCC1pi);
   fChain->SetBranchAddress("mc_wgt_1_NonRESBGvnCC2pi", &mc_wgt_1_NonRESBGvnCC2pi, &b_mc_wgt_1_NonRESBGvnCC2pi);
   fChain->SetBranchAddress("mc_wgt_1_NonRESBGvnNC1pi", &mc_wgt_1_NonRESBGvnNC1pi, &b_mc_wgt_1_NonRESBGvnNC1pi);
   fChain->SetBranchAddress("mc_wgt_1_NonRESBGvnNC2pi", &mc_wgt_1_NonRESBGvnNC2pi, &b_mc_wgt_1_NonRESBGvnNC2pi);
   fChain->SetBranchAddress("mc_wgt_1_NonRESBGvpCC1pi", &mc_wgt_1_NonRESBGvpCC1pi, &b_mc_wgt_1_NonRESBGvpCC1pi);
   fChain->SetBranchAddress("mc_wgt_1_NonRESBGvpCC2pi", &mc_wgt_1_NonRESBGvpCC2pi, &b_mc_wgt_1_NonRESBGvpCC2pi);
   fChain->SetBranchAddress("mc_wgt_1_NonRESBGvpNC1pi", &mc_wgt_1_NonRESBGvpNC1pi, &b_mc_wgt_1_NonRESBGvpNC1pi);
   fChain->SetBranchAddress("mc_wgt_1_NonRESBGvpNC2pi", &mc_wgt_1_NonRESBGvpNC2pi, &b_mc_wgt_1_NonRESBGvpNC2pi);
   fChain->SetBranchAddress("mc_wgt_1_NormCCMEC", &mc_wgt_1_NormCCMEC, &b_mc_wgt_1_NormCCMEC);
   fChain->SetBranchAddress("mc_wgt_1_NormNCMEC", &mc_wgt_1_NormNCMEC, &b_mc_wgt_1_NormNCMEC);
   fChain->SetBranchAddress("evtwgt_genie_pm1_nfunc", &evtwgt_genie_pm1_nfunc, &b_evtwgt_genie_pm1_nfunc);
   fChain->SetBranchAddress("evtwgt_genie_pm1_funcname", &evtwgt_genie_pm1_funcname, &b_evtwgt_genie_pm1_funcname);
   fChain->SetBranchAddress("evtwgt_genie_pm1_nweight", &evtwgt_genie_pm1_nweight, &b_evtwgt_genie_pm1_nweight);
   fChain->SetBranchAddress("evtwgt_genie_pm1_weight", &evtwgt_genie_pm1_weight, &b_evtwgt_genie_pm1_weight);
   fChain->SetBranchAddress("evtwgt_genie_multisim_nfunc", &evtwgt_genie_multisim_nfunc, &b_evtwgt_genie_multisim_nfunc);
   fChain->SetBranchAddress("evtwgt_genie_multisim_funcname", &evtwgt_genie_multisim_funcname, &b_evtwgt_genie_multisim_funcname);
   fChain->SetBranchAddress("evtwgt_genie_multisim_nweight", &evtwgt_genie_multisim_nweight, &b_evtwgt_genie_multisim_nweight);
   fChain->SetBranchAddress("evtwgt_genie_multisim_weight", &evtwgt_genie_multisim_weight, &b_evtwgt_genie_multisim_weight);
   fChain->SetBranchAddress("evtwgt_flux_multisim_nfunc", &evtwgt_flux_multisim_nfunc, &b_evtwgt_flux_multisim_nfunc);
   fChain->SetBranchAddress("evtwgt_flux_multisim_funcname", &evtwgt_flux_multisim_funcname, &b_evtwgt_flux_multisim_funcname);
   fChain->SetBranchAddress("evtwgt_flux_multisim_nweight", &evtwgt_flux_multisim_nweight, &b_evtwgt_flux_multisim_nweight);
   fChain->SetBranchAddress("evtwgt_flux_multisim_weight", &evtwgt_flux_multisim_weight, &b_evtwgt_flux_multisim_weight);
   fChain->SetBranchAddress("mc_nupdg", &mc_nupdg, &b_mc_nupdg);
   fChain->SetBranchAddress("mc_n_muon", &mc_n_muon, &b_mc_n_muon);
   fChain->SetBranchAddress("mc_n_proton", &mc_n_proton, &b_mc_n_proton);
   fChain->SetBranchAddress("mc_n_pionpm", &mc_n_pionpm, &b_mc_n_pionpm);
   fChain->SetBranchAddress("mc_n_pion0", &mc_n_pion0, &b_mc_n_pion0);
   fChain->SetBranchAddress("mc_n_electron", &mc_n_electron, &b_mc_n_electron);
   fChain->SetBranchAddress("mc_n_neutron", &mc_n_neutron, &b_mc_n_neutron);
   fChain->SetBranchAddress("mc_n_threshold_muon", &mc_n_threshold_muon, &b_mc_n_threshold_muon);
   fChain->SetBranchAddress("mc_n_threshold_proton", &mc_n_threshold_proton, &b_mc_n_threshold_proton);
   fChain->SetBranchAddress("mc_n_threshold_pionpm", &mc_n_threshold_pionpm, &b_mc_n_threshold_pionpm);
   fChain->SetBranchAddress("mc_n_threshold_pion0", &mc_n_threshold_pion0, &b_mc_n_threshold_pion0);
   fChain->SetBranchAddress("mc_n_threshold_electron", &mc_n_threshold_electron, &b_mc_n_threshold_electron);
   fChain->SetBranchAddress("mc_n_threshold_neutron", &mc_n_threshold_neutron, &b_mc_n_threshold_neutron);
   fChain->SetBranchAddress("mc_g4_mom_all", &mc_g4_mom_all, &b_mc_g4_mom_all);
   fChain->SetBranchAddress("mc_g4_mom_muon", &mc_g4_mom_muon, &b_mc_g4_mom_muon);
   fChain->SetBranchAddress("mc_g4_mom_proton", &mc_g4_mom_proton, &b_mc_g4_mom_proton);
   fChain->SetBranchAddress("mc_g4_mom_electron", &mc_g4_mom_electron, &b_mc_g4_mom_electron);
   fChain->SetBranchAddress("mc_g4_mom_pionpm", &mc_g4_mom_pionpm, &b_mc_g4_mom_pionpm);
   fChain->SetBranchAddress("mc_g4_mom_pion0", &mc_g4_mom_pion0, &b_mc_g4_mom_pion0);
   fChain->SetBranchAddress("mc_g4_mom_neutron", &mc_g4_mom_neutron, &b_mc_g4_mom_neutron);
   fChain->SetBranchAddress("mc_g4_E", &mc_g4_E, &b_mc_g4_E);
   fChain->SetBranchAddress("mc_g4_p", &mc_g4_p, &b_mc_g4_p);
   fChain->SetBranchAddress("mc_g4_mass", &mc_g4_mass, &b_mc_g4_mass);
   fChain->SetBranchAddress("mc_g4_phi", &mc_g4_phi, &b_mc_g4_phi);
   fChain->SetBranchAddress("mc_g4_theta", &mc_g4_theta, &b_mc_g4_theta);
   fChain->SetBranchAddress("mc_g4_pdg", &mc_g4_pdg, &b_mc_g4_pdg);
   fChain->SetBranchAddress("mc_g4_start_x", &mc_g4_start_x, &b_mc_g4_start_x);
   fChain->SetBranchAddress("mc_g4_start_y", &mc_g4_start_y, &b_mc_g4_start_y);
   fChain->SetBranchAddress("mc_g4_start_z", &mc_g4_start_z, &b_mc_g4_start_z);
   fChain->SetBranchAddress("mc_g4_end_x", &mc_g4_end_x, &b_mc_g4_end_x);
   fChain->SetBranchAddress("mc_g4_end_y", &mc_g4_end_y, &b_mc_g4_end_y);
   fChain->SetBranchAddress("mc_g4_end_z", &mc_g4_end_z, &b_mc_g4_end_z);
   fChain->SetBranchAddress("mc_g4_start_x_sce", &mc_g4_start_x_sce, &b_mc_g4_start_x_sce);
   fChain->SetBranchAddress("mc_g4_start_y_sce", &mc_g4_start_y_sce, &b_mc_g4_start_y_sce);
   fChain->SetBranchAddress("mc_g4_start_z_sce", &mc_g4_start_z_sce, &b_mc_g4_start_z_sce);
   fChain->SetBranchAddress("mc_g4_end_x_sce", &mc_g4_end_x_sce, &b_mc_g4_end_x_sce);
   fChain->SetBranchAddress("mc_g4_end_y_sce", &mc_g4_end_y_sce, &b_mc_g4_end_y_sce);
   fChain->SetBranchAddress("mc_g4_end_z_sce", &mc_g4_end_z_sce, &b_mc_g4_end_z_sce);
   fChain->SetBranchAddress("is_from_nu_slice", &is_from_nu_slice, &b_is_from_nu_slice);
   fChain->SetBranchAddress("is_primary", &is_primary, &b_is_primary);
   fChain->SetBranchAddress("is_contained", &is_contained, &b_is_contained);
   fChain->SetBranchAddress("mc_pdg", &mc_pdg, &b_mc_pdg);
   fChain->SetBranchAddress("mc_primary", &mc_primary, &b_mc_primary);
   fChain->SetBranchAddress("mc_origin", &mc_origin, &b_mc_origin);
   fChain->SetBranchAddress("mc_length", &mc_length, &b_mc_length);
   fChain->SetBranchAddress("mc_start_x", &mc_start_x, &b_mc_start_x);
   fChain->SetBranchAddress("mc_start_y", &mc_start_y, &b_mc_start_y);
   fChain->SetBranchAddress("mc_start_z", &mc_start_z, &b_mc_start_z);
   fChain->SetBranchAddress("mc_end_x", &mc_end_x, &b_mc_end_x);
   fChain->SetBranchAddress("mc_end_y", &mc_end_y, &b_mc_end_y);
   fChain->SetBranchAddress("mc_end_z", &mc_end_z, &b_mc_end_z);
   fChain->SetBranchAddress("mc_start_x_sce", &mc_start_x_sce, &b_mc_start_x_sce);
   fChain->SetBranchAddress("mc_start_y_sce", &mc_start_y_sce, &b_mc_start_y_sce);
   fChain->SetBranchAddress("mc_start_z_sce", &mc_start_z_sce, &b_mc_start_z_sce);
   fChain->SetBranchAddress("mc_end_x_sce", &mc_end_x_sce, &b_mc_end_x_sce);
   fChain->SetBranchAddress("mc_end_y_sce", &mc_end_y_sce, &b_mc_end_y_sce);
   fChain->SetBranchAddress("mc_end_z_sce", &mc_end_z_sce, &b_mc_end_z_sce);
   fChain->SetBranchAddress("mc_theta", &mc_theta, &b_mc_theta);
   fChain->SetBranchAddress("mc_phi", &mc_phi, &b_mc_phi);
   fChain->SetBranchAddress("mc_ke", &mc_ke, &b_mc_ke);
   fChain->SetBranchAddress("mc_mom", &mc_mom, &b_mc_mom);
   fChain->SetBranchAddress("n_pfp", &n_pfp, &b_n_pfp);
   fChain->SetBranchAddress("n_trk", &n_trk, &b_n_trk);
   fChain->SetBranchAddress("id_pfp", &id_pfp, &b_id_pfp);
   fChain->SetBranchAddress("n_shower", &n_shower, &b_n_shower);
   fChain->SetBranchAddress("parentPDG", &parentPDG, &b_parentPDG);
   fChain->SetBranchAddress("trk_score", &trk_score, &b_trk_score);
   fChain->SetBranchAddress("KE_len", &KE_len, &b_KE_len);
   fChain->SetBranchAddress("dislen_ratio", &dislen_ratio, &b_dislen_ratio);
   fChain->SetBranchAddress("reco_q2", &reco_q2, &b_reco_q2);
   fChain->SetBranchAddress("top_score", &top_score, &b_top_score);
   fChain->SetBranchAddress("n_daughters", &n_daughters, &b_n_daughters);
   fChain->SetBranchAddress("vtx_n_pfp", &vtx_n_pfp, &b_vtx_n_pfp);
   fChain->SetBranchAddress("has_shower", &has_shower, &b_has_shower);
   fChain->SetBranchAddress("reco_nu_vtxx", &reco_nu_vtxx, &b_reco_nu_vtxx);
   fChain->SetBranchAddress("reco_nu_vtxy", &reco_nu_vtxy, &b_reco_nu_vtxy);
   fChain->SetBranchAddress("reco_nu_vtxz", &reco_nu_vtxz, &b_reco_nu_vtxz);
   fChain->SetBranchAddress("reco_length", &reco_length, &b_reco_length);
   fChain->SetBranchAddress("reco_start_x", &reco_start_x, &b_reco_start_x);
   fChain->SetBranchAddress("reco_start_y", &reco_start_y, &b_reco_start_y);
   fChain->SetBranchAddress("reco_start_z", &reco_start_z, &b_reco_start_z);
   fChain->SetBranchAddress("reco_end_x", &reco_end_x, &b_reco_end_x);
   fChain->SetBranchAddress("reco_end_y", &reco_end_y, &b_reco_end_y);
   fChain->SetBranchAddress("reco_end_z", &reco_end_z, &b_reco_end_z);
   fChain->SetBranchAddress("reco_theta", &reco_theta, &b_reco_theta);
   fChain->SetBranchAddress("reco_phi", &reco_phi, &b_reco_phi);
   fChain->SetBranchAddress("reco_ke", &reco_ke, &b_reco_ke);
   fChain->SetBranchAddress("reco_mom", &reco_mom, &b_reco_mom);
   fChain->SetBranchAddress("reco_mom_muon", &reco_mom_muon, &b_reco_mom_muon);
   fChain->SetBranchAddress("reco_mom_proton", &reco_mom_proton, &b_reco_mom_proton);
   fChain->SetBranchAddress("reco_mom_pion", &reco_mom_pion, &b_reco_mom_pion);
   fChain->SetBranchAddress("nhits_0", &nhits_0, &b_nhits_0);
   fChain->SetBranchAddress("nhits_1", &nhits_1, &b_nhits_1);
   fChain->SetBranchAddress("nhits_2", &nhits_2, &b_nhits_2);
   fChain->SetBranchAddress("chi2_p_0", &chi2_p_0, &b_chi2_p_0);
   fChain->SetBranchAddress("chi2_p_1", &chi2_p_1, &b_chi2_p_1);
   fChain->SetBranchAddress("chi2_p_2", &chi2_p_2, &b_chi2_p_2);
   fChain->SetBranchAddress("chi2_mu_0", &chi2_mu_0, &b_chi2_mu_0);
   fChain->SetBranchAddress("chi2_mu_1", &chi2_mu_1, &b_chi2_mu_1);
   fChain->SetBranchAddress("chi2_mu_2", &chi2_mu_2, &b_chi2_mu_2);
   fChain->SetBranchAddress("LL3", &LL3, &b_LL3);
   fChain->SetBranchAddress("LL_p_0", &LL_p_0, &b_LL_p_0);
   fChain->SetBranchAddress("LL_p_1", &LL_p_1, &b_LL_p_1);
   fChain->SetBranchAddress("LL_p_2", &LL_p_2, &b_LL_p_2);
   fChain->SetBranchAddress("LL_mip_0", &LL_mip_0, &b_LL_mip_0);
   fChain->SetBranchAddress("LL_mip_1", &LL_mip_1, &b_LL_mip_1);
   fChain->SetBranchAddress("LL_mip_2", &LL_mip_2, &b_LL_mip_2);
   fChain->SetBranchAddress("LL_mu_0", &LL_mu_0, &b_LL_mu_0);
   fChain->SetBranchAddress("LL_mu_1", &LL_mu_1, &b_LL_mu_1);
   fChain->SetBranchAddress("LL_mu_2", &LL_mu_2, &b_LL_mu_2);
   fChain->SetBranchAddress("LL_back_p_0", &LL_back_p_0, &b_LL_back_p_0);
   fChain->SetBranchAddress("LL_back_p_1", &LL_back_p_1, &b_LL_back_p_1);
   fChain->SetBranchAddress("LL_back_p_2", &LL_back_p_2, &b_LL_back_p_2);
   fChain->SetBranchAddress("LL_back_mu_0", &LL_back_mu_0, &b_LL_back_mu_0);
   fChain->SetBranchAddress("LL_back_mu_1", &LL_back_mu_1, &b_LL_back_mu_1);
   fChain->SetBranchAddress("LL_back_mu_2", &LL_back_mu_2, &b_LL_back_mu_2);
   fChain->SetBranchAddress("TM_dedx_0", &TM_dedx_0, &b_TM_dedx_0);
   fChain->SetBranchAddress("TM_dedx_1", &TM_dedx_1, &b_TM_dedx_1);
   fChain->SetBranchAddress("TM_dedx_2", &TM_dedx_2, &b_TM_dedx_2);
   fChain->SetBranchAddress("PIDA_0", &PIDA_0, &b_PIDA_0);
   fChain->SetBranchAddress("PIDA_1", &PIDA_1, &b_PIDA_1);
   fChain->SetBranchAddress("PIDA_2", &PIDA_2, &b_PIDA_2);
   fChain->SetBranchAddress("start_dedx_0", &start_dedx_0, &b_start_dedx_0);
   fChain->SetBranchAddress("start_dedx_1", &start_dedx_1, &b_start_dedx_1);
   fChain->SetBranchAddress("start_dedx_2", &start_dedx_2, &b_start_dedx_2);
   fChain->SetBranchAddress("end_dedx_0", &end_dedx_0, &b_end_dedx_0);
   fChain->SetBranchAddress("end_dedx_1", &end_dedx_1, &b_end_dedx_1);
   fChain->SetBranchAddress("end_dedx_2", &end_dedx_2, &b_end_dedx_2);
   fChain->SetBranchAddress("ratio_dedx_0", &ratio_dedx_0, &b_ratio_dedx_0);
   fChain->SetBranchAddress("ratio_dedx_1", &ratio_dedx_1, &b_ratio_dedx_1);
   fChain->SetBranchAddress("ratio_dedx_2", &ratio_dedx_2, &b_ratio_dedx_2);
   fChain->SetBranchAddress("avg_dedx_0", &avg_dedx_0, &b_avg_dedx_0);
   fChain->SetBranchAddress("avg_dedx_1", &avg_dedx_1, &b_avg_dedx_1);
   fChain->SetBranchAddress("avg_dedx_2", &avg_dedx_2, &b_avg_dedx_2);
   fChain->SetBranchAddress("total_dedx_0", &total_dedx_0, &b_total_dedx_0);
   fChain->SetBranchAddress("total_dedx_1", &total_dedx_1, &b_total_dedx_1);
   fChain->SetBranchAddress("total_dedx_2", &total_dedx_2, &b_total_dedx_2);
   fChain->SetBranchAddress("KE_calo_0", &KE_calo_0, &b_KE_calo_0);
   fChain->SetBranchAddress("KE_calo_1", &KE_calo_1, &b_KE_calo_1);
   fChain->SetBranchAddress("KE_calo_2", &KE_calo_2, &b_KE_calo_2);
   Notify();
}

Bool_t twoproton_unfiltered::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void twoproton_unfiltered::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t twoproton_unfiltered::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef twoproton_unfiltered_cxx
