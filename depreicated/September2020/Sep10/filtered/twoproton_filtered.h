//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr  9 11:36:39 2020 by ROOT version 6.12/06
// from TTree tree/
// found on file: /uboone/data/users/sfehlber/2020/Sep2020/Sep10/overlay_filtered.root
//////////////////////////////////////////////////////////

#ifndef twoproton_filtered_h
#define twoproton_filtered_h

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

class twoproton_filtered {
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
   vector<float>   *chi2p_3D;
   vector<float>   *chi2mu_3D;
   vector<float>   *chi2pi_3D;
   vector<float>   *chi2K_3D;
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
   TBranch        *b_chi2p_3D;
   TBranch        *b_chi2mu_3D;
   TBranch        *b_chi2pi_3D;
   TBranch        *b_chi2K_3D;
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

   twoproton_filtered(TTree *tree=0);
   virtual ~twoproton_filtered();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     In_FV(float FV_edge, float x, float y, float z);
   virtual void     Check_Containment(float FV_edge, int shortest_trk_index, int longest_trk_index, vector<float> good_trk, vector<float> good_trk_start_x, vector<float> good_trk_start_y, vector<float> good_trk_start_z, vector<float> good_trk_end_x, vector<float> good_trk_end_y, vector<float> good_trk_end_z);
   virtual void     Define_Histograms();
   virtual void     Fill_Histograms_Mine(int i);
   virtual void     Fill_Histograms_Raquel(int i);
   virtual void     Fill_Histograms_Particles(int mu, int p1, int p2);
   virtual void     Fill_Mine(int i, int j);
   virtual void     Fill_Raquel(int i, int j);
   virtual void     Fill_Particles(int j, int mu, int p1, int p2);
   virtual void     Fill_Particle_Plots(int y, float value);
   virtual void     Fill_3D_chi2(int x,float value);
   virtual void     Write_Histograms();

   //this is where you can define all your constants that you want
 private:

   //Boolean for inFV, second shortest track, and shortest track
   bool fv;
   bool second_trk;
   bool short_trk;

   //POT Weight:
   //const double wgt = 0.0289;
   const double wgt = 0.0350;

   //Per Plane Vairables
   float chi2p;
   float chi2mu;
   float dEdx;

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
   static const int  number=4; //before and after selection                                                                                                      
   static const int  number2 = 11; //categories I defined                                                                                                        
   static const int  number3 = 10; //categories raquel defined
   const char * point[number] ={"_before_selection","_after_3goodtrks","_after_containment","_after_PID"}; //this defines histograms before and after the selection                        
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

   //Chi2 Plots
   static const int num_plane = 3;
   const char* plane[num_plane] = {"_Plane_0", "_Plane_1","_Plane_2"};
   static const int num_part = 6;
   const char* particle[num_part] = {"_total","_proton","_muon","_pion","_electron","_other"};
   TH1D* h_chi2p[num_plane][num_part];
   TH1D* h_chi2mu[num_plane][num_part];
   TH1D* h_dEdx_total[num_plane][num_part];

   //3D chi2 plots
   static const int num_3D = 2;
   const char* point_3D[num_3D] = {"_before_selection","_after_selection"};
   TH1D* h_chi2p_3D[num_3D][num_part]; //3D PID
   TH1D* h_chi2mu_3D[num_3D][num_part]; //3D PID
   TH1D* h_chi2pi_3D[num_3D][num_part];

   //All the single particle plots
   static const int num_var = 4;
   const char* var[num_var] = {"_mom","_E","_theta","_phi"};
   int num_bins[num_var] = {50,50,50,10};
   double xlim_low[num_var] = {0.0,0.0,-1.5,-3.15};
   double xlim_high[num_var] = {2.0,3.0,1.5,3.15};
   double xlim_high_muon[num_var]={2.0,3.0,1.5,3.15};
   const char* xlabel[num_var] ={"P [GeV/c]","E [GeV]","cos(#theta)","#phi [Rad]"};
   TH1D* h_muon[num_var][number2];
   TH1D* h_recoil[num_var][number2];
   TH1D* h_leading[num_var][number2];
   TH1D* h_opening_angle_protons[number2];
   TH1D* h_opening_angle_mu_leading[number2];
   TH1D* h_delta_PT[number2];
   TH1D* h_delta_alphaT[number2];
   TH1D* h_delta_phiT[number2];
   TH1D* h_cos_gamma_cm_raw[number2];

   //List of all the 1D histograms
   vector<TH1*> h_list;

   //number of generated event/channel                                                                                                                                                               
   int cc0p0pi[number+1] = {0};
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

#ifdef twoproton_filtered_cxx
twoproton_filtered::twoproton_filtered(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/uboone/data/users/sfehlber/2020/Sep2020/Sep10/overlay_filtered.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/uboone/data/users/sfehlber/2020/Sep2020/Sep10/overlay_filtered.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/uboone/data/users/sfehlber/2020/Sep2020/Sep10/overlay_filtered.root:/TwoProtonAna");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

twoproton_filtered::~twoproton_filtered()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void twoproton_filtered::Define_Histograms(){

  //Total Histograms                                                                                                                                                                                    
  for(int i=0; i < num; i++){
    h_overlay[i] = new TH1D(Form("h_%s_overlay",total[i]),Form("h_%s_overlay",total[i]),10,0,10);
    h_list.push_back(h_overlay[i]);
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

      h_list.push_back(h_vtx_x[i][j]);
      h_list.push_back(h_vtx_y[i][j]);
      h_list.push_back(h_vtx_z[i][j]);
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

  for(int k=0; k < num_plane; k++){
    for(int l=0; l < num_part; l++){
      h_chi2p[k][l] = new TH1D(Form("h_chi2p%s%s",plane[k],particle[l]),Form("h_chi2p%s%s",plane[k],particle[l]),50,0,400);
      h_chi2mu[k][l] = new TH1D(Form("h_chi2mu%s%s",plane[k],particle[l]),Form("h_chi2mu%s%s",plane[k],particle[l]),50,0,120);
      h_dEdx_total[k][l] = new TH1D(Form("h_dEdx_total%s%s",plane[k],particle[l]),Form("h_dEdx_total%s%s",plane[k],particle[l]),20,0,1000);
      h_list.push_back(h_chi2p[k][l]);
      h_list.push_back(h_chi2mu[k][l]);
      h_list.push_back(h_dEdx_total[k][l]);
    }
  }

  for(int j=0; j < num_3D; j++){ 
    for(int k=0; k < num_part; k++){
      h_chi2p_3D[j][k] = new TH1D(Form("h_chi2p_3D%s%s",point_3D[j],particle[k]),Form("h_chi2p_3D%s%s",point_3D[j],particle[k]),50,0,400);
      h_chi2mu_3D[j][k] = new TH1D(Form("h_chi2mu_3D%s%s",point_3D[j],particle[k]),Form("h_chi2mu_3D%s%s",point_3D[j],particle[k]),50,0,120);
      h_chi2pi_3D[j][k] = new TH1D(Form("h_chi2pi_3D%s%s",point_3D[j],particle[k]),Form("h_chi2pi_3D%s%s",point_3D[j],particle[k]),50,0,120);
      h_list.push_back(h_chi2p_3D[j][k]);
      h_list.push_back(h_chi2mu_3D[j][k]);
      h_list.push_back(h_chi2pi_3D[j][k]);
    }
  }

  //particle specific plots
  for(int j = 0; j < num_var; j++){
    for(int k = 0; k < number2; k++){
      h_muon[j][k] = new TH1D(Form("h_muon%s%s",var[j],channel[k]),Form(" h_muon%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_muon[j]);
      h_recoil[j][k] = new TH1D(Form("h_recoil%s%s",var[j],channel[k]),Form("h_recoil%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high[j]);
      h_leading[j][k] = new TH1D(Form("h_leading%s%s",var[j],channel[k]),Form("h_leading%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high[j]);
      h_list.push_back(h_muon[j][k]);
      h_list.push_back(h_recoil[j][k]);
      h_list.push_back(h_leading[j][k]);
    }
  }

  //more particle specific plots
  for(int i = 0; i < number; i++){
    h_opening_angle_protons[i] = new TH1D(Form("h_opening_angle_protons%s",channel[i]),Form("h_opening_angle_protons%s; Opening Angle btwn Two Protons; Counts",channel[i]),35,-1.5,1.5); //50, 0, 1.5                                                                            
    h_opening_angle_mu_leading[i] = new TH1D(Form("h_opening_angle_mu_leading%s",channel[i]),Form("h_opening_angle_mu_leading%s;Opening Angle btwn Muon and Leading Proton; Counts",channel[i]),35,-1.5,1.5);
    h_delta_PT[i] = new TH1D(Form("h_delta_PT%s",channel[i]),Form("h_deltaPT%s;#delta P_{T} [GeV/c];Counts",channel[i]),10,0,1);
    h_delta_alphaT[i] = new TH1D(Form("h_delta_alphaT%s",channel[i]),Form("h_delta_alphaT%s; #delta #alpha_{T} [Deg.];Counts",channel[i]),10,-1,1); //0,180                 
    h_delta_phiT[i] = new TH1D(Form("h_delta_phiT%s",channel[i]),Form("h_delta_phiT%s; #delta #phi_{T} [Deg.];Counts",channel[i]),50,-3.15,3.15); //0,180                   
    h_cos_gamma_cm_raw[i] = new TH1D(Form("h_cos_gamma_cm_raw%s",channel[i]),Form("h_cos_gamma_cm_raw%s; cos(#gamma_{COM}); Counts",channel[i]),35,-1.5,1.5);

    h_list.push_back(h_cos_gamma_cm_raw[i]);
    h_list.push_back(h_opening_angle_protons[i]);
    h_list.push_back(h_opening_angle_mu_leading[i]);
    h_list.push_back(h_delta_PT[i]);
    h_list.push_back(h_delta_alphaT[i]);
    h_list.push_back(h_delta_phiT[i]);
  }

}

void twoproton_filtered::Fill_Mine(int i, int j){
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

void twoproton_filtered::Fill_Raquel(int i, int j){
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

void twoproton_filtered::Fill_Particles(int j, int mu, int p1, int p2){
  //first index indicates which variable is being filled: mom, energy, theta, phi                                                           
  //second index represents what channel we are filling                                                                                  
  h_muon[0][j]->Fill(reco_mom_muon->at(mu),wgt);
  h_leading[0][j]->Fill(reco_mom_proton->at(p1),wgt);
  h_recoil[0][j]->Fill(reco_mom_proton->at(p2),wgt);

  h_muon[1][j]->Fill(std::sqrt(std::pow(reco_mom_muon->at(mu),2)+std::pow(0.10565837,2)),wgt);
  h_leading[1][j]->Fill(std::sqrt(std::pow(reco_mom_proton->at(p1),2)+std::pow(0.93827208,2)),wgt);
  h_recoil[1][j]->Fill(std::sqrt(std::pow(reco_mom_proton->at(p2),2)+std::pow(0.93827208,2)),wgt);

  h_muon[2][j]->Fill(cos(reco_theta->at(mu)),wgt);
  h_leading[2][j]->Fill(cos(reco_theta->at(p1)),wgt);
  h_recoil[2][j]->Fill(cos(reco_theta->at(p2)),wgt);

  h_muon[3][j]->Fill(reco_phi->at(mu),wgt);
  h_leading[3][j]->Fill(reco_phi->at(p1),wgt);
    h_recoil[3][j]->Fill(reco_phi->at(p2),wgt);

    /*
  //Specific parameters
  ////////////////////////////////
  double muon_x = reco_mom_muon->at(mu)*cos(reco_phi->at(mu));
  double muon_y = reco_mom_muon->at(mu)*sin(reco_theta->at(mu));
  double muon_z = reco_mom_muon->at(mu)*cos(reco_theta->at(mu));
  TVector3 vMuon(muon_x,muon_y,muon_z);
  
  double p1_x = reco_mom_proton->at(p1)*cos(reco_phi->at(p1));
  double p1_y = reco_mom_proton->at(p1)*sin(reco_theta->at(p1));
  double p1_z = reco_mom_proton->at(p1)*cos(reco_theta->at(p1);
  double ELead = (std::sqrt(std::pow(reco_mom_proton->at(p1),2)+std::pow(0.93827208,2));
  TVector3 vLead(p1_x,p1_y,p1_z);
  TLorentzVector lead(vLead[0],vLead[1],vLead[2],ELead);//leading proton TLorentzVector    
		  
  double p2_x = reco_mom_proton->at(p2)*cos(reco_phi->at(p2));
  double p2_y = reco_mom_proton->at(p2)*sin(reco_theta->at(p2));
  double p2_z = reco_mom_proton->at(p2)*cos(reco_theta->at(p2));
  double ERec = (std::sqrt(std::pow(reco_mom_proton->at(p2),2)+std::pow(0.93827208,2));
  TVector3 vRec(p2_x,p2_y,p2_z);
  TLorentzVector rec(vRec[0],vRec[1],vRec[2],ERec); //Recoil proton TLorentzVector   

  //Beam Stuff
  double Eneutrino = EMuon + + 30.4;
  TVector3 vBeam(0.,0.,Eneutrino); // z-direction is defined along the neutrino direction                                                                         
  TVector3 vq = vBeam - vMuon; // Momentum transfer                                                                                                               
  TVector3 vmiss = vLead - vq; // Missing momentum         

  open_angle = ((vLead[0]*vRec[0])+(vLead[1]*vRec[1])+(vLead[2]*vRec[2]))/(vLead.Mag()*vRec.Mag()); //note this is the cos(opening angle)                             
  open_angle_mu = ((vLead[0]*vMuon[0])+(vLead[1]*vMuon[1])+(vLead[2]*vMuon[2]))/(vLead.Mag()*vMuon.Mag()); //note this is the cos(opening angle)   
  delta_pT = (vMuon + vLead).Perp();
  delta_phiT = std::acos( (-vMuon.X()*vLead.X() - vMuon.Y()*vLead.Y()) / (vMuon.XYvector().Mod() * vLead.XYvector().Mod()) );
  delta_alphaT = std::acos( (-vMuon.X()*delta_pT_vec.X()- vMuon.Y()*delta_pT_vec.Y()) / (vMuon.XYvector().Mod() * delta_pT_vec.Mod()) );

  TLorentzVector betacm(vmiss[0] + vRec[0] + vBeam[0],vmiss[1] + vRec[1] + vBeam[1], vmiss[2] + vRec[2]+ vBeam[2], En + ERec + Eneutrino); //beta for CM              
  TVector3 boost = betacm.BoostVector(); //the boost vector                                                                                                           
  lead.Boost(-boost); //boost leading proton                                                                                                                          
  rec.Boost(-boost); //boost recoil proton                                                                                                                            
  cos_gamma_cm = cos(lead.Angle(rec.Vect())); //uses Lorentz Vectors                                                                                                  

  //Some more specific plots
  h_opening_angle_protons[i]->Fill(open_angle,wgt);
  h_opening_angle_mu_leading[i]->Fill(open_angle_mu,wgt);
  h_delta_PT[i]->Fill(delta_pT,wgt);
  h_delta_alphaT[i]->Fill(cos(delta_alphaT),wgt);
  h_delta_phiT[i]->Fill(delta_phiT,wgt);
  h_cos_gamma_cm_raw[i]->Fill(cos_gamma_cm,wgt); //plots of the cos of the opening angle in the CM frame for my sakeo                             */                 

}

void twoproton_filtered::Fill_Histograms_Mine(int i){
  Fill_Mine(i,0);
  //cc0p0pi                                                                                                                                  
  if(mc_ccnc == 0 && mc_nupdg == 14 && mc_n_proton == 0 && mc_n_pion0 == 0 && mc_n_pionpm == 0 && fv == true){
    Fill_Mine(i,1);
    cc0p0pi[i]++;
    //cc1p0pi                                                                                                                                  
  } else if(mc_ccnc == 0 && mc_nupdg == 14 && mc_n_proton == 1 && mc_n_pion0 == 0 && mc_n_pionpm == 0 && fv == true){
    Fill_Mine(i,2);
    cc1p0pi[i]++;
    //cc2p0pi                                                                                                                                    
  } else if (mc_ccnc == 0 && mc_nupdg == 14 && mc_n_proton == 2 && mc_n_pion0 == 0 && mc_n_pionpm == 0 && fv == true){
    Fill_Mine(i,3);
    cc2p0pi[i]++;
    //ccNp0pi                                                                                                                                  
  } else if (mc_ccnc == 0 && mc_nupdg == 14 && mc_n_proton > 2 && mc_n_pion0 == 0 && mc_n_pionpm == 0 && fv == true){
    Fill_Mine(i,4);
    ccNp0pi[i]++;
    //ccNp1pi                                                                                                                                    
  } else if(mc_ccnc == 0 && mc_nupdg == 14 && mc_n_proton >= 0 && (mc_n_pion0 == 1 || mc_n_pionpm == 1) && fv == true){
    Fill_Mine(i,5);
    ccNp1pi[i]++;
    //ccNpNpi                                                                                                                                   
  } else if(mc_ccnc == 0 && mc_nupdg == 14 && mc_n_proton >= 0 && (mc_n_pion0 > 1 || mc_n_pionpm > 1) && fv == true){
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

void twoproton_filtered::Fill_Histograms_Raquel(int i){
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

void twoproton_filtered::Fill_Histograms_Particles(int mu, int p1, int p2){
  Fill_Particles(0, mu, p1, p2);
  //cc0p0pi                                                                                                                               
  if(mc_ccnc == 0 && mc_nupdg == 14 && mc_n_proton == 0 && mc_n_pion0 == 0 && mc_n_pionpm == 0 && fv == true){
    Fill_Particles(1, mu, p1, p2);
    cc0p0pi[number]++;
    //cc1p0pi                                                                                                                            
  } else if(mc_ccnc == 0 && mc_nupdg == 14 && mc_n_proton == 1 && mc_n_pion0 == 0 && mc_n_pionpm == 0 && fv == true){
    Fill_Particles(2, mu, p1, p2);
    cc1p0pi[number]++;
    //cc2p0pi                                                                                                                            
  } else if (mc_ccnc == 0 && mc_nupdg == 14 && mc_n_proton == 2 && mc_n_pion0 == 0 && mc_n_pionpm == 0 && fv == true){
    Fill_Particles(3, mu, p1, p2);
    cc2p0pi[number]++;
    //ccNp0pi                                                                                                                            
  } else if (mc_ccnc == 0 && mc_nupdg == 14 && mc_n_proton > 2 && mc_n_pion0 == 0 && mc_n_pionpm == 0 && fv == true){
    Fill_Particles(4, mu, p1, p2);
    ccNp0pi[number]++;
    //ccNp1pi                                                                                                                            
  } else if(mc_ccnc == 0 && mc_nupdg == 14 && mc_n_proton >= 0 && (mc_n_pion0 == 1 || mc_n_pionpm == 1) && fv == true){
    Fill_Particles(5, mu, p1, p2);
    ccNp1pi[number]++;
    //ccNpNpi                                                                                                                            
  } else if(mc_ccnc == 0 && mc_nupdg == 14 && mc_n_proton >= 0 && (mc_n_pion0 > 1 || mc_n_pionpm > 1) && fv == true){
    Fill_Particles(6, mu, p1, p2);
    ccNpNpi[number]++;
    //CC NUE                                                                                                                             
  } else if(mc_ccnc == 0 && mc_nupdg == 13 && fv == true){
    Fill_Particles(7, mu, p1, p2);
    ccnue[number]++;
  //OUT OF FV                                                                                                                            
  } else if(fv == false){                                                                                                                     Fill_Particles(8, mu, p1, p2);
    outfv[number]++;
    //NC                                                                                                                                 
  } else if(mc_ccnc == 1 && fv == true){
    Fill_Particles(9, mu, p1, p2);
    nc[number]++;    
    //else                                                                                                                               
  } else{
    Fill_Particles(10, mu, p1, p2);
    other[number]++;
  }
}

void twoproton_filtered::Write_Histograms(){

  for(int i=0; i< num2d; i++){
    h_correlation_overlay[i]->Write();
  }

  for(int i = 0; i < h_list.size(); i++){
    h_list[i]->Write();
  }
}

void twoproton_filtered::In_FV(float FV_edge, float x, float y, float z){

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

} //end of In_FV

void twoproton_filtered::Check_Containment(float FV_edge, int shortest_trk_index, int longest_trk_index, vector<float> good_trk, vector<float> good_trk_start_x, vector<float> good_trk_start_y, vector<float> good_trk_start_z, vector<float> good_trk_end_x, vector<float> good_trk_end_y, vector<float> good_trk_end_z){

  //Just defining the x,y, and z locations of the FV                                                                                      
  float_t xmin = 0.0 + FV_edge;
  float_t xmax = 256.35 - FV_edge;
  float_t ymin = -116.5 + FV_edge;
  float_t ymax = 116.5 - FV_edge;
  float_t zmin = 0.0 + FV_edge;
  float_t zmax = 1036.8 - FV_edge;
  
  for(int j=0; j < good_trk.size(); j++){  
    if(j != shortest_trk_index && j != longest_trk_index){ //second shortest track/second longest track
      if((good_trk_start_x.at(j) <= xmin || good_trk_start_x.at(j) >= xmax || good_trk_end_x.at(j) <= xmin || good_trk_end_x.at(j) >= xmax)\
	 || (good_trk_start_y.at(j) <= ymin || good_trk_start_y.at(j) >= ymax || good_trk_end_y.at(j) <= ymin || good_trk_end_y.at(j) >= ymax) 
	 || (good_trk_start_z.at(j) <= zmin || good_trk_start_z.at(j) >= zmax || good_trk_end_z.at(j) <= zmin || good_trk_end_z.at(j) >= zmax)) {
	second_trk = false;
      } else{
	second_trk = true;
      }//end of the else
      
    }else if( j == shortest_trk_index) {//shortest track
      if((good_trk_start_x.at(j) <= xmin || good_trk_start_x.at(j) >= xmax || good_trk_end_x.at(j) <= xmin || good_trk_end_x.at(j) >= xmax)
	 || (good_trk_start_y.at(j) <= ymin || good_trk_start_y.at(j) >= ymax || good_trk_end_y.at(j) <= ymin || good_trk_end_y.at(j) >= ymax)
	 || (good_trk_start_z.at(j) <= zmin || good_trk_start_z.at(j) >= zmax || good_trk_end_z.at(j) <= zmin || good_trk_end_z.at(j) >= zmax)) {
	short_trk = false;
      } else{
	short_trk = true;
      } //end of the else
    } //end of the shortest track loop
  } //end of the for loop
} //end of check_containment


void twoproton_filtered::Fill_Particle_Plots(int y, float value){
  
  h_chi2p[y][0]->Fill(chi2p,wgt);
  h_chi2mu[y][0]->Fill(chi2mu,wgt);
  h_dEdx_total[y][0]->Fill(dEdx,wgt);
	    
  if(mc_pdg->at(value) == 2212 || mc_pdg->at(value) == -2212){
    h_chi2p[y][1]->Fill(chi2p,wgt);
    h_chi2mu[y][1]->Fill(chi2mu,wgt);
    h_dEdx_total[y][1]->Fill(dEdx,wgt);
  } else if(mc_pdg->at(value) == 211 || mc_pdg->at(value) == -211 || mc_pdg->at(value) == 111) {
    h_chi2p[y][3]->Fill(chi2p,wgt);
    h_chi2mu[y][3]->Fill(chi2mu,wgt);
    h_dEdx_total[y][3]->Fill(dEdx,wgt);
  } else if(mc_pdg->at(value) == 13 || mc_pdg->at(value) == -13){
    h_chi2p[y][2]->Fill(chi2p,wgt);
    h_chi2mu[y][2]->Fill(chi2mu,wgt);
    h_dEdx_total[y][2]->Fill(dEdx,wgt);
  } else if(mc_pdg->at(value) == 11 || mc_pdg->at(value) == -11){
    h_chi2p[y][4]->Fill(chi2p,wgt);
    h_chi2mu[y][4]->Fill(chi2mu,wgt);
    h_dEdx_total[y][4]->Fill(dEdx,wgt);    
  } else {
    h_chi2p[y][5]->Fill(chi2p,wgt);
    h_chi2mu[y][5]->Fill(chi2mu,wgt);
    h_dEdx_total[y][5]->Fill(dEdx,wgt);
  }
}

void twoproton_filtered::Fill_3D_chi2(int x, float value){

  h_chi2p_3D[x][0]->Fill(chi2p_3D->at(value),wgt);
  h_chi2mu_3D[x][0]->Fill(chi2mu_3D->at(value),wgt);
  h_chi2pi_3D[x][0]->Fill(chi2pi_3D->at(value),wgt);

  if(mc_pdg->at(value) == 2212 || mc_pdg->at(value) == -2212){
    h_chi2p_3D[x][1]->Fill(chi2p_3D->at(value),wgt);
    h_chi2mu_3D[x][1]->Fill(chi2mu_3D->at(value),wgt);
    h_chi2pi_3D[x][1]->Fill(chi2pi_3D->at(value),wgt);
  } else if(mc_pdg->at(value) == 211 || mc_pdg->at(value) == -211 || mc_pdg->at(value) == 111) {
    h_chi2p_3D[x][3]->Fill(chi2p_3D->at(value),wgt);
    h_chi2mu_3D[x][3]->Fill(chi2mu_3D->at(value),wgt);
    h_chi2pi_3D[x][3]->Fill(chi2pi_3D->at(value),wgt);
  } else if(mc_pdg->at(value) == 13 || mc_pdg->at(value) == -13){
    h_chi2p_3D[x][2]->Fill(chi2p_3D->at(value),wgt);
    h_chi2mu_3D[x][2]->Fill(chi2mu_3D->at(value),wgt);
    h_chi2pi_3D[x][2]->Fill(chi2pi_3D->at(value),wgt);
  } else if(mc_pdg->at(value) == 11 || mc_pdg->at(value) == -11){
    h_chi2p_3D[x][4]->Fill(chi2p_3D->at(value),wgt);
    h_chi2mu_3D[x][4]->Fill(chi2mu_3D->at(value),wgt);
    h_chi2pi_3D[x][4]->Fill(chi2pi_3D->at(value),wgt);
  } else {
    h_chi2p_3D[x][5]->Fill(chi2p_3D->at(value),wgt);
    h_chi2mu_3D[x][5]->Fill(chi2mu_3D->at(value),wgt);
    h_chi2pi_3D[x][5]->Fill(chi2pi_3D->at(value),wgt);
    std::cout<<"PDG Value of These Tracks: "<<mc_pdg->at(value)<<std::endl;
  }
       
}

Int_t twoproton_filtered::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t twoproton_filtered::LoadTree(Long64_t entry)
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

void twoproton_filtered::Init(TTree *tree)
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
   chi2p_3D = 0;
   chi2mu_3D = 0;
   chi2pi_3D = 0;
   chi2K_3D = 0;
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
   fChain->SetBranchAddress("chi2p_3D", &chi2p_3D, &b_chi2p_3D);
   fChain->SetBranchAddress("chi2mu_3D", &chi2mu_3D, &b_chi2mu_3D);
   fChain->SetBranchAddress("chi2pi_3D", &chi2pi_3D, &b_chi2pi_3D);
   fChain->SetBranchAddress("chi2K_3D", &chi2K_3D, &b_chi2K_3D);
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

Bool_t twoproton_filtered::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void twoproton_filtered::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t twoproton_filtered::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef twoproton_filtered_cxx
