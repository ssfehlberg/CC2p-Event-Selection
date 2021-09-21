#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <fstream> 
#include "vector"
#include "variables.h"
#include "constants.h"
using namespace Constants;

class histogram_funcs
{

 public:
  virtual void Define_Histograms(const char* sample, bool Overlay); //defines histograms. works for all samples
  virtual void MC_Event_Type(int ccnc, int interaction, int nu_pdg, int mc_n_threshold_muon, int mc_n_threshold_proton, int mc_n_threshold_pion0,int mc_n_threshold_pionpm, bool fv); //fills event_type_mine and event_type_raquel. ONLY OVERLAY
  virtual void Fill_Histograms(bool Overlay,int i, TVector3 reco_nu_vertex, TVector3 true_nu_vtx, TVector3 true_nu_vtx_sce,double CosmicIP, double topological_score, double wgt); //only fills the bnb, ext, & dirt. i indicates cut
  virtual void Fill_Track_Plots(int which_cut, bool Overlay, float trk_score_v, float trk_distance_v, float trk_len_v, float trk_llr_pid_score_v, double wgt, int pdg = 0, bool contained_start = false, bool contained_end = false); //fills the track score, distance, length, and pid score
  virtual void Fill_Particles(bool Overlay, TVector3 vMuon, TLorentzVector muon, TVector3 vLead, TLorentzVector lead, TVector3 vRec, TLorentzVector rec, double wgt,double nu_e = 0.); //fills all the interesting physics plots
  virtual void Write_Histograms(); //writes histograms. works for all samples
  
  //Number of Cuts to be Applied
  //////////////////////////////////////////////////
  static const int  number = 8; //number cuts
  const char * point[number] ={"_before_selection","_after_fv","_after_three_pfps","_after_track_cut","_after_connection_cut","_after_pid","_after_containment","_after_reco_mom"}; //this defines histograms after each cut

  //Defining the two different channel types for the Overlay
  /////////////////////////////////////////

  //mine
  static const int  number2 = 11; //categories I defined                                                            
  const char * channel[number2]={"_total","_cc0p0pi","_cc1p0pi","_cc2p0pi","_ccNp0pi",
				 "_ccNp1pi","_ccNpNpi","_ccnue","_outfv","_nc","_other"}; //these are the channels I defined        

  //raquel
  static const int  number3 = 10; //categories raquel defined   
  const char * channel2[number3] = {"_total","_ccQE","_ccCOH","_ccMEC","_ccRES","_ccDIS",
				    "_ccNue","_nc","_outfv","_other"};//these are the channels that raquel defined=

  //PFP Histograms                                                                                                           
  ///////////////////////
  static const int num = 4;
  const char * total[num] = {"npfp","vtx_npfp","ntrack","nshower"};
  TH1D* h_pfp[num]; //bnb, ext, dirt
  TH1D* h_pfp_overlay[num];//overlay

  //Vertex
  ////////////////////
  TH1D* h_vtx_x[number]; //reco x: bnb, ext, dirt                                                                                   
  TH1D* h_vtx_y[number]; //reco y: bnb, ext, dirt                                                                                  
  TH1D* h_vtx_z[number]; //reco z: bnb, ext, dirt   
  TH1D* h_topological_score[number]; //bnb,ext,dirt
  TH1D* h_cosmic_impact_parameter[number];

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

  //Correlation histogramms for vertex
  static const int num2d = 3;
  const char * total2d[num2d] = {"reco","truth","truth_sce"};
  const char * labelx[num2d] = {"reco x","true x","true x + sce"};
  const char * labely[num2d] = {"reco y","true y","true y + sce"};
  TH2D *h_correlation_overlay[num2d];

  //Track related variables
  //////////////////////////////
  static const int num_track = 4;
  const char* variable[num_track] = {"_track_score","_track_vertex_distance","_track_length","_track_pid"};
  static const int track_cut = 3;
  const char* which_track_cut[track_cut] ={"_after_3_pfps","_after_track_score","_after_distance_cut"};
  static const int num_part = 11;
  const char* particle[num_part] = {"_total","_proton_contained","_proton_uncontained","_muon_contained","_muon_uncontained","_pionpm","_pion0","_electron","_gamma","_kaon","_other"};
  int num_bins_track[num_track] = {30,10,50,50};
  double xlim_low_track[num_track] = {0.0,0.0,0.0,-1.0};
  double xlim_high_track[num_track] = {1.0,10.0,50.0,1.0};
  TH1D* h_track[num_track][track_cut]; //bnb, ext, dirt
  TH1D* h_track_overlay[num_track][track_cut][num_part]; //overlay

  //All the single particle plots
  //////////////////////////////
  static const int num_var = 4;
  const char* var[num_var] = {"_mom","_E","_theta","_phi"};
  int num_bins[num_var] = {40,50,30,10}; //50 first bin, 10 last bin
  double xlim_low_muon[num_var] = {0.1,0,-1.5,-3.15}; //0.2 normally first -1.5
  double xlim_low_proton[num_var] = {0.3,0,-1.5,-3.15};
  double xlim_high_recoil[num_var] = {1.0,0.35,1.5,3.15};
  double xlim_high_leading[num_var] = {1.0,0.6,1.5,3.15}; //1.5 normally in first, 1.2
  double xlim_high_muon[num_var]={1.2,1,1.5,3.15};
  const char* xlabel[num_var] ={"P [GeV/c]","E [GeV]","cos(#theta)","#phi [Rad]"};

  TH1D* h_muon[num_var]; //bnb,ext,dirt
  TH1D* h_muon_overlay[num_var][number2];
  TH1D* h_muon_raquel[num_var][number3];

  TH1D* h_leading[num_var]; //bnb,ext,dirt
  TH1D* h_leading_overlay[num_var][number2];
  TH1D* h_leading_raquel[num_var][number3];

  TH1D* h_recoil[num_var]; //bnb,ext,dirt
  TH1D* h_recoil_overlay[num_var][number2];
  TH1D* h_recoil_raquel[num_var][number3];

  //Interesting Physics Plots
  //////////////////////////////////
  TH1D* h_opening_angle_protons_lab; //opening angle between the protons in lab frame
  TH1D* h_opening_angle_protons_lab_overlay[number2];
  TH1D* h_opening_angle_protons_lab_raquel[number3];

  TH1D* h_opening_angle_protons_com; //opening angle between the protons in com frame
  TH1D* h_opening_angle_protons_com_overlay[number2];
  TH1D* h_opening_angle_protons_com_raquel[number3];

  TH1D* h_opening_angle_mu_leading; //opening angle between the muon and lead proton
  TH1D* h_opening_angle_mu_leading_overlay[number2];
  TH1D* h_opening_angle_mu_leading_raquel[number3];

  TH1D* h_opening_angle_mu_both; //opening angle between both protons and the muon
  TH1D* h_opening_angle_mu_both_overlay[number2];
  TH1D* h_opening_angle_mu_both_raquel[number3];

  TH1D* h_delta_PT; //stvs delta PT
  TH1D* h_delta_PT_overlay[number2];  
  TH1D* h_delta_PT_raquel[number3];

  TH1D* h_delta_alphaT; //stvs delta alphaT
  TH1D* h_delta_alphaT_overlay[number2];
  TH1D* h_delta_alphaT_raquel[number3];

  TH1D* h_delta_phiT; //stvs delta phiT
  TH1D* h_delta_phiT_overlay[number2];
  TH1D* h_delta_phiT_raquel[number3];

  TH1D* h_pn; //neutron momentum estimate
  TH1D* h_pn_overlay[number2];
  TH1D* h_pn_raquel[number3];

  TH1D* h_nu_E; //neutrino energy estimate
  TH1D* h_nu_E_overlay[number2];
  TH1D* h_nu_E_raquel[number3];

  TH1D* h_mom_struck_nuc; //estimate for neutron momentum: NOT GOOD
  TH1D* h_mom_struck_nuc_overlay[number2];
  TH1D* h_mom_struck_nuc_raquel[number3];

  TH1D* h_tot_pz; //total pz of the system
  TH1D* h_tot_pz_overlay[number2];
  TH1D* h_tot_pz_raquel[number3];

  TH1D* h_tot_E; //total energy of the system
  TH1D* h_tot_E_overlay[number2];
  TH1D* h_tot_E_raquel[number3];

  TH1D* h_tot_E_minus_beam; //the total energy minus the beam
  TH1D* h_tot_E_minus_beam_overlay[number2];
  TH1D* h_tot_E_minus_beam_raquel[number3];

  TH1D* h_PT_squared; //the PT squared quantitiy in Raquel's note
  TH1D* h_PT_squared_overlay[number2];
  TH1D* h_PT_squared_raquel[number3];

  //Overlay Specific Plots
  //////////////////////
  TGraph* eff_graph = new TGraph(number); //efficiency as function of cuts                                                                                                                                                           
  TGraph* pur_graph = new TGraph(number); //efficiency as function of purity 

  TH1D* h_E_resolution_overlay[number2]; //energy resolution
  TH1D* h_E_resolution_raquel[number3];

  //Histogram Lists
  /////////////////
  vector<TH1*> h_list; //list of all the 1D histograms
  vector<TH2*> h_list_2D; //vector of all the 2D histograms

  variables variables; //variables class

  //Other parameters:    
  int event_type_mine; //type of event using my definitions
  int event_type_raquel; //type of event using raquel's definition
  double open_angle; //note this is the cos(opening angle)                                                                                                                 
  double open_angle_mu; //note this is the cos(opening angle)                                                                                                              
  double open_angle_mu_proton; //cos(opening angle)
  double delta_pT; //stv delta_pT                                                                                                                                          
  double delta_alphaT; //stv delta_alphaT                                                                                                                                  
  double delta_phiT; //stv delta_phiT                                                                                                                                      
  double cos_gamma_lab; //cos(opening angle) in lab                                                                                                                        
  double cos_gamma_cm; //cos(opening angle) in cm                                                                                                                         
  double En; //energy of struck nucleon                                                                                                                                    
  double p_struck_nuc; //momentum of the struck nucleon                                                                                                                    
  double pz_tot;

}; //end of class

//DEFINE THE HISTOGRAMS FOR BNB, EXT AND DIRT                                                                                                                                                                                                       
////////////////////////////////////////////////////   
void histogram_funcs::Define_Histograms(const char* sample, bool Overlay){

  if(Overlay == true){
  
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
	    h_muon_overlay[j][k] = new TH1D(Form("h_muon%s%s",var[j],channel[k]),Form(" h_muon%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low_muon[j],xlim_high_muon[j]);
            h_recoil_overlay[j][k] = new TH1D(Form("h_recoil%s%s",var[j],channel[k]),Form("h_recoil%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low_proton[j],xlim_high_recoil[j]);
            h_leading_overlay[j][k] = new TH1D(Form("h_leading%s%s",var[j],channel[k]),Form("h_leading%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low_proton[j],xlim_high_leading[j]);
	  }
	} else if (use_xsec_binning == false){
	  h_muon_overlay[j][k] = new TH1D(Form("h_muon%s%s",var[j],channel[k]),Form(" h_muon%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low_muon[j],xlim_high_muon[j]);
          h_recoil_overlay[j][k] = new TH1D(Form("h_recoil%s%s",var[j],channel[k]),Form("h_recoil%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low_proton[j],xlim_high_recoil[j]);
          h_leading_overlay[j][k] = new TH1D(Form("h_leading%s%s",var[j],channel[k]),Form("h_leading%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low_proton[j],xlim_high_leading[j]);
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
	    h_muon_raquel[j][k] = new TH1D(Form("h_muon_raquel%s%s",var[j],channel2[k]),Form(" h_muon_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low_muon[j],xlim_high_muon[j]);
	    h_recoil_raquel[j][k] = new TH1D(Form("h_recoil_raquel%s%s",var[j],channel2[k]),Form("h_recoil_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low_proton[j],xlim_high_recoil[j]);
	    h_leading_raquel[j][k] = new TH1D(Form("h_leading_raquel%s%s",var[j],channel2[k]),Form("h_leading_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low_proton[j],xlim_high_leading[j]);
	  }
	}else if(use_xsec_binning == false){
	  h_muon_raquel[j][k] = new TH1D(Form("h_muon_raquel%s%s",var[j],channel2[k]),Form(" h_muon_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low_muon[j],xlim_high_muon[j]);
	  h_recoil_raquel[j][k] = new TH1D(Form("h_recoil_raquel%s%s",var[j],channel2[k]),Form("h_recoil_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low_proton[j],xlim_high_recoil[j]);
	  h_leading_raquel[j][k] = new TH1D(Form("h_leading_raquel%s%s",var[j],channel2[k]),Form("h_leading_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low_proton[j],xlim_high_leading[j]);
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

    //now for all the other data products: BNB, EXT, DIRT
    /////////////////////////////////////////
  } else {

    //PFP Histograms
    for(int i=0; i < num; i++){
      h_pfp[i] = new TH1D(Form("h_%s_%s",total[i],sample),Form("h_%s_%s",total[i],sample),10,0,10);
      h_list.push_back(h_pfp[i]);
    }

    //Cut Values and Parameters
    for(int i=0; i< number; i++){
      h_vtx_x[i]=new TH1D(Form("h_vtx_x%s_%s",point[i],sample),Form("h_vtx_x%s_%s",point[i],sample),26,0,260);
      h_vtx_y[i]=new TH1D(Form("h_vtx_y%s_%s",point[i],sample),Form("h_vtx_y%s_%s",point[i],sample),24,-120,120);
      h_vtx_z[i]=new TH1D(Form("h_vtx_z%s_%s",point[i],sample),Form("h_vtx_z%s_%s",point[i],sample),104,0,1040);
      h_cosmic_impact_parameter[i] = new TH1D(Form("h_cosmic_impact_parameter%s_%s",point[i],sample),Form("h_cosmic_impact_parameter%s_%s; Cosmic Impact Distance (cm); No. Events",point[i],sample),20,0,200);
      h_topological_score[i] = new TH1D(Form("h_topological_score%s_%s",point[i],sample),Form("h_topological_score%s_%s; Topological Score; No. Events",point[i],sample),50,0.0,1.0); //30 for Wouter, 50 for steven
      h_list.push_back(h_cosmic_impact_parameter[i]);
      h_list.push_back(h_topological_score[i]);
      h_list.push_back(h_vtx_x[i]);
      h_list.push_back(h_vtx_y[i]);
      h_list.push_back(h_vtx_z[i]);
    }

    //Track Values at Various Points
    for(int j=0; j < num_track; j++){
      for(int k=0; k < track_cut; k++){
	h_track[j][k] = new TH1D(Form("h_track%s%s",variable[j],which_track_cut[k]),Form("h_track%s%s",variable[j],which_track_cut[k]),num_bins_track[j],xlim_low_track[j],xlim_high_track[j]);
	h_list.push_back(h_track[j][k]);
      }
    }

    for(int j = 0; j < num_var; j++){
      if(use_xsec_binning == true){
	if(j == 0){ //momentum
	  //muon
	  const Int_t bins_muon = 6;
	  Double_t edges_muon[bins_muon+1] = {0.1,0.2,0.3,0.5,0.7,1.3,2.5};
	  h_muon[j] = new TH1D(Form("h_muon%s_%s",var[j],sample),Form(" h_muon%s_%s ;%s; Counts",var[j],xlabel[j],sample),bins_muon,edges_muon);
	  //recoil
	  const Int_t bins_recoil = 5;
	  Double_t edges_recoil[bins_recoil+1] = {0.25,0.35,0.45,0.55,0.65,1.2};
	  h_recoil[j] = new TH1D(Form("h_recoil%s_%s",var[j],sample),Form("h_recoil%s_%s ;%s; Counts",var[j],xlabel[j],sample),bins_recoil,edges_recoil);
	  //leading
	  const Int_t bins_leading = 6;
	  Double_t edges_leading[bins_leading+1] = {0.25,0.35,0.45,0.55,0.65,0.75,1.2};
	  h_leading[j] = new TH1D(Form("h_leading%s_%s",var[j],sample),Form("h_leading%s_%s ;%s; Counts",var[j],xlabel[j],sample),bins_leading,edges_leading);

	}else if (j == 2){ //theta
	  const Int_t bins_theta = 10;
	  Double_t edges_theta[bins_theta+1] = {-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0};
	  h_muon[j] = new TH1D(Form("h_muon%s_%s",var[j],sample),Form(" h_muon%s_%s ;%s; Counts",var[j],xlabel[j],sample),bins_theta,edges_theta);
	  h_recoil[j] = new TH1D(Form("h_recoil%s_%s",var[j],sample),Form("h_recoil%s_%s ;%s; Counts",var[j],xlabel[j],sample),bins_theta,edges_theta);
	  h_leading[j] = new TH1D(Form("h_leading%s_%s",var[j],sample),Form("h_leading%s_%s ;%s; Counts",var[j],xlabel[j],sample),bins_theta,edges_theta);
	}else{ //phi and energy
	  h_muon[j] = new TH1D(Form("h_muon%s_%s",var[j],sample),Form(" h_muon%s_%s ;%s; Counts",var[j],xlabel[j],sample),num_bins[j],xlim_low_muon[j],xlim_high_muon[j]);
	  h_recoil[j] = new TH1D(Form("h_recoil%s_%s",var[j],sample),Form("h_recoil%s_%s ;%s; Counts",var[j],xlabel[j],sample),num_bins[j],xlim_low_proton[j],xlim_high_recoil[j]);
	  h_leading[j] = new TH1D(Form("h_leading%s_%s",var[j],sample),Form("h_leading%s_%s ;%s; Counts",var[j],xlabel[j],sample),num_bins[j],xlim_low_proton[j],xlim_high_leading[j]);
	} //end of phi
      }else if (use_xsec_binning == false){
	h_muon[j] = new TH1D(Form("h_muon%s_%s",var[j],sample),Form(" h_muon%s_%s ;%s; Counts",var[j],xlabel[j],sample),num_bins[j],xlim_low_muon[j],xlim_high_muon[j]);
	h_recoil[j] = new TH1D(Form("h_recoil%s_%s",var[j],sample),Form("h_recoil%s_%s ;%s; Counts",var[j],xlabel[j],sample),num_bins[j],xlim_low_proton[j],xlim_high_recoil[j]);
	h_leading[j] = new TH1D(Form("h_leading%s_%s",var[j],sample),Form("h_leading%s_%s ;%s; Counts",var[j],xlabel[j],sample),num_bins[j],xlim_low_proton[j],xlim_high_leading[j]);
      } //end of false statement
      h_list.push_back(h_muon[j]);
      h_list.push_back(h_recoil[j]);
      h_list.push_back(h_leading[j]);
    } //end of loop over variables

    if(use_xsec_binning == true){
      const Int_t bins = 10;
      Double_t edges[bins+1] = {-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0};
      h_opening_angle_protons_lab = new TH1D(Form("h_opening_angle_protons_lab_%s",sample),Form("h_opening_angle_protons_lab_%s; Opening Angle btwn Two Protons; Counts",sample),bins,edges);  
      h_opening_angle_protons_com = new TH1D(Form("h_opening_angle_protons_com_%s",sample),Form("h_opening_angle_protons_com_%s;cos(#gamma_{COM});Counts",sample),bins,edges);
      h_opening_angle_mu_leading = new TH1D(Form("h_opening_angle_mu_leading_%s",sample),Form("h_opening_angle_mu_leading_%s;Opening Angle btwn Muon and Leading Proton; Counts",sample),bins,edges);
      h_opening_angle_mu_both = new TH1D(Form("h_opening_angle_mu_both_%s",sample),Form("h_opening_angle_mu_both_%s; Opening Angle btwn Muon and Total Proton Momentum; Counts",sample),bins,edges);

      const Int_t bins_stv_mom = 4;
      Double_t edges_stv_mom[bins_stv_mom+1] = {0,0.2,0.4,0.6,1.0};
      h_delta_PT = new TH1D(Form("h_delta_PT_%s",sample),Form("h_deltaPT_%s;#delta P_{T} [GeV/c];Counts",sample),bins_stv_mom,edges_stv_mom);
      
      const Int_t bins_stv_angles = 6;
      Double_t edges_stv_angles[bins_stv_angles+1] = {0,30,60,90,120,150,180};
      h_delta_alphaT = new TH1D(Form("h_delta_alphaT_%s",sample),Form("h_delta_alphaT_%s; #delta #alpha_{T} [Deg.];Counts",sample),bins_stv_angles,edges_stv_angles); 
      h_delta_phiT = new TH1D(Form("h_delta_phiT_%s",sample),Form("h_delta_phiT_%s; #delta #phi_{T} [Deg.];Counts",sample),bins_stv_angles,edges_stv_angles);

      const Int_t bins_neutron = 20;
      Double_t edges_neutron[bins_neutron+1] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
      h_pn = new TH1D(Form("h_pn_%s",sample),Form("h_pn_%s; p_{n} [GeV/c];Counts",sample),bins_neutron,edges_neutron);
  
      const Int_t bins_nuE = 6;
      Double_t edges_nuE[bins_nuE+1] = {0,0.3,0.5,0.7,0.9,1.2,4.0};
      h_nu_E = new TH1D(Form("h_nu_E_%s",sample),Form("h_nu_E_%s; Total Energy; Counts;",sample),bins_nuE,edges_nuE);

      //basic bitches
      h_mom_struck_nuc = new TH1D(Form("h_mom_struck_nuc_%s",sample),Form("h_mom_struck_nuc_%s; P_{Init}; Counts", sample),30, 0, 1);
      h_tot_pz = new TH1D(Form("h_tot_pz_%s",sample),Form("h_tot_pz_%s; P_{Z}^{Total}; Counts",sample), 20, 0, 2);
      h_tot_E = new TH1D(Form("h_tot_E_%s",sample),Form("h_tot_E_%s; Total Energy; Counts;",sample),50,0,2.5);
      h_tot_E_minus_beam = new TH1D(Form("h_tot_E_minus_beam_%s",sample),Form("h_tot_E_minus_beam_%s; Total Energy Remaining (MeV/c); Counts;",sample),100,-100,0);
      h_PT_squared = new TH1D(Form("h_PT_squared_%s",sample),Form("h_PT_squared_%s; P_{T}^{2}; Counts", sample),50,0,5);
    
    } else if (use_xsec_binning == false){
      h_opening_angle_protons_lab = new TH1D(Form("h_opening_angle_protons_lab_%s",sample),Form("h_opening_angle_protons_lab_%s; Opening Angle btwn Two Protons; Counts",sample),30,-1.5,1.5); //50, 0, 1.5                       
      h_opening_angle_protons_com = new TH1D(Form("h_opening_angle_protons_com_%s",sample),Form("h_opening_angle_protons_com_%s;cos(#gamma_{COM});Counts",sample),30,-1.5,1.5);
      h_opening_angle_mu_leading = new TH1D(Form("h_opening_angle_mu_leading_%s",sample),Form("h_opening_angle_mu_leading_%s;Opening Angle btwn Muon and Leading Proton; Counts",sample),30,-1.5,1.5);
      h_opening_angle_mu_both = new TH1D(Form("h_opening_angle_mu_both_%s",sample),Form("h_opening_angle_mu_both_%s; Opening Angle btwn Muon and Total Proton Momentum; Counts",sample),30,-1.5,1.5);
      h_delta_PT = new TH1D(Form("h_delta_PT_%s",sample),Form("h_deltaPT_%s;#delta P_{T} [GeV/c];Counts",sample),15,0,1); //normally 10 bins                                                                                      
      h_delta_alphaT = new TH1D(Form("h_delta_alphaT_%s",sample),Form("h_delta_alphaT_%s; #delta #alpha_{T} [Deg.];Counts",sample),10,0,180); //0,180                                                                            
      h_delta_phiT = new TH1D(Form("h_delta_phiT_%s",sample),Form("h_delta_phiT_%s; #delta #phi_{T} [Deg.];Counts",sample),10,0,180); //0,180                                                                                     
      h_pn = new TH1D(Form("h_pn_%s",sample),Form("h_pn_%s; p_{n} [GeV/c];Counts",sample),50,0,1.0);
      h_nu_E = new TH1D(Form("h_nu_E_%s",sample),Form("h_nu_E_%s; Total Energy; Counts;",sample),50,0,2.5);
      h_mom_struck_nuc = new TH1D(Form("h_mom_struck_nuc_%s",sample),Form("h_mom_struck_nuc_%s; P_{Init}; Counts", sample),30, 0, 1);
      h_tot_pz = new TH1D(Form("h_tot_pz_%s",sample),Form("h_tot_pz_%s; P_{Z}^{Total}; Counts",sample), 20, 0, 2);
      h_tot_E = new TH1D(Form("h_tot_E_%s",sample),Form("h_tot_E_%s; Total Energy; Counts;",sample),50,0,2.5);
      h_tot_E_minus_beam = new TH1D(Form("h_tot_E_minus_beam_%s",sample),Form("h_tot_E_minus_beam_%s; Total Energy Remaining (MeV/c); Counts;",sample),100,-100,0);
      h_PT_squared = new TH1D(Form("h_PT_squared_%s",sample),Form("h_PT_squared_%s; P_{T}^{2}; Counts", sample),50,0,5);
    } //end loop over false statemment

    h_list.push_back(h_PT_squared);
    h_list.push_back(h_opening_angle_mu_both);
    h_list.push_back(h_nu_E);
    h_list.push_back(h_tot_E);
    h_list.push_back(h_tot_E_minus_beam);
    h_list.push_back(h_opening_angle_protons_lab);
    h_list.push_back(h_opening_angle_protons_com);
    h_list.push_back(h_opening_angle_mu_leading);
    h_list.push_back(h_delta_PT);
    h_list.push_back(h_delta_alphaT);
    h_list.push_back(h_delta_phiT);
    h_list.push_back(h_pn);
    h_list.push_back(h_mom_struck_nuc);
    h_list.push_back(h_tot_pz);
    
    //make sure to handle the weights correctly
    for (int i = 0; i < h_list.size(); i++){
      h_list[i]->Sumw2();
    }
    //for(int i = 0; i < h_list_2D.size(); i++){
    // h_list_2D[i]->Sumw2();
    // }
  } //end of else other data products

} //end of define histograms


//Function to determine what type of mC event this event is
void histogram_funcs::MC_Event_Type(int ccnc, int interaction, int nu_pdg, int mc_n_threshold_muon, int mc_n_threshold_proton, int mc_n_threshold_pion0,int mc_n_threshold_pionpm, bool fv){

  //My Definitions
  /////////////////
  event_type_mine = 0; //total

  //cc0p0pi                                                                                                                                  
  if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 0 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    event_type_mine = 1;
    //cc1p0pi                                                                                                                                  
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 1 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    event_type_mine = 1;

    //cc2p0pi                                                                                                           
  } else if (ccnc == 0 && nu_pdg == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 2 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    event_type_mine = 3;

    //ccNp0pi                                                                                                                                  
  } else if (ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton > 2 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    event_type_mine = 4;

    //ccNp1pi                                                                                                                                   
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton >= 0 && (mc_n_threshold_pion0 == 1 || mc_n_threshold_pionpm == 1) && fv == true){
    event_type_mine = 5;

    //ccNpNpi                                                                                                                                   
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton >= 0 && (mc_n_threshold_pion0 > 1 || mc_n_threshold_pionpm > 1) && fv == true){
    event_type_mine = 6;

    //CC NUE                                                                                                                                   
  } else if(ccnc == 0 && abs(nu_pdg) == 12 && fv == true){
    event_type_mine = 7;    

    //OUT OF FV                                                                                                                                
  } else if(fv == false){                                                                                                                  
    event_type_mine = 8;    

    //NC                                                                                                                                       
  } else if(ccnc == 1 && fv == true){
    event_type_mine = 9;   

    //else                                                                                                                                     
  } else{
    event_type_mine = 10;
  }

  //Raquel's Definnitions
  ///////////////////////
  event_type_raquel = 0; //total

  //CCQE
  if(ccnc == 0 && interaction == 0 && abs(nu_pdg) == 14 && fv==true){
    event_type_raquel = 1;
 
    //CCCoh                                                                                                                                                                                             
  } else if(ccnc == 0 && interaction == 3 && abs(nu_pdg) == 14 && fv == true){
    event_type_raquel = 2;

    //CCMEC                                                                                                                                                                                             
  } else if(ccnc == 0 && interaction == 10 && abs(nu_pdg) == 14 && fv==true){
    event_type_raquel = 3;

    //CCRES                                                                                                                                                                                             
  } else if(ccnc == 0 && interaction == 1 && abs(nu_pdg) == 14 && fv==true){
    event_type_raquel = 4;

    //CCDIS                                                                                                                                                                                             
  } else if(ccnc == 0 && interaction == 2 && abs(nu_pdg) == 14 && fv==true){
    event_type_raquel = 5;

    //CCNue                                                                                                                                                                                             
  } else if(ccnc == 0  && abs(nu_pdg) == 12 && fv ==true){
    event_type_raquel = 6;

    //NC                                                                                                                                                                                                
  } else if(ccnc == 1 && fv == true){
    event_type_raquel = 7;

    //OUT OF FV                                                                                                                                                                                          
  } else if(fv == false){
    event_type_raquel = 8;

    //Other                                                                                                                                                                                             
  }else{
    event_type_raquel = 9;
  } 


} //end of which MC event

//Fills the vertex, cosmicIP, and topo score.
// i = which cut.
// wgt = pot_wgt. this is mc_wgt*pot_wgt for Dirt and Overlay
///////////////////////////////////////////////////////////////////////////
void histogram_funcs::Fill_Histograms(bool Overlay,int i, TVector3 reco_nu_vertex, TVector3 true_nu_vtx, TVector3 true_nu_vtx_sce,double CosmicIP, double topological_score, double wgt){ // which cut, reco vertex, wgt to apply

  //BNB, EXT, Dirt
  if(Overlay == false){
    h_vtx_x[i]->Fill(reco_nu_vertex[0],wgt);
    h_vtx_y[i]->Fill(reco_nu_vertex[1],wgt);
    h_vtx_z[i]->Fill(reco_nu_vertex[2],wgt);
    h_topological_score[i]->Fill(topological_score,wgt);
    h_cosmic_impact_parameter[i]->Fill(CosmicIP,wgt);

  } else {

    //index i indicates at which point the histograms are being filled 
    //index j represents what channel we are filling                                                                
    h_vtx_x_overlay[i][event_type_mine]->Fill(reco_nu_vertex[0],wgt);
    h_vtx_y_overlay[i][event_type_mine]->Fill(reco_nu_vertex[1],wgt);
    h_vtx_z_overlay[i][event_type_mine]->Fill(reco_nu_vertex[2],wgt);

    h_vtx_x_mc[i][event_type_mine]->Fill(true_nu_vtx[0],wgt);
    h_vtx_y_mc[i][event_type_mine]->Fill(true_nu_vtx[1],wgt);
    h_vtx_z_mc[i][event_type_mine]->Fill(true_nu_vtx[2],wgt);

    h_vtx_x_mc_sce[i][event_type_mine]->Fill(true_nu_vtx_sce[0],wgt);
    h_vtx_y_mc_sce[i][event_type_mine]->Fill(true_nu_vtx_sce[1],wgt);
    h_vtx_z_mc_sce[i][event_type_mine]->Fill(true_nu_vtx_sce[2],wgt);

    h_topological_score_overlay[i][event_type_mine]->Fill(topological_score,wgt); //remember to add POT weight
    h_cosmic_impact_parameter_overlay[i][event_type_mine]->Fill(CosmicIP,wgt); //remember to add POT weight
    //h_q2[i][event_type_mine]->Fill(mc_q2,wgt);
    //h_X[i][event_type_mine]->Fill(mc_X,wgt);
    //h_Y[i][event_type_mine]->Fill(mc_Y,wgt);
    //h_Pt[i][event_type_mine]->Fill(mc_Pt,wgt);

    h_vtx_x_raquel[i][event_type_raquel]->Fill(reco_nu_vertex[0],wgt);
    h_vtx_y_raquel[i][event_type_raquel]->Fill(reco_nu_vertex[1],wgt);
    h_vtx_z_raquel[i][event_type_raquel]->Fill(reco_nu_vertex[2],wgt);

    h_vtx_x_mc_raquel[i][event_type_raquel]->Fill(true_nu_vtx[0],wgt);
    h_vtx_y_mc_raquel[i][event_type_raquel]->Fill(true_nu_vtx[1],wgt);
    h_vtx_z_mc_raquel[i][event_type_raquel]->Fill(true_nu_vtx[2],wgt);

    h_vtx_x_mc_sce_raquel[i][event_type_raquel]->Fill(true_nu_vtx_sce[0],wgt);
    h_vtx_y_mc_sce_raquel[i][event_type_raquel]->Fill(true_nu_vtx_sce[1],wgt);
    h_vtx_z_mc_sce_raquel[i][event_type_raquel]->Fill(true_nu_vtx_sce[2],wgt);

    h_topological_score_raquel[i][event_type_raquel]->Fill(topological_score,wgt); //remember to add POT weight
    h_cosmic_impact_parameter_raquel[i][event_type_raquel]->Fill(CosmicIP,wgt); //remember to add POT weight
    //h_q2_raquel[i][event_type_raquel]->Fill(mc_q2,wgt);
    ///h_X_raquel[i][event_type_raquel]->Fill(mc_X,wgt);
    //h_Y_raquel[i][event_type_raquel]->Fill(mc_Y,wgt);
    //h_Pt_raquel[i][event_type_raquel]->Fill(mc_Pt,wgt);

  } //end of else -> Overlay

} //end of fill histograms

//Fills the track plots, like track score, distance to the vertex, track length, and LLR.
//which_cut is an integer indicating after which cut: 0 = after 3 pfps, 1 = after track score, 2 = after vertex distance 
//trk_score, trk_distance,trk_len_v,and trk_llr_pid_score are the track score, distance to vertex, length, and PID score at a specific PFP
//wgt is the POT weight. This is POT_wgt*mc_wgt for dirt and overlay
//pdg is the pdg of the particular PFP (Overlay only)
//Contained_Start and contained_end indicate if the track is contained (Overlay only)
void histogram_funcs::Fill_Track_Plots(int which_cut, bool Overlay, float trk_score_v, float trk_distance_v, float trk_len_v, float trk_llr_pid_score_v, double wgt, int pdg = 0, bool contained_start = false, bool contained_end = false){ 

  //BNB, EXT, Dirt
  ////////////////////
  h_track[0][which_cut]->Fill(trk_score_v,wgt);
  h_track[1][which_cut]->Fill(trk_distance_v,wgt);
  h_track[2][which_cut]->Fill(trk_len_v,wgt);
  h_track[3][which_cut]->Fill(trk_llr_pid_score_v,wgt);

  //Overlay
  //////////////
  if(Overlay == true){
    h_track_overlay[0][which_cut][0]->Fill(trk_score_v,wgt); //fills the total
    h_track_overlay[1][which_cut][0]->Fill(trk_distance_v,wgt);
    h_track_overlay[2][which_cut][0]->Fill(trk_len_v,wgt);
    h_track_overlay[3][which_cut][0]->Fill(trk_llr_pid_score_v,wgt);  

    if(pdg == 2212 || pdg == -2212){
      //total_protons++;
      if(contained_start == true && contained_end == true){
	h_track_overlay[0][which_cut][1]->Fill(trk_score_v,wgt); //fills the contained protons
	h_track_overlay[1][which_cut][1]->Fill(trk_distance_v,wgt);
	h_track_overlay[2][which_cut][1]->Fill(trk_len_v,wgt);
	h_track_overlay[3][which_cut][1]->Fill(trk_llr_pid_score_v,wgt);  
	//contain++;
      } else if (contained_start == true && contained_end == false){
	h_track_overlay[0][which_cut][2]->Fill(trk_score_v,wgt); //fills the uncontained protons                                                                                                                 
	h_track_overlay[1][which_cut][2]->Fill(trk_distance_v,wgt);
	h_track_overlay[2][which_cut][2]->Fill(trk_len_v,wgt);
	h_track_overlay[3][which_cut][2]->Fill(trk_llr_pid_score_v,wgt);
	//uncontain++;
      }

    } else if(pdg == 13 || pdg == -13){
      if(contained_start == true && contained_end == true){ //contained muons
	h_track_overlay[0][which_cut][3]->Fill(trk_score_v,wgt); 
	h_track_overlay[1][which_cut][3]->Fill(trk_distance_v,wgt);
	h_track_overlay[2][which_cut][3]->Fill(trk_len_v,wgt);
	h_track_overlay[3][which_cut][3]->Fill(trk_llr_pid_score_v,wgt);  

      }else if(contained_start == true && contained_end == false){ //uncontained muons
	h_track_overlay[0][which_cut][4]->Fill(trk_score_v,wgt);
	h_track_overlay[1][which_cut][4]->Fill(trk_distance_v,wgt);
	h_track_overlay[2][which_cut][4]->Fill(trk_len_v,wgt);
	h_track_overlay[3][which_cut][4]->Fill(trk_llr_pid_score_v,wgt);
      }

    } else if(pdg == 211 || pdg == -211) {
      h_track_overlay[0][which_cut][5]->Fill(trk_score_v,wgt); //fills the pionpm
      h_track_overlay[1][which_cut][5]->Fill(trk_distance_v,wgt);
      h_track_overlay[2][which_cut][5]->Fill(trk_len_v,wgt);
      h_track_overlay[3][which_cut][5]->Fill(trk_llr_pid_score_v,wgt);  

    } else if(pdg == 111) {
      h_track_overlay[0][which_cut][6]->Fill(trk_score_v,wgt); //fills the pion0
      h_track_overlay[1][which_cut][6]->Fill(trk_distance_v,wgt);
      h_track_overlay[2][which_cut][6]->Fill(trk_len_v,wgt);
      h_track_overlay[3][which_cut][6]->Fill(trk_llr_pid_score_v,wgt);  

    } else if(pdg == 11 || pdg == -11){
      h_track_overlay[0][which_cut][7]->Fill(trk_score_v,wgt); //fills the electron
      h_track_overlay[1][which_cut][7]->Fill(trk_distance_v,wgt);
      h_track_overlay[2][which_cut][7]->Fill(trk_len_v,wgt);
      h_track_overlay[3][which_cut][7]->Fill(trk_llr_pid_score_v,wgt);  

    } else if(pdg == 22){
      h_track_overlay[0][which_cut][8]->Fill(trk_score_v,wgt); //fills the gamma
      h_track_overlay[1][which_cut][8]->Fill(trk_distance_v,wgt);
      h_track_overlay[2][which_cut][8]->Fill(trk_len_v,wgt);
      h_track_overlay[3][which_cut][8]->Fill(trk_llr_pid_score_v,wgt);  

    } else if(pdg == 321 || pdg == -321 || pdg == 311){
      h_track_overlay[0][which_cut][9]->Fill(trk_score_v,wgt); //fills the kaon
      h_track_overlay[1][which_cut][9]->Fill(trk_distance_v,wgt);
      h_track_overlay[2][which_cut][9]->Fill(trk_len_v,wgt);
      h_track_overlay[3][which_cut][9]->Fill(trk_llr_pid_score_v,wgt);  

    } else {
      if(_debug) std::cout<<"[FILL_TRACK_PLOTS] Here is the Value of the PDG in the Else Loop: "<<pdg<<std::endl;
      //other_else++;
      if(pdg == 2112){
	//neutron++;
      }
      if( pdg == 14 || pdg == -14){
	//neutrino++;
      }
      if(pdg == 0){
	//zeros++;
      }
      h_track_overlay[0][which_cut][10]->Fill(trk_score_v,wgt); //fills the else
      h_track_overlay[1][which_cut][10]->Fill(trk_distance_v,wgt);
      h_track_overlay[2][which_cut][10]->Fill(trk_len_v,wgt);
      h_track_overlay[3][which_cut][10]->Fill(trk_llr_pid_score_v,wgt);  
    }

  } //end of if overlay = true
} //end of fill track plots

//Fill the interesting physics plots
//nu_e is the truth level neutrino energy from Overlay
/////////////////////////////////////////////////////
void histogram_funcs::Fill_Particles(bool Overlay, TVector3 vMuon, TLorentzVector muon, TVector3 vLead, TLorentzVector lead, TVector3 vRec, TLorentzVector rec, double wgt, double nu_e = 0.){

  // Run the Calculate Variables function inside of variables.h. Returns following:
  // 1) vector: momenta(muon_mom,lead_mom,rec_mom);
  // 2) vector: Energies(KE_muon, TotE_muon, KE_Lead, TotE_Lead, KE_Rec, TotE_Rec);    
  // 3) vector: detector_angles(muon_theta,muon_phi,lead_theta,lead_phi,recoil_theta,recoil_phi);
  // 4) vector: opening_angles(opening_angle_protons_lab,opening_angle_protons_mu_leading,opening_angle_protons_mu_both);  // 5) double: opening_angle_protons_COM 
  // 6) vector: STVS(delta_pT,delta_alphaT,delta_phiT); 
  // 7) double: calculated_nu_E

  variables.Calculate_Variables(vMuon,vLead,vRec,add_protons);

  //muon
  double muon_mom = variables.momenta[0];
  double muon_theta = variables.detector_angles[0]; //cosine applied
  double muon_phi = variables.detector_angles[1];
  double EMuon = variables.Energies[0]; //KE
  
  //lead proton
  double lead_mom = variables.momenta[1];
  double lead_theta = variables.detector_angles[2]; //cosine applied
  double lead_phi = variables.detector_angles[3];
  double ELead = variables.Energies[2]; //KE

  //recoil proton
  double recoil_mom = variables.momenta[2];
  double recoil_theta = variables.detector_angles[4]; //cosine applied
  double recoil_phi =  variables.detector_angles[5];
  double ERec = variables.Energies[4]; //KE

  //opening angles
  double opening_angle_protons_lab = variables.opening_angles[0]; //cosine applied
  double opening_angle_protons_COM = variables.opening_angle_protons_COM; //cosine applied
  double opening_angle_protons_mu_leading = variables.opening_angles[1]; //cosine applied
  double opening_angle_protons_mu_both = variables.opening_angles[2]; //cosine applied

  //stvs
  double delta_PT = variables.stvs[0];
  double delta_alphaT = variables.stvs[1]; //degrees
  double delta_phiT = variables.stvs[2]; //degrees
  double pn = variables.stvs[3]; //GeV/c

  //neutrino energy & PT Miss
  TVector3 PT_miss(vMuon[0]+vLead[0]+vRec[0],vMuon[1]+vRec[1]+vLead[1],0);
  double Eneutrino = variables.calculated_nu_E;

  //Beam Stuff
  double E_tot = (EMuon + MASS_MUON) + ELead + ERec;
  TVector3 vBeam(0.,0.,Eneutrino); // z-direction is defined along the neutrino direction                            
  TVector3 vq = vBeam - vMuon; // Momentum transfer                                                                  
  TVector3 vmiss = vLead - vq; // Missing momentum        
  double E_tot_minus_beam = (E_tot - Eneutrino) * 1000;
  
  if(_debug) std::cout<<"[HISTOGRAM_FUNCS] Value of PT_miss magnitude: "<<PT_miss.Mag()<<std::endl;
  if(_debug) std::cout<<"[HISTOGRAM_FUNCS] Value of PT_miss magnitude2: "<<PT_miss.Mag2()<<std::endl;
  if(_debug) std::cout<<"[HISTOGRAM_FUNCS] Value of PT_miss magnitude2 divided by 2*35.37: "<<(PT_miss.Mag2())/(2.0*35.37)<<std::endl;
  if(_debug) std::cout<<"[HISTOGRAM_FUNCS] Value of Eneutrino: "<<Eneutrino<<std::endl;
  if(_debug) std::cout<<"[HISTOGRAM_FUNCS] Value of E_tot_minus_beam: "<<E_tot_minus_beam<<std::endl;

  //Struck nucleon Momentum:                                                                                            
  TVector3 vector_sum(vMuon[0] + vLead[0] + vRec[0], vMuon[1] + vLead[1] + vRec[1], vMuon[2] + vLead[2] + vRec[2]);
  TVector3 p_struck_nuc_vector(vector_sum[0], vector_sum[1], 0);
  double p_struck_nuc = p_struck_nuc_vector.Mag();
  double pz_tot = vLead[2] + vRec[2];

  //Energy of Struck Nucleon:
  double En = std::sqrt(std::pow(MASS_NEUTRON,2) + vmiss.Mag2()); //energy of struck nucleon 

  //Now to fill all the variables
  /////////////////////////////////

  //BNB, EXT, Dirt
  if(Overlay == false){

    h_muon[0]->Fill(muon_mom,wgt);
    h_leading[0]->Fill(lead_mom,wgt);
    h_recoil[0]->Fill(recoil_mom,wgt);

    h_muon[1]->Fill(EMuon,wgt);
    h_leading[1]->Fill(ELead,wgt);
    h_recoil[1]->Fill(ERec,wgt);

    h_muon[2]->Fill(muon_theta,wgt);
    h_leading[2]->Fill(lead_theta,wgt);
    h_recoil[2]->Fill(recoil_theta,wgt);

    h_muon[3]->Fill(muon_phi,wgt);
    h_leading[3]->Fill(lead_phi,wgt);
    h_recoil[3]->Fill(recoil_phi,wgt);

    h_opening_angle_protons_lab->Fill(opening_angle_protons_lab,wgt);
    h_opening_angle_protons_com->Fill(opening_angle_protons_COM,wgt);
    h_opening_angle_mu_leading->Fill(opening_angle_protons_mu_leading,wgt);
    h_opening_angle_mu_both->Fill(opening_angle_protons_mu_both,wgt);

    h_delta_PT->Fill(delta_PT,wgt);
    h_delta_alphaT->Fill(delta_alphaT,wgt);
    h_delta_phiT->Fill(delta_phiT,wgt);
    h_pn->Fill(pn,wgt);
    
    h_nu_E->Fill(Eneutrino,wgt);
    h_mom_struck_nuc->Fill(p_struck_nuc,wgt);
    h_tot_pz->Fill(pz_tot,wgt);
    h_tot_E->Fill(E_tot,wgt);
    h_tot_E_minus_beam->Fill(E_tot_minus_beam,wgt);
    h_PT_squared->Fill(PT_miss.Mag2(),wgt);

    //Overlay
  } else {
    
    //Mine
    h_muon_overlay[0][event_type_mine]->Fill(variables.momenta[0],wgt); //muon mom
    h_muon_overlay[1][event_type_mine]->Fill(variables.Energies[0],wgt); //muon energy
    h_muon_overlay[2][event_type_mine]->Fill(variables.detector_angles[0],wgt); //muon theta
    h_muon_overlay[3][event_type_mine]->Fill(variables.detector_angles[1],wgt); //muon phi

    h_leading_overlay[0][event_type_mine]->Fill(variables.momenta[1],wgt); //leading momentum
    h_leading_overlay[1][event_type_mine]->Fill(variables.Energies[2],wgt); //leading energy
    h_leading_overlay[2][event_type_mine]->Fill(variables.detector_angles[2],wgt); //leading theta
    h_leading_overlay[3][event_type_mine]->Fill(variables.detector_angles[3],wgt); //leading phi

    h_recoil_overlay[0][event_type_mine]->Fill(variables.momenta[2],wgt); //recoil momentum
    h_recoil_overlay[1][event_type_mine]->Fill(variables.Energies[4],wgt); //recoil energy
    h_recoil_overlay[2][event_type_mine]->Fill(variables.detector_angles[4],wgt); //recoil theta
    h_recoil_overlay[3][event_type_mine]->Fill(variables.detector_angles[5],wgt); //recoil phi

    h_opening_angle_protons_lab_overlay[event_type_mine]->Fill(variables.opening_angles[0],wgt);
    h_opening_angle_protons_com_overlay[event_type_mine]->Fill(variables.opening_angle_protons_COM,wgt);  
    h_opening_angle_mu_leading_overlay[event_type_mine]->Fill(variables.opening_angles[1],wgt);
    h_opening_angle_mu_both_overlay[event_type_mine]->Fill(variables.opening_angles[2],wgt);
    h_delta_PT_overlay[event_type_mine]->Fill(variables.stvs[0],wgt);
    h_delta_alphaT_overlay[event_type_mine]->Fill(variables.stvs[1],wgt);
    h_delta_phiT_overlay[event_type_mine]->Fill(variables.stvs[2],wgt);
    h_pn_overlay[event_type_mine]->Fill(variables.stvs[3],wgt);
    h_nu_E_overlay[event_type_mine]->Fill(Eneutrino,wgt);  
    h_mom_struck_nuc_overlay[event_type_mine]->Fill(p_struck_nuc,wgt);
    h_tot_pz_overlay[event_type_mine]->Fill(pz_tot,wgt);
    h_tot_E_overlay[event_type_mine]->Fill(E_tot,wgt);
    h_tot_E_minus_beam_overlay[event_type_mine]->Fill(E_tot_minus_beam,wgt);
    h_E_resolution_overlay[event_type_mine]->Fill(Eneutrino - double(nu_e) ,wgt);
    h_PT_squared_overlay[event_type_mine]->Fill(PT_miss.Mag2(),wgt);

    //Raquel
    h_muon_raquel[0][event_type_raquel]->Fill(variables.momenta[0],wgt);
    h_leading_raquel[0][event_type_raquel]->Fill(variables.momenta[1],wgt);
    h_recoil_raquel[0][event_type_raquel]->Fill(variables.momenta[2],wgt);
  
    h_muon_raquel[1][event_type_raquel]->Fill(variables.Energies[0],wgt);
    h_leading_raquel[1][event_type_raquel]->Fill(variables.Energies[2],wgt);
    h_recoil_raquel[1][event_type_raquel]->Fill(variables.Energies[4],wgt);

    h_muon_raquel[2][event_type_raquel]->Fill(variables.detector_angles[0],wgt);
    h_leading_raquel[2][event_type_raquel]->Fill(variables.detector_angles[2],wgt);
    h_recoil_raquel[2][event_type_raquel]->Fill(variables.detector_angles[4],wgt);

    h_muon_raquel[3][event_type_raquel]->Fill(variables.detector_angles[1],wgt);
    h_leading_raquel[3][event_type_raquel]->Fill(variables.detector_angles[3],wgt);
    h_recoil_raquel[3][event_type_raquel]->Fill(variables.detector_angles[5],wgt);

    h_opening_angle_protons_lab_raquel[event_type_raquel]->Fill(variables.opening_angles[0],wgt);
    h_opening_angle_protons_com_raquel[event_type_raquel]->Fill(variables.opening_angle_protons_COM,wgt);  
    h_opening_angle_mu_leading_raquel[event_type_raquel]->Fill(variables.opening_angles[1],wgt);
    h_opening_angle_mu_both_raquel[event_type_raquel]->Fill(variables.opening_angles[2],wgt);
    h_delta_PT_raquel[event_type_raquel]->Fill(variables.stvs[0],wgt);
    h_delta_alphaT_raquel[event_type_raquel]->Fill(variables.stvs[1],wgt);
    h_delta_phiT_raquel[event_type_raquel]->Fill(variables.stvs[2],wgt);
    h_pn_raquel[event_type_raquel]->Fill(variables.stvs[3],wgt);
    h_nu_E_raquel[event_type_raquel]->Fill(Eneutrino,wgt);
    h_mom_struck_nuc_raquel[event_type_raquel]->Fill(p_struck_nuc,wgt);
    h_tot_pz_raquel[event_type_raquel]->Fill(pz_tot,wgt);
    h_tot_E_raquel[event_type_raquel]->Fill(E_tot,wgt);
    h_tot_E_minus_beam_raquel[event_type_raquel]->Fill(E_tot_minus_beam,wgt);
    h_E_resolution_raquel[event_type_raquel]->Fill(Eneutrino - double(nu_e) ,wgt);
    h_PT_squared_raquel[event_type_raquel]->Fill(PT_miss.Mag2(),wgt);

  } //end of overlaay

  //have to make sure to clear the variables before you leave
  //////////////////////////////////////////////////////////
  variables.momenta.clear();
  variables.detector_angles.clear();
  variables.opening_angles.clear();
  variables.stvs.clear();
  variables.Energies.clear();

}

//Function to write Histograms.
/////////////////////////////////
void histogram_funcs::Write_Histograms(){ 
  for(int i = 0; i < h_list.size(); i++){
    h_list[i]->Write();
  }
}
