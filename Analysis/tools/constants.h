#ifndef CONSTANTS_H
#define CONSTANTS_H
namespace Constants{

 /////////////////////////////////////////
  //General variables that need to be shared:
  //////////////////////////////////////////
  
  //stuff for date and time
  time_t now = time(0);
  tm *ltm = localtime(&now);
  int Day = ltm->tm_mday;
  int Month = ltm->tm_mon + 1;
  int Year = ltm->tm_year + 1900;

  //Latex
  TLatex* t = new TLatex();

  //dummy line
  TLine* a1;
  TLine* a3;

  //POT Num and sample name:
  char const* pot_num;
  char const* sample_name;

  //number of runs
  static const int num_runs = 6;
  int run_num;

  int title_size = 60; //ROOT font size
  int title_font = 43; //ROOT font title 4 = hevelatica normal 3 = precision  
  double x_label_offset = 1.3; //RoOT how far to offset x axis label
  double y_label_offset = 1.4; //ROOT how far to offset y axis label
  
 ////////////////////////////////////////////////////////////////////
  //PLOTS THAT DON'T HAVE EITHER RAQUEL'S OR MY MONTE-CARLO BREAKDOWNS
  ////////////////////////////////////////////////////////////////////

  //Plots of variables concerning the pfp particles:
  /////////////////////////////////////////////////
  static const int num_group = 4;
  const char* group[num_group] = {"npfp","vtx_npfp","ntrack","nshower"};
  const char* titles_pfp[num_group] = {"Number of PFP","Number of PFP Attached to the Vertex","Number of Tracks","Number of Showers"};
  int ylim_pfp[num_runs][num_group] = {{6000,6000,6000,6000},
				       {6000,6000,6000,6000},
				       {6000,6000,6000,6000},
				       {6000,6000,6000,6000},
				       {6000,6000,6000,6000},
				       {6000,6000,6000,6000}}; 
  TH1D* h_overlay_pfp[num_group];
  TH1D* h_bnb_pfp[num_group];
  TH1D* h_ext_pfp[num_group];
  THStack* h_pfp[num_group];
  TCanvas* canv_pfp[num_group];
  TLegend* legend_pfp[num_group];

  //Plots of 2D Variables:
  //////////////////////////
  static const int num_group2d = 3;
  const char* group2d[num_group2d] = {"reco","truth","truth_sce"};
  TH2D* h_overlay2D[num_group2d];
  TCanvas* canv_2d[num_group2d];
  
  //Plots of Efficiency
  ///////////////////////
  static const int num_eff = 9;
  const char* eff[num_eff] = {"muon_all","muon_contained","muon_uncontained","proton_all","proton_leading","proton_recoil","pion_plus","pion_minus","pion0"};
  const char* titles_eff[num_eff] = {"True Muon Momentum (GeV/c)","True Contained #mu Momentum (GeV/c)","True Uncontained #mu Momentum","True Proton Momentum (GeV/c)","True Leading Proton Momentum (GeV/c)","True Recoil Proton Momentum (GeV/c)","True #pi^{+} Momentum (GeV/c)","True #pi^{-} Momentum (GeV/c)","True #pi^{0} Momentum (GeV/c)"};
  int ylim_eff[num_runs][num_eff] = {{150,150,150,150,150,150,20,20,20},//jan 
				     {150,150,150,150,150,150,20,20,20},//run1
				     {150,150,150,150,150,150,20,20,20},///run2
				     {150,150,150,150,150,150,20,20,20},//run3
				     {150,150,150,150,150,150,20,20,20}, //runs 1+2+# normal binning
				     {150,150,150,150,150,150,20,20,20}};//runs1+2+3 xsec binning


  double xlim_eff[num_runs][num_eff] = {{0.1,0.1,0.1,0.3,0.3,0.3,0.065,0.065,0.065},//jan
					{0.1,0.1,0.1,0.3,0.3,0.3,0.065,0.065,0.065},//run1
					{0.1,0.1,0.1,0.3,0.3,0.3,0.065,0.065,0.065},//run2
					{0.1,0.1,0.1,0.3,0.3,0.3,0.065,0.065,0.065},//run3
					{0.1,0.1,0.1,0.3,0.3,0.3,0.065,0.065,-999}, //runs 1 +2+3 normal binning
					{0.1,0.1,0.1,0.3,0.3,0.3,0.065,0.065,-999}};//runs 1+2+3  xsec binning
  
  double xlim_up_eff[num_runs][num_eff] = {{1.2,1.2,1.2,1.0,1.0,1.0,1.0,1.0,1.0},//jan
					   {1.2,1.2,1.2,1.0,1.0,1.0,1.0,1.0,1.0},//run1
					   {1.2,1.2,1.2,1.0,1.0,1.0,1.0,1.0,1.0},//run2
					   {1.2,1.2,1.2,1.0,1.0,1.0,1.0,1.0,1.0},//run3
					   {1.2,1.2,1.2,1.0,1.0,1.0,1.0,1.0,1.0}, //runs 1+2+3 normal binning
					   {1.2,1.2,1.2,1.0,1.0,1.0,1.0,1.0,1.0},};//runs 1+2+3  xsec binning
  

  
  TLine* a[num_eff];
  TLine* a_up[num_eff];
  TH1D* h_num[num_eff]; //numerator
  TH1D* h_num0[num_eff]; //numerator clone
  TH1D* h_num1[num_eff]; //numerator clone
  TH1D* h_denom[num_eff]; //denominator
  TH1D* h_denom0[num_eff]; //denominator clone
  TH1D* h_denom1[num_eff]; //denominator clone
  TCanvas* canv_eff[num_eff];
  TCanvas* canv_both[num_eff];
  TLegend* legend_eff[num_eff];
  TGraph* eff_graph; //efficiency aas function of cuts
  TGraph* pur_graph; //purity as function of cuts


  static const int num_particles_eff_plots = 3;                                                                                                                                                                                        
  const char* particles_eff_var[num_particles_eff_plots] = {"_mom","_costheta","_phi"};
  const char* particles_eff_var_titles[num_particles_eff_plots] = {"True Momentum (GeV/c)","True cos(#theta)","True #phi (Rad.)"};
  static const int num_particles_eff = 5;                                                                                                                                                                                              
  const char* particles_eff[num_particles_eff] = {"_muon_all","_muon_contained","_muon_uncontainied","_lead_proton","_recoil_proton"};
  const char* particles_eff_titles[num_particles_eff] = {"All Muons","Contained Muons","Uncontained Muons","Leading Proton","Recoil Proton"};
  TH1D* h_particle_num[num_particles_eff][num_particles_eff_plots]; //particles for                                                                                                                                                    
  TH1D* h_particle_denom[num_particles_eff][num_particles_eff_plots];     
  TCanvas* canv_particle_eff[num_particles_eff][num_particles_eff_plots];
  
  //efficiency plots of other variables to determine potential xsec candidates
  static const int num_other_eff = 8;
  const char* other_eff[num_other_eff] = {"_opening_angle_protons_lab","_opening_angle_protons_com","_opening_angle_mu_leading","_opening_angle_mu_both","_delta_PT","_delta_alphaT","_delta_phiT","_nu_E"};
  const char* other_eff_titles[num_other_eff] = {"cos(#gamma_{p_{L},p_{R}})","cos(#gamma_{1#mu2p COM}","cos(#gamma_{#mu,p_{L}})","cos(#gamma_{#mu,p_{sum})","#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)","True Neutrino Energy (GeV)"};
  TH1D* h_other_eff_num[num_other_eff];
  TH1D* h_other_eff_denom[num_other_eff];
  TCanvas* canv_other_eff[num_other_eff];
  
  //Various Track Variables for Cut Tracking
  /////////////////////////////////////////
  static const int num_track = 4;
  static const int num_particles = 11;
  static const int track_cut = 3;
  const char* variable[num_track] = {"_track_score","_track_vertex_distance","_track_length","_track_pid"};
  const char* which_track_cut[track_cut] = {"_after_3_pfps","_after_track_score","_after_distance_cut"};
  const char* particles[num_particles] = {"_total",
					  "_proton_contained",
					  "_proton_uncontained",
					  "_muon_contained",
					  "_muon_uncontained",
					  "_pionpm",
					  "_pion0",
					  "_electron",
					  "_gamma",
					  "_kaon",
					  "_other"};
  std::vector<const char*> channel_legend_chi2 = {"Total Overlay","Protons:Contained","Protons:Uncontained","Muon:Contained","Muon:Uncontained","Pion#pm","Pion0,","Electron","#gamma","Kaon","Other"};
  TH1D* h_track_overlay[num_track][track_cut][num_particles];
  std::vector<TH1D*> h_track_overlay_vec;
  TH1D* h_track_overlay0[num_track][track_cut][3]; 
  TH1D* h_track_bnb[num_track][track_cut][2];
  TH1D* h_track_ext[num_track][track_cut][3];
  TH1D* h_track_dirt[num_track][track_cut][3];
  THStack* h_track[num_track][track_cut];
  TCanvas* canv_track[num_track][track_cut];
  TLegend* legend_track[num_track][track_cut];
  TPad* pad_track[num_track][track_cut];
  TPad* pad0_track[num_track][track_cut];
  int ymin_track[num_track] = {0,0,0,0};
  int ymax_track[num_track] = {50000,60000,700,4500};
  Double_t xlim_track[num_track] = {0.8, 4.0 ,0.0 ,0.2}; 
  TLine* a_track[num_track];
  const char* titles_track[num_track] = {"Track Score","Track Vertex Distance (cm)","Track Length (cm)","Track PID"};
  const char* titles_track_cuts[track_cut] = {"After 3PFPs","After Track Score","After Distance Cut"};
  
  ///////////////////////////////////////
  //PLOTS WITH MY MONTEC-CARLO DEFINITIONS
  ////////////////////////////////////////

  //General
  /////////
  static const int num_cuts = 3; //number of applied cuts
  static const int num_channels = 11; //number of various overlay channels
  const char* cut[num_cuts] = {"_before_selection","_after_fv","_after_three_pfps"};//,"_after_track_cut","_after_coonnection_cut","_after_pid"}; 
  //const char* channel[num_channels] = {"_total", "_cc2p0pi","_ccNp1pi","_ccNp0pi","_cc1p0pi","_nc","_ccNpNpi","_cc0p0pi","_ccnue","_outfv","_other"};
  //std::vector<const char*> channel_legend = {"Total Overlay","CC2p0#pi (Signal)","CC(N>=0)p1#pi","CC(N>2)p0#pi","CC1p0#pi","NC","CC(N>=0)p(N>1)#pi","CC0p#pi","CC#nu_{e}","Out of FV","Other"};
  /*const char* channel[num_channels] = {"_total",
				       "_cc2p0pi",
				       "_ccNp1pi",
				       "_cc0p0pi",
				       "_ccNp0pi",
				       "_ccnue",
				       "_cc1p0pi",
				       "_outfv",
				       "_nc",
				       "_other",
				       "_ccNpNpi"};

  std::vector<const char*> channel_legend = {"Total Overlay",
					     "CC2p0#pi (Signal)",
					     "CC(N>=0)p(N>1)#pi",
					     "CC(N>=0)p1#pi",
					     "CC1p0#pi",
					     "CC(N>2)p0#pi",
					     "CC#nu_{e}",
					     "CC0p0#pi",
					     "Out of FV",
					     "NC",
					     "Other"};*/



  
  const char* channel[num_channels] = {"_total",
				       "_cc2p0pi",
				       "_ccNpNpi",
				       "_ccNp1pi",
				       "_cc0p0pi",
				       "_ccNp0pi",
				       "_ccnue",
				       "_cc1p0pi",
				       "_outfv",
				       "_nc",
				       "_other"};
  std::vector<const char*> channel_legend = {"Total Overlay",
					     "CC2p0#pi (Signal)",
					     "CC(N>=0)p(N>1)#pi",
					     "CC(N>=0)p1#pi",
					     "CC0p0#pi",
					     "CC(N>2)p0#pi",
					     "CC#nu_{e}",
					     "CC1p0#pi",
					     "Out of FV",
					     "NC",
					     "Other"};

  const char* cut_titles[num_cuts] = {"Before Selection", "After FV", "After PFP Cut"};//,"After Track Score Cut","After Vertex Attachment Cut","After PID"};
  int z = 0; //dumb indices
  int f = 0; //dumb indices

  //Plots of the reconstructed vertex: MC, Dirt, BNB, & EXT
  /////////////////////////////////////////////////////////
  static const int num_variables = 5; //number of various variables to plot
  const char* plots[num_variables] = {"_vtx_x","_vtx_y","_vtx_z","_cosmic_impact_parameter","_topological_score"};
  const char* titles[num_variables] = {"Reconstructed X of Vertex (cm)", "Reconstructed Y of Vertex (cm)", "Reconstructed Z of Vertex (cm)","Cosmic Impact Parameter","Topological Score"};
  int ylim[num_runs][num_cuts][num_variables] = {
						 { {1000,1000,1000,2000,3500},{1000,1000,600,600,3500},{1000,1000,600,450,3500} }, //jan
						 { {1000,1000,1000,2000,3500},{1000,1000,600,600,3500},{1000,1000,600,450,3500} }, //run1
						 { {1000,1000,1000,2000,3500},{1000,1000,600,600,3500},{1000,1000,600,450,3500} }, //run2
						 { {1000,1000,1000,2000,3500},{1000,1000,600,600,3500},{1000,1000,600,450,3500} }, //run3
						 { {80000,160000,35000,2000,3500},{35000,35000,10000,600,3500},{10000,10000,1000,450,3500} }, //runs 1+2+3 normal binning
						 { {80000,160000,35000,2000,3500},{100000,100000,100000,600,3500},{1000,1000,1000,450,3500} }//runs 1+2+3 xsec binning
                                                };

  int ymin[num_runs][num_cuts][num_variables] = {
						 { {-3,-3,-4,-4,-4}, {-4,-3,-3,-3,-3}, {-2,-2,-2,-4,-4} }, //jan
						 { {-3,-3,-4,-4,-4}, {-4,-3,-3,-3,-3}, {-2,-2,-2,-4,-4} },//run1
						 { {-3,-3,-4,-4,-4}, {-4,-3,-3,-3,-3}, {-2,-2,-2,-4,-4} }, //run2
						 { {-3,-3,-4,-4,-4}, {-4,-3,-3,-3,-3}, {-2,-2,-2,-4,-4} }, //run3
						 { {-3,-3,-4,-4,-4}, {-4,-3,-3,-3,-3}, {-2,-2,-2,-4,-4} }, //runs 1+2+3 normal binning
						  { {-3,-3,-4,-4,-4}, {-4,-3,-3,-3,-3}, {-2,-2,-2,-4,-4} }//runs 1+2+3 xsec binning
                                                };
  Double_t xmin_vtx[num_variables] = {10.0,-110,10.0};
  Double_t xmax_vtx[num_variables] = {250,110,1030};
  bool plot_total = true; //sanity check on total overlay
  bool plot_ccnue = true; //do you want to plot ccnue?
  bool flip_legend = false;
  std::vector<TH1D*> h_overlay_vec; //overlay plots
  TH1D* h_overlay[num_cuts][num_variables][num_channels];
  TH1D* h_overlay0[num_cuts][num_variables][3]; //copies of the overlay plots for statistics
  TH1D* h_dirt[num_cuts][num_variables][3];
  TH1D* h_bnb[num_cuts][num_variables][2]; //two copies
  TH1D* h_ext[num_cuts][num_variables][3]; //two copies
  THStack* h[num_cuts][num_variables];
  TCanvas* canv[num_cuts][num_variables];
  TLegend* legend[num_cuts][num_variables];
  TPad* pad[num_cuts][num_variables];
  TPad* pad0[num_cuts][num_variables];

 //Plots of reconstructed vertex: Truth Only
 /////////////////////////////////////////////
 static const int num_truth = 9;//10;
 const char* truth[num_truth] = {"_vtx_x_mc","_vtx_y_mc","_vtx_z_mc","_vtx_x_mc_sce","_vtx_y_mc_sce","_vtx_z_mc_sce","_q2","_X","_Y"};
 const char* truth_titles[num_truth] = {"True Vertex X (cm)","True Vertex Y (cm)","True Vertex Z (cm)",
					"True Vertex w/ SCE X (cm)","True Vertex w/ SCE Y (cm)","True Vertex w/ SCE Z (cm)",
					"True Q^{2} (GeV^{2})","True Bjorken X","True Bjorken Y"};
 int ylim_truth[num_runs][num_cuts][num_truth] = {
						  { {2000,2000,2000,2000,2000,2000,2000,2000,2000},{6000,60,60,60,60,60,200,200,200},{6000,60,60,60,60,60,200,200,200} }, //jan
						  { {2000,2000,2000,2000,2000,2000,2000,2000,2000},{6000,60,60,60,60,60,200,200,200},{6000,60,60,60,60,60,200,200,200} }, //run1
						  { {2000,2000,2000,2000,2000,2000,2000,2000,2000},{6000,60,60,60,60,60,200,200,200},{6000,60,60,60,60,60,200,200,200} }, //run2
						  { {2000,2000,2000,2000,2000,2000,2000,2000,2000},{6000,60,60,60,60,60,200,200,200},{6000,60,60,60,60,60,200,200,200} }, //run3
						  { {2000,2000,2000,2000,2000,2000,2000,2000,2000},{6000,60,60,60,60,60,200,200,200},{6000,60,60,60,60,60,200,200,200} }, //runs 1+2+3 normal binning
						  { {2000,2000,2000,2000,2000,2000,2000,2000,2000},{6000,60,60,60,60,60,200,200,200},{6000,60,60,60,60,60,200,200,200} }//runs 1+2+3 xsec binning
                                                 };
 bool plot_total_truth = true;
 bool plot_ccnue_truth = true;
 TH1D* h_mc[num_cuts][num_truth][num_channels];
 TCanvas* canv_truth[num_cuts][num_truth];
 THStack* h_truth[num_cuts][num_truth];
 TLegend* legend_truth[num_cuts][num_truth];

 //Plots of the individual particle quantities
 ////////////////////////////////////////////
  static const int num_var = 3;
  const char* var[num_var] = {"_mom","_theta","_phi"};//,"_E"};
  const char* titles_var[num_var] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};//,"Energy(GeV)"};
  
 //Muon
  TH1D* h_muon_overlay[num_var][num_channels]; //actual overlay
  std::vector<TH1D*> h_muon_overlay_vec;
  TH1D* h_muon_overlay0[num_var][3]; //copies for statistics
  TH1D* h_muon_stat[num_var]; //statistical uncertainty
  TH1D* h_muon_systematic[num_var]; //systematic uncertainty
  TH1D* h_muon_bnb[num_var][2];
  TH1D* h_muon_ext[num_var][3];
  TH1D* h_muon_dirt[num_var][3];
  THStack* h_muon[num_var];
  TCanvas* canv_muon[num_var];
  TLegend* legend_muon[num_var];
  TPad* pad_muon[num_var];
  TPad* pad0_muon[num_var];
  std::vector<bool> flip_muon = {false,true,false};
  int muon_ylim[num_runs][num_var] = {{40,160,100}, //jan
				      {160,450,275},//run1
				      {250,650,450},//run2
				      {300,1000,600},//run3
				      {300,1100,800}, //runs 1+2+3 normal binning
				      {2700,2500,1450}}; //runs 1+2+3 xsec binning
  
  int muon_ymin[num_runs][num_var] = {{-2,-2,-2},//jan
				      {-2,-2,-2},//run1
				      {-2,-2,-2},//run2
				      {-2,-2,-2},//run3
				      {-7,-2,-2}, //runs 1+2+3 normal binning
				      {-2,-2,-2}};//runs 1+2+3 xsec binning

  double muon_bnb_ymin[num_runs][num_var] = {{2.7,2.7,2.7},//jan
				      {2.7,2.7,2.7},//run1
				      {2.7,2.7,2.7},//run2
				      {2.7,2.7,2.7},//run3
				      {1.5,2.7,2.7}, //runs 1+2+3 normal binning
				      {1.5,2.7,2.7}};//runs 1+2+3 xsec binning

 //Leading Proton
  TH1D* h_leading_overlay[num_var][num_channels];
  std::vector<TH1D*> h_leading_overlay_vec;
  TH1D* h_leading_overlay0[num_var][3];
  TH1D* h_leading_stat[num_var]; //statistical uncertainty
  TH1D* h_leading_systematic[num_var]; //systematic uncertainty
  TH1D* h_leading_bnb[num_var][2];
  TH1D* h_leading_ext[num_var][3];
  TH1D* h_leading_dirt[num_var][3];
  THStack* h_leading[num_var];
  TCanvas* canv_leading[num_var];
  TLegend* legend_leading[num_var];
  TPad* pad_leading[num_var];
  TPad* pad0_leading[num_var];
  std::vector<bool> flip_lead = {false,true,false};
  int leading_ylim[num_runs][num_var] = {{35,160,100},//jan
					 {130,400,275},//run1
					 {200,600,500},//run2
					 {200,600,500},//run3
					 {400,800,800},//runs 1+2+3 normal binning
					 {2500,1800,1500}};//runs 1+2+3 xsec binning
  
  double leading_ymin[num_runs][num_var] = {{-2,-2,-2},//jan
					 {-2,-2,-2},//run1
					 {-2,-2,-2},//run2
					 {-2,-2,-2},//run3
					 {-7,-2,-2},//runs 1+2+3 normal binning
					 {-2,-2,-2}};//runs 1+2+3 xsec binning
 
  //Recoil Proton
  TH1D* h_recoil_overlay[num_var][num_channels];
  std::vector<TH1D*> h_recoil_overlay_vec;
  TH1D* h_recoil_overlay0[num_var][3];
  TH1D* h_recoil_stat[num_var]; //statistical uncertainty
  TH1D* h_recoil_systematic[num_var]; //systematic uncertainty
  TH1D* h_recoil_bnb[num_var][2];
  TH1D* h_recoil_ext[num_var][3];
  TH1D* h_recoil_dirt[num_var][3];
  THStack* h_recoil[num_var];
  TCanvas* canv_recoil[num_var];
  TLegend* legend_recoil[num_var];
  TPad* pad_recoil[num_var];
  TPad* pad0_recoil[num_var];
  std::vector<bool> flip_recoil = {false,true,false};
  int recoil_ylim[num_runs][num_var] = {{60,80,150},//jan
					{125,150,400},//run1
					{200,300,575},//run2
					{200,300,575},//run3
					{400,400,850}, //runs 1+2+3 normal binning
					{2200,1300,1300}};//runs 1+2+3 xsec binning
  
  int recoil_ymin[num_runs][num_var] = {{-2,-2,-2},//jan
					{-2,-2,-2},//run1
					{-2,-2,-2},//run2
					{-2,-2,-2}, //run3
					{-2,-2,-2},//runs 1+2+3 normal binning
					{-2,-2,-2},};//runs 1+2+3 xsec binning

  
  //Random Physics Variables
 ///////////////////////////
  static const int num_phys = 2; //NOTE: These values are the same for Raquel's plots. If this changes in the future, be sure to chage this!
  const char* physics[num_phys] = {"_opening_angle_protons_lab","_opening_angle_mu_both"};//,"_PT_squared","_mom_struck_nuc","_tot_pz","_tot_E","_tot_E_minus_beam",};
  const char* physics_titles[num_phys] = {"cos(#gamma_{p_{L},p_{R}})","cos(#gamma_{#mu,p_{sum}})"};//,"P^{T}_{Miss}^{2} (GeV^{2}/c^{2})","Momentum of Struck Nucleon (GeV)","Total P_{z} of the Two Protons", "Total Kinetic Energy of System (GeV/c)","Total Kinetic Energy - Beam Energy (MeV/c)"};

  /* int phys_ylim[num_runs][num_phys] = {{60,60,60,80,50,150,50,250},//jan
				       {150,300,60,80,100,450,125,1000},//run1
				       {250,500,250,80,200,800,200,1500},//run2
				       {250,800,250,250,250,800,200,1500},//run3
				       {400,1200,250,250,250,800,500,3500}, //runs 1+2+3 normal binning
				       {800,1500,250,250,250,800,2800,3500}}; //runs 1+2+3 xsec binning
  
  int phys_ymin[num_runs][num_phys] =  {{-2,-2,-2,-2,-2,-2,-2},//jan
					{-2,-2,-2,-2,-2,-2,-2},//run1
					{-2,-2,-2,-2,-2,-2,-2},//run2
					{-2,-2,-2,-2,-2,-2,-2},//run3
					{-2,-10,-2,-2,-2,-2,-2},//runs 1+2+3 normal binning
					{-2,-2,-2,-2,-2,-2,-2}};//runs 1+2+3 xsec binning

  double phys_bnb_ymin[num_runs][num_phys] =  {{2.7,2.7,2.7,2.7,2.7,2.7,2.7},//jan
					{2.7,2.7,2.7,2.7,2.7,2.7,2.7},//run1
					{2.7,2.7,2.7,2.7,2.7,2.7,2.7},//run2
					{2.7,2.7,2.7,2.7,2.7,2.7,2.7},//run3
					{2.7,2.7,2.7,2.7,2.7,2.7,2.7},//runs 1+2+3 normal binning
					{2.7,2.7,2.7,2.7,2.7,2.7,2.7}};//runs 1+2+3 xsec binning    

  */

  int phys_ylim[num_runs][num_phys] = {{60,60},//jan
				       {150,300},//run1
				       {250,500},//run2
				       {250,800},//run3
				       {400,1200}, //runs 1+2+3 normal binning
				       {1700,1700}}; //runs 1+2+3 xsec binning
  
  int phys_ymin[num_runs][num_phys] =  {{-2,-2},//jan
					{-2,-2},//run1
					{-2,-2},//run2
					{-2,-2},//run3
					{-2,-10},//runs 1+2+3 normal binning
					{-2,-2}};//runs 1+2+3 xsec binning

  double phys_bnb_ymin[num_runs][num_phys] =  {{2.7,2.7},//jan
					{2.7,2.7},//run1
					{2.7,2.7},//run2
					{2.7,2.7},//run3
					{2.7,2.7},//runs 1+2+3 normal binning
					{2.7,2.7}};//runs 1+2+3 xsec binning    
  
  TH1D* h_phys_overlay[num_phys][num_channels];
  std::vector<TH1D*> h_phys_overlay_vec;
  TH1D* h_phys_overlay0[num_phys][3];
  TH1D* h_phys_stat[num_phys]; //statistical uncertainty
  TH1D* h_phys_systematic[num_phys]; //systematic uncertainty
  TH1D* h_phys_bnb[num_phys][2];
  TH1D* h_phys_ext[num_phys][3];
  TH1D* h_phys_dirt[num_phys][3];
  THStack* h_phys[num_phys];
  TCanvas* canv_phys[num_phys];
  TLegend* legend_phys[num_phys];
  TPad* pad_phys[num_phys];
  TPad* pad0_phys[num_phys];

  //STVs
  //////////////////////////
  static const int num_stv =4;
  const char* stv[num_stv] = {"_delta_PT","_delta_phiT","_delta_alphaT","_pn"};
  const char* stv_titles[num_stv] = {"#delta P_{T} (GeV/c)","#delta #phi_{T} (Deg.)","#delta #alpha_{T} (Deg.)","p_{n} (GeV/c)"};
  int stv_ylim[num_runs][num_stv] = {{100,150,120,120},//jan
				     {300,600,425,150},//run1
				     {500,900,900,200},//run2
				     {500,900,900,200},//run3
				     {800,1200,1000,300}, //runs 1+2+3 normal binning
				     {1600,3000,2000,1500}};//runs 1+2+3 xsec binning
  int stv_ymin[num_runs][num_stv] = {{-2,-2,-2,-2},//jan
				     {-2,-2,-2,-2},//run1
				     {-2,-2,-2,-2},//run2
				     {-2,-2,-2,-2},//run3
				     {-12,-30,-20,-2},//runs 1+2+3 normal binning
				     {-2,-2,-2,-2}};//runs 1+2+3 xsec binning

  double stv_bnb_ymin[num_runs][num_stv] = {{2.7,2.7,2.7,2.7},//jan
					    {2.7,2.7,2.7,2.7},//run1
					    {2.7,2.7,2.7,2.7},//run2
					    {2.7,2.7,2.7,2.7},//run3
					    {1.7,1.5,1.5,2.7},//runs 1+2+3 normal binning
					    {2.7,2.7,2.7,2.7}};//runs 1+2+3 xsec binning 
  
  TH1D* h_stv_overlay[num_stv][num_channels];
  std::vector<TH1D*> h_stv_overlay_vec;
  TH1D* h_stv_overlay0[num_stv][3];
  TH1D* h_stv_stat[num_stv]; //statistical uncertainty
  TH1D* h_stv_systematic[num_stv]; //systematic uncertainty
  TH1D* h_stv_bnb[num_stv][2];
  TH1D* h_stv_ext[num_stv][3];
  TH1D* h_stv_dirt[num_stv][3];
  THStack* h_stv[num_stv];
  TCanvas* canv_stv[num_stv];
  TLegend* legend_stv[num_stv];
  TPad* pad_stv[num_stv];
  TPad* pad0_stv[num_stv];

 //////////////////////////////////////////////////
 //PLOTS WITH RAQUEL'S MOTE-CARLO DEFINITIONS
 //////////////////////////////////////////////////

 //General
 //////////////////////////
  bool plot_total_raquel = false; //sanity check on total overlay
  bool plot_ccnue_raquel = true; //do you want to plot ccnue? 
  int z_raquel=0; //dumb indices                                                                                                                                                                                                        
  int f_raquel = 0; //dumb indices 
  static const int num_channels_raquel = 8;
  //const char* channel_raquel[num_channels_raquel] = {"_total", "_ccMEC", "_ccRES","_ccQE","_nc","_ccDIS","_outfv","_ccNue"};
  //std::vector<const char*> channel_legend_raquel = {"Total Overlay","CCMEC","CCRES","CCQE","NC","CCDIS","Out of FV","CC#nu_{e}"};

  const char* channel_raquel[num_channels_raquel] = {"_total", "_ccMEC", "_ccDIS","_ccRES","_ccNue","_ccQE","_outfv","_nc"};
  std::vector<const char*> channel_legend_raquel = {"Total Overlay","CCMEC","CCDIS","CCRES","CC#nu_{e}","CCQE","Out of FV","NC"};
  
 //Plots of the individual particle quantities
 //Note: we are taking all the variables from my definintions since they are the same for Raquel's plots
  std::vector<TH1D*> h_overlay_raquel_vec; //overlay plots
  TH1D* h_overlay_raquel[num_cuts][num_variables][num_channels];
  TH1D* h_overlay0_raquel[num_cuts][num_variables][3]; //copies of the overlay plots for statistics
  THStack* h_raquel[num_cuts][num_variables];
  TCanvas* canv_raquel[num_cuts][num_variables];
  TLegend* legend_raquel[num_cuts][num_variables];
  TPad* pad_raquel[num_cuts][num_variables];
  TPad* pad0_raquel[num_cuts][num_variables];

  
 //We are just defining new cavases, legends, pads, and overlay histograms 
 ////////////////////////////////////////////

 //Muon
 TH1D* h_muon_overlay_raquel[num_var][num_channels_raquel]; //actual overlay
  std::vector<TH1D*> h_muon_overlay_raquel_vec;
 TH1D* h_muon_overlay0_raquel[num_var][3]; //copies for statistics
 THStack* h_muon_raquel[num_var];
 TCanvas* canv_muon_raquel[num_var];
 TLegend* legend_muon_raquel[num_var];
 TPad* pad_muon_raquel[num_var];
 TPad* pad0_muon_raquel[num_var];

 //Recoil Proton
 TH1D* h_recoil_overlay_raquel[num_var][num_channels_raquel];
  std::vector<TH1D*> h_recoil_overlay_raquel_vec;
 TH1D* h_recoil_overlay0_raquel[num_var][3];
 THStack* h_recoil_raquel[num_var];
 TCanvas* canv_recoil_raquel[num_var];
 TLegend* legend_recoil_raquel[num_var];
 TPad* pad_recoil_raquel[num_var];
 TPad* pad0_recoil_raquel[num_var];

 //Leading Proton
 TH1D* h_leading_overlay_raquel[num_var][num_channels_raquel];
  std::vector<TH1D*> h_leading_overlay_raquel_vec;
 TH1D* h_leading_overlay0_raquel[num_var][3];
 THStack* h_leading_raquel[num_var];
 TCanvas* canv_leading_raquel[num_var];
 TLegend* legend_leading_raquel[num_var];
 TPad* pad_leading_raquel[num_var];
 TPad* pad0_leading_raquel[num_var];

 //Rando plots
 //Note: we are taking all the variables from my definintions since they are the same for Raquel's plots
 //We are just defining new cavases, legends, pads, and overlay histograms 
 ////////////////////////////////////////////////////////
  TH1D* h_phys_overlay_raquel[num_phys][num_channels_raquel];
  std::vector<TH1D*> h_phys_overlay_raquel_vec;
 TH1D* h_phys_overlay0_raquel[num_phys][3];
 THStack* h_phys_raquel[num_phys];
 TCanvas* canv_phys_raquel[num_phys];
 TLegend* legend_phys_raquel[num_phys];
 TPad* pad_phys_raquel[num_phys];
 TPad* pad0_phys_raquel[num_phys];

 //STVs
 //Note: we are taking all the variables from my definintions since they are the same for Raquel's plots
 //We are just defining new cavases, legends, pads, and overlay histograms
 ////////////////////////////   
 TH1D* h_stv_overlay_raquel[num_stv][num_channels_raquel];
  std::vector<TH1D*> h_stv_overlay_raquel_vec;
 TH1D* h_stv_overlay0_raquel[num_stv][3];
 THStack* h_stv_raquel[num_stv];
 TCanvas* canv_stv_raquel[num_stv];
 TLegend* legend_stv_raquel[num_stv];
 TPad* pad_stv_raquel[num_stv];
 TPad* pad0_stv_raquel[num_stv];
  
}

#endif
