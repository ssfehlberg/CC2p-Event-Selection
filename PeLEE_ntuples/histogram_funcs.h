#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <fstream> 
#include "vector"

class histogram_funcs
{

 public:
  virtual void Define_Histograms(const char* sample); //defines histograms. works for all samples
  virtual void Fill_Histograms(int i, TVector3 reco_nu_vertex, double wgt); //only fills the bnb, ext, & dirt. i indicates cut
  //virtual void Fill_Particles(int mu, int p1, int p2, double wgt); //only fills the bnb, ext, & dirt. Muon ID, Leding ID, Recoil ID
  virtual void Write_Histograms(bool truth); //writes histograms. works for all samples

  //Total Histograms                                                                                                                
  static const int num = 4;
  const char * total[num] = {"npfp","vtx_npfp","ntrack","nshower"};
  TH1D * h_pfp_overlay[num]; //overlay
  TH1D* h_pfp[num]; //bnb, ext, dirt

  //Correlation Histograms                                                                                                          
  static const int num2d = 3;
  const char * total2d[num2d] = {"reco","truth","truth_sce"};
  const char * labelx[num2d] = {"reco x","true x","true x + sce"};
  const char * labely[num2d] = {"reco y","true y","true y + sce"};
  TH2D *h_correlation_overlay[num2d];

  //Now to define all the specific channels and their histograms                                                                   
  //Since I am basically a plot factory now, I am going to try and do this the smart way                                           
  ///////////////////////////////////////////////////////////////////////////////////////                                          
  static const int  number=5; //before and after selection                                                                         
  static const int  number2 = 11; //categories I defined                                                                            
  static const int  number3 = 10; //categories raquel defined
  const char * point[number] ={"_before_selection","_after_fv","_after_pfp_containment","_after_topo","_after_cosmicIP"}; //this defines histograms before and after the selection                        
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
  TH1D* h_vtx_x[number]; //reco x: bnb, ext, dirt                                                                                   
  TH1D* h_vtx_y[number]; //reco y: bnb, ext, dirt                                                                                  
  TH1D* h_vtx_z[number]; //reco z: bnb, ext, dirt   

  TH1D* h_topological_score_overlay[number2]; //overlay
  TH1D* h_topological_score_raquel[number3]; //raquel
  TH1D* h_topological_score; //bnb,ext,dirt
  TH1D* h_cosmic_impact_parameter_overlay[number2]; //overlay
  TH1D* h_cosmic_impact_parameter_raquel[number3]; //raquel
  TH1D* h_cosmic_impact_parameter; //bnb,ext,dirt
  
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

  //Chi2 Plots
  static const int num_plane = 3;
  const char* plane[num_plane] = {"_Plane_0", "_Plane_1","_Plane_2"};
  static const int num_part = 6;
  const char* particle[num_part] = {"_total","_proton","_muon","_pion","_electron","_other"};
  TH1D* h_chi2p_overlay[num_plane][num_part];
  TH1D* h_chi2mu_overlay[num_plane][num_part];
  TH1D* h_dEdx_total_overlay[num_plane][num_part];
  TH1D* h_chi2p[num_plane];
  TH1D* h_chi2mu[num_plane];
  TH1D* h_dEdx_total[num_plane];

  //3D chi2 plots
  static const int num_3D = 2;
  const char* point_3D[num_3D] = {"_before_selection","_after_selection"};
  TH1D* h_chi2p_3D_overlay[num_3D][num_part]; //3D PID
  TH1D* h_chi2mu_3D_overlay[num_3D][num_part]; //3D PID
  TH1D* h_chi2pi_3D_overlay[num_3D][num_part];
  TH1D* h_chi2p_3D[num_3D]; //3D PID                                                                                                                                                                                                                
  TH1D* h_chi2mu_3D[num_3D]; //3D PID                                                                                                                                                                                                             
  TH1D* h_chi2pi_3D[num_3D];

  //All the single particle plots
  static const int num_var = 4;
  const char* var[num_var] = {"_mom","_E","_theta","_phi"};
  int num_bins[num_var] = {50,50,30,10};
  double xlim_low[num_var] = {0.0,0.0,-1.5,-3.15};
  double xlim_high_recoil[num_var] = {0.8,0.35,1.5,3.15};
  double xlim_high_leading[num_var] = {1.5,0.6,1.5,3.15};
  double xlim_high_muon[num_var]={1.2,1,1.5,3.15};
  const char* xlabel[num_var] ={"P [GeV/c]","E [GeV]","cos(#theta)","#phi [Rad]"};
  TH1D* h_muon_overlay[num_var][number2]; //overlay
  TH1D* h_recoil_overlay[num_var][number2];
  TH1D* h_leading_overlay[num_var][number2];
  TH1D* h_opening_angle_protons_overlay[number2];
  TH1D* h_opening_angle_mu_leading_overlay[number2];
  TH1D* h_delta_PT_overlay[number2];
  TH1D* h_delta_alphaT_overlay[number2];
  TH1D* h_delta_phiT_overlay[number2];
  TH1D* h_cos_gamma_cm_overlay[number2];

  TH1D* h_muon[num_var]; //bnb,ext,dirt
  TH1D* h_recoil[num_var];
  TH1D* h_leading[num_var];
  TH1D* h_opening_angle_protons;
  TH1D* h_opening_angle_mu_leading;
  TH1D* h_delta_PT;
  TH1D* h_delta_alphaT;
  TH1D* h_delta_phiT;
  TH1D* h_cos_gamma_cm;

  TH1D* h_muon_raquel[num_var][number3];
  TH1D* h_recoil_raquel[num_var][number3];
  TH1D* h_leading_raquel[num_var][number3];
  TH1D* h_opening_angle_protons_raquel[number3];
  TH1D* h_opening_angle_mu_leading_raquel[number3];
  TH1D* h_delta_PT_raquel[number3];
  TH1D* h_delta_alphaT_raquel[number3];
  TH1D* h_delta_phiT_raquel[number3];
  TH1D* h_cos_gamma_cm_raquel[number3];
   
  vector<TH1*> h_list; //list of all the 1D histograms

}; //end of class

void histogram_funcs::Define_Histograms(const char* sample){
  
  //THIS IS TO DEFINE THE OVERLAY HISTOGRAMS
  //////////////////////////////////////
  if(strncmp(sample,"overlay",7) == 0){

    //Total Histograms                                                                                                       
    for(int i=0; i < num; i++){
      h_pfp_overlay[i] = new TH1D(Form("h_%s_overlay",total[i]),Form("h_%s_overlay",total[i]),10,0,10);
      h_list.push_back(h_pfp_overlay[i]);
    }

    //Correlation Histograms                                                                                                           
    for(int i=0; i < num2d; i++){
      h_correlation_overlay[i] = new TH2D(Form("h_correlation_overlay_%s",total2d[i]),Form(";%s ;%s",labelx[i],labely[i]),40,0,275,40,-125,-125);
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
	h_chi2p_overlay[k][l] = new TH1D(Form("h_chi2p%s%s",plane[k],particle[l]),Form("h_chi2p%s%s",plane[k],particle[l]),50,0,400);
	h_chi2mu_overlay[k][l] = new TH1D(Form("h_chi2mu%s%s",plane[k],particle[l]),Form("h_chi2mu%s%s",plane[k],particle[l]),50,0,120);
	h_dEdx_total_overlay[k][l] = new TH1D(Form("h_dEdx_total%s%s",plane[k],particle[l]),Form("h_dEdx_total%s%s",plane[k],particle[l]),20,0,1000);
	h_list.push_back(h_chi2p_overlay[k][l]);
	h_list.push_back(h_chi2mu_overlay[k][l]);
	h_list.push_back(h_dEdx_total_overlay[k][l]);
      }
    }

    for(int j=0; j < num_3D; j++){ 
      for(int k=0; k < num_part; k++){
	h_chi2p_3D_overlay[j][k] = new TH1D(Form("h_chi2p_3D%s%s",point_3D[j],particle[k]),Form("h_chi2p_3D%s%s",point_3D[j],particle[k]),50,0,350);
	h_chi2mu_3D_overlay[j][k] = new TH1D(Form("h_chi2mu_3D%s%s",point_3D[j],particle[k]),Form("h_chi2mu_3D%s%s",point_3D[j],particle[k]),50,0,120);
	h_chi2pi_3D_overlay[j][k] = new TH1D(Form("h_chi2pi_3D%s%s",point_3D[j],particle[k]),Form("h_chi2pi_3D%s%s",point_3D[j],particle[k]),50,0,120);
	h_list.push_back(h_chi2p_3D_overlay[j][k]);
	h_list.push_back(h_chi2mu_3D_overlay[j][k]);
	h_list.push_back(h_chi2pi_3D_overlay[j][k]);
      }
    }

    //particle specific plots
    for(int j = 0; j < num_var; j++){
      for(int k = 0; k < number2; k++){
	h_muon_overlay[j][k] = new TH1D(Form("h_muon%s%s",var[j],channel[k]),Form(" h_muon%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_muon[j]);
	h_recoil_overlay[j][k] = new TH1D(Form("h_recoil%s%s",var[j],channel[k]),Form("h_recoil%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_recoil[j]);
	h_leading_overlay[j][k] = new TH1D(Form("h_leading%s%s",var[j],channel[k]),Form("h_leading%s%s ;%s; Counts",var[j],channel[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_leading[j]);
	h_list.push_back(h_muon_overlay[j][k]);
	h_list.push_back(h_recoil_overlay[j][k]);
	h_list.push_back(h_leading_overlay[j][k]);
      }
    }

    for(int j = 0; j < num_var; j++){
      for(int k = 0; k < number3; k++){
	h_muon_raquel[j][k] = new TH1D(Form("h_muon_raquel%s%s",var[j],channel2[k]),Form(" h_muon_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_muon[j]);
	h_recoil_raquel[j][k] = new TH1D(Form("h_recoil_raquel%s%s",var[j],channel2[k]),Form("h_recoil_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_recoil[j]);
	h_leading_raquel[j][k] = new TH1D(Form("h_leading_raquel%s%s",var[j],channel2[k]),Form("h_leading_raquel%s%s ;%s; Counts",var[j],channel2[k],xlabel[j]),num_bins[j],xlim_low[j],xlim_high_leading[j]);
	h_list.push_back(h_muon_raquel[j][k]);
	h_list.push_back(h_recoil_raquel[j][k]);
	h_list.push_back(h_leading_raquel[j][k]);
      }
    }

    //more particle specific plots
    for(int i = 0; i < number2; i++){
      h_opening_angle_protons_overlay[i] = new TH1D(Form("h_opening_angle_protons%s",channel[i]),Form("h_opening_angle_protons%s; Opening Angle btwn Two Protons; Counts",channel[i]),30,-1.5,1.5); //50, 0, 1.5        
      h_opening_angle_mu_leading_overlay[i] = new TH1D(Form("h_opening_angle_mu_leading%s",channel[i]),Form("h_opening_angle_mu_leading%s;Opening Angle btwn Muon and Leading Proton; Counts",channel[i]),30,-1.5,1.5);
      h_delta_PT_overlay[i] = new TH1D(Form("h_delta_PT%s",channel[i]),Form("h_deltaPT%s;#delta P_{T} [GeV/c];Counts",channel[i]),10,0,1);
      h_delta_alphaT_overlay[i] = new TH1D(Form("h_delta_alphaT%s",channel[i]),Form("h_delta_alphaT%s; #delta #alpha_{T} [Deg.];Counts",channel[i]),10,0,180); //0,180                 
      h_delta_phiT_overlay[i] = new TH1D(Form("h_delta_phiT%s",channel[i]),Form("h_delta_phiT%s; #delta #phi_{T} [Deg.];Counts",channel[i]),10,0,180); //0,180                   
      h_cos_gamma_cm_overlay[i] = new TH1D(Form("h_cos_gamma_cm%s",channel[i]),Form("h_cos_gamma_cm%s; cos(#gamma_{COM}); Counts",channel[i]),30,-1.5,1.5);
      h_cosmic_impact_parameter_overlay[i] = new TH1D(Form("h_cosmic_impact_parameter%s",channel[i]),Form("h_cosmic_impact_parameter%s; Cosmic Impact Distance (cm); No. Events",channel[i]),100,0,300);
      h_topological_score_overlay[i] = new TH1D(Form("h_topological_score%s",channel[i]),Form("h_topological_score%s; Topological Score; No. Events",channel[i]),100,0.0,1.0);

      h_list.push_back(h_topological_score_overlay[i]);
      h_list.push_back(h_cosmic_impact_parameter_overlay[i]);
      h_list.push_back(h_cos_gamma_cm_overlay[i]);
      h_list.push_back(h_opening_angle_protons_overlay[i]);
      h_list.push_back(h_opening_angle_mu_leading_overlay[i]);
      h_list.push_back(h_delta_PT_overlay[i]);
      h_list.push_back(h_delta_alphaT_overlay[i]);
      h_list.push_back(h_delta_phiT_overlay[i]);
    }

    for(int i = 0; i < number3; i++){
      h_opening_angle_protons_raquel[i] = new TH1D(Form("h_opening_angle_protons_raquel%s",channel2[i]),Form("h_opening_angle_protons_raquel%s; Opening Angle btwn Two Protons; Counts",channel2[i]),30,-1.5,1.5); //50, 0, 1.5        
      h_opening_angle_mu_leading_raquel[i] = new TH1D(Form("h_opening_angle_mu_leading_raquel%s",channel2[i]),Form("h_opening_angle_mu_leading_raquel%s;Opening Angle btwn Muon and Leading Proton; Counts",channel2[i]),30,-1.5,1.5);
      h_delta_PT_raquel[i] = new TH1D(Form("h_delta_PT_raquel%s",channel2[i]),Form("h_deltaPT_raquel%s;#delta P_{T} [GeV/c];Counts",channel2[i]),10,0,1);
      h_delta_alphaT_raquel[i] = new TH1D(Form("h_delta_alphaT_raquel%s",channel2[i]),Form("h_delta_alphaT_raquel%s; #delta #alpha_{T} [Deg.];Counts",channel2[i]),10,0,180); //0,180                 
      h_delta_phiT_raquel[i] = new TH1D(Form("h_delta_phiT_raquel%s",channel2[i]),Form("h_delta_phiT_raquel%s; #delta #phi_{T} [Deg.];Counts",channel2[i]),10,0,180); //0,180                   
      h_cos_gamma_cm_raquel[i] = new TH1D(Form("h_cos_gamma_cm_raquel%s",channel2[i]),Form("h_cos_gamma_cm_raquel%s; cos(#gamma_{COM}); Counts",channel2[i]),30,-1.5,1.5);
      h_cosmic_impact_parameter_raquel[i] = new TH1D(Form("h_cosmic_impact_parameter_raquel%s",channel2[i]),Form("h_cosmic_impact_parameter_raquel%s; Cosmic Impact Distance (cm); No. Events",channel2[i]),100,0,300);
      h_topological_score_raquel[i] = new TH1D(Form("h_topological_score_raquel%s",channel2[i]),Form("h_topological_score_raquel%s; Topological Score; No. Events",channel2[i]),100,0.0,1.0);

      h_list.push_back(h_topological_score_raquel[i]);
      h_list.push_back(h_cosmic_impact_parameter_raquel[i]);
      h_list.push_back(h_cos_gamma_cm_raquel[i]);
      h_list.push_back(h_opening_angle_protons_raquel[i]);
      h_list.push_back(h_opening_angle_mu_leading_raquel[i]);
      h_list.push_back(h_delta_PT_raquel[i]);
      h_list.push_back(h_delta_alphaT_raquel[i]);
      h_list.push_back(h_delta_phiT_raquel[i]);
    }

    //make sure to handle the weights correcly
    for (int i = 0; i < h_list.size(); i++){
      h_list[i]->Sumw2();
    }
    ///for(int i = 0; i < h_list_2D.size(); i++){
    // h_list_2D[i]->Sumw2();
    // }

  //NOW TO DEFINE THE HISTOGRAMS FOR BNB, EXT AND DIRT
  ////////////////////////////////////////////////////
  }else{

    //Total Histograms                                                                                                        
    for(int i=0; i < num; i++){
      h_pfp[i] = new TH1D(Form("h_%s_%s",total[i],sample),Form("h_%s_%s",total[i],sample),10,0,10);
      h_list.push_back(h_pfp[i]);
    }

    //Other histograms
    for(int i=0; i< number; i++){
      h_vtx_x[i]=new TH1D(Form("h_vtx_x%s_%s",point[i],sample),Form("h_vtx_x%s_%s",point[i],sample),50,0,250);
      h_vtx_y[i]=new TH1D(Form("h_vtx_y%s_%s",point[i],sample),Form("h_vtx_y%s_%s",point[i],sample),50,-125,125);
      h_vtx_z[i]=new TH1D(Form("h_vtx_z%s_%s",point[i],sample),Form("h_vtx_z%s_%s",point[i],sample),50,0,1050);
      h_list.push_back(h_vtx_x[i]);
      h_list.push_back(h_vtx_y[i]);
      h_list.push_back(h_vtx_z[i]);
    }

    for(int k=0; k < num_plane; k++){
      h_chi2p[k] = new TH1D(Form("h_chi2p%s_%s",plane[k],sample),Form("h_chi2p%s_%s",plane[k],sample),50,0,400);
      h_chi2mu[k] = new TH1D(Form("h_chi2mu%s_%s",plane[k],sample),Form("h_chi2mu%s_%s",plane[k],sample),50,0,120);
      h_dEdx_total[k] = new TH1D(Form("h_dEdx_total%s_%s",plane[k],sample),Form("h_dEdx_total%s_%s",plane[k],sample),20,0,1000);
      h_list.push_back(h_chi2p[k]);
      h_list.push_back(h_chi2mu[k]);
      h_list.push_back(h_dEdx_total[k]);
    }

    for(int j=0; j < num_3D; j++){ 
      h_chi2p_3D[j] = new TH1D(Form("h_chi2p_3D%s_%s",point_3D[j],sample),Form("h_chi2p_3D%s_%s",point_3D[j],sample),50,0,350);
      h_chi2mu_3D[j] = new TH1D(Form("h_chi2mu_3D%s_%s",point_3D[j],sample),Form("h_chi2mu_3D%s_%s",point_3D[j],sample),50,0,120);
      h_chi2pi_3D[j] = new TH1D(Form("h_chi2pi_3D%s_%s",point_3D[j],sample),Form("h_chi2pi_3D%s_%s",point_3D[j],sample),50,0,120);
      h_list.push_back(h_chi2p_3D[j]);
      h_list.push_back(h_chi2mu_3D[j]);
      h_list.push_back(h_chi2pi_3D[j]);
    }

    for(int j = 0; j < num_var; j++){
      h_muon[j] = new TH1D(Form("h_muon%s_%s",var[j],sample),Form(" h_muon%s_%s ;%s; Counts",var[j],xlabel[j],sample),num_bins[j],xlim_low[j],xlim_high_muon[j]);
      h_recoil[j] = new TH1D(Form("h_recoil%s_%s",var[j],sample),Form("h_recoil%s_%s ;%s; Counts",var[j],xlabel[j],sample),num_bins[j],xlim_low[j],xlim_high_recoil[j]);
      h_leading[j] = new TH1D(Form("h_leading%s_%s",var[j],sample),Form("h_leading%s_%s ;%s; Counts",var[j],xlabel[j],sample),num_bins[j],xlim_low[j],xlim_high_leading[j]);
      h_list.push_back(h_muon[j]);
      h_list.push_back(h_recoil[j]);
      h_list.push_back(h_leading[j]);
    }

    h_opening_angle_protons = new TH1D(Form("h_opening_angle_protons_%s",sample),Form("h_opening_angle_protons_%s; Opening Angle btwn Two Protons; Counts",sample),30,-1.5,1.5); //50, 0, 1.5                                      
    h_opening_angle_mu_leading = new TH1D(Form("h_opening_angle_mu_leading_%s",sample),Form("h_opening_angle_mu_leading_%s;Opening Angle btwn Muon and Leading Proton; Counts",sample),30,-1.5,1.5);
    h_delta_PT = new TH1D(Form("h_delta_PT_%s",sample),Form("h_deltaPT_%s;#delta P_{T} [GeV/c];Counts",sample),10,0,1);
    h_delta_alphaT = new TH1D(Form("h_delta_alphaT_%s",sample),Form("h_delta_alphaT_%s; #delta #alpha_{T} [Deg.];Counts",sample),10,0,180); //0,180 
    h_delta_phiT = new TH1D(Form("h_delta_phiT_%s",sample),Form("h_delta_phiT_%s; #delta #phi_{T} [Deg.];Counts",sample),10,0,180); //0,180     
    h_cos_gamma_cm = new TH1D(Form("h_cos_gamma_cm_%s",sample),Form("h_cos_gamma_cm_%s;cos(#gamma_{COM});Counts",sample),30,-1.5,1.5);
    h_cosmic_impact_parameter = new TH1D(Form("h_cosmic_impact_parameter_%s",sample),Form("h_cosmic_impact_parameter_%s; Cosmic Impact Distance (cm); No. Events",sample),100,0,300);
    h_topological_score = new TH1D(Form("h_topological_score_%s",sample),Form("h_topological_score_%s; Topological Score; No. Events",sample),100,0.0,1.0);

    h_list.push_back(h_cosmic_impact_parameter);
    h_list.push_back(h_cos_gamma_cm);
    h_list.push_back(h_opening_angle_protons);
    h_list.push_back(h_opening_angle_mu_leading);
    h_list.push_back(h_delta_PT);
    h_list.push_back(h_delta_alphaT);
    h_list.push_back(h_delta_phiT);

    //make sure to handle the weights correctly
    for (int i = 0; i < h_list.size(); i++){
      h_list[i]->Sumw2();
    }
    //for(int i = 0; i < h_list_2D.size(); i++){
    // h_list_2D[i]->Sumw2();
    // }

  } //else loop for bnb, dirt, ext
} //end of define histograms

void histogram_funcs::Fill_Histograms(int i, TVector3 reco_nu_vertex, double wgt){ // which cut, reco vertex, wgt to apply
  h_vtx_x[i]->Fill(reco_nu_vertex[0],wgt);
  h_vtx_y[i]->Fill(reco_nu_vertex[1],wgt);
  h_vtx_z[i]->Fill(reco_nu_vertex[2],wgt);
}

/*void twoproton_filtered_bnb::Fill_Particles(int mu, int p1, int p2){
  //first index indicates which variable is being filled: mom, energy, theta, phi                                                        
  h_muon[0]->Fill(reco_mom_muon->at(mu),wgt);
  h_leading[0]->Fill(reco_mom_proton->at(p1),wgt);
  h_recoil[0]->Fill(reco_mom_proton->at(p2),wgt);

  h_muon[1]->Fill(std::sqrt(std::pow(reco_mom_muon->at(mu),2)+std::pow(0.10565837,2))-0.10565837,wgt);
  h_leading[1]->Fill(std::sqrt(std::pow(reco_mom_proton->at(p1),2)+std::pow(0.93827208,2))-0.93827208,wgt);
  h_recoil[1]->Fill(std::sqrt(std::pow(reco_mom_proton->at(p2),2)+std::pow(0.93827208,2))-0.93827208,wgt);

  h_muon[2]->Fill(cos(reco_theta->at(mu)),wgt);
  h_leading[2]->Fill(cos(reco_theta->at(p1)),wgt);
  h_recoil[2]->Fill(cos(reco_theta->at(p2)),wgt);

  h_muon[3]->Fill(reco_phi->at(mu),wgt);
  h_leading[3]->Fill(reco_phi->at(p1),wgt);
  h_recoil[3]->Fill(reco_phi->at(p2),wgt);

  //if((cos(reco_theta->at(p1)) <= -0.72 && cos(reco_theta->at(p1)) >= -0.78) ||(cos(reco_theta->at(p2)) <= -0.72 && cos(reco_theta->at(p2)) >= -0.78)){
    //if(cos(reco_theta->at(p2)) == cos(reco_theta->at(p1))){
      //std::cout<<"RSE: "<<run<<" ,  "<<subrun<<" , "<<event<<std::endl;
      //} 
    //}


  //Specific parameters
  ////////////////////////////////
  TVector3 vMuon(1,1,1);
  double EMuon = std::sqrt(std::pow(reco_mom_muon->at(mu),2)+std::pow(0.10565837,2));
  vMuon.SetMag(reco_mom_muon->at(mu));
  vMuon.SetTheta(reco_theta->at(mu));
  vMuon.SetPhi(reco_phi->at(mu));
  TLorentzVector muon(vMuon[0],vMuon[1],vMuon[2],EMuon); //Recoil proton TLorentzVector    

  TVector3 vLead(1,1,1);
  double ELead = std::sqrt(std::pow(reco_mom_proton->at(p1),2)+std::pow(0.93827208,2));
  vLead.SetMag(reco_mom_proton->at(p1));
  vLead.SetTheta(reco_theta->at(p1));
  vLead.SetPhi(reco_phi->at(p1));
  TLorentzVector lead(vLead[0],vLead[1],vLead[2],ELead);//leading proton TLorentzVector    
		  
  TVector3 vRec(1,1,1);
  double ERec = std::sqrt(std::pow(reco_mom_proton->at(p2),2)+std::pow(0.93827208,2));
  vRec.SetMag(reco_mom_proton->at(p2));
  vRec.SetTheta(reco_theta->at(p2));
  vRec.SetPhi(reco_phi->at(p2));
  TLorentzVector rec(vRec[0],vRec[1],vRec[2],ERec); //Recoil proton TLorentzVector   

  //Beam Stuff
  double PT_miss = vMuon.Perp() + vRec.Perp() + vLead.Perp();
  double Eneutrino = EMuon + (ELead-0.93827208) + (ERec-0.93827208) +(std::pow(PT_miss,2)/(2*353.7)) + 0.304;
  TVector3 vBeam(0.,0.,Eneutrino); // z-direction is defined along the neutrino direction                            
  TVector3 vq = vBeam - vMuon; // Momentum transfer                                                                  
  TVector3 vmiss = vLead - vq; // Missing momentum        
  open_angle = ((vLead[0]*vRec[0])+(vLead[1]*vRec[1])+(vLead[2]*vRec[2]))/(vLead.Mag()*vRec.Mag()); //vLead.Angle(vRec);  //note this is the cos(opening angle)                             
  open_angle_mu = ((vLead[0]*vMuon[0])+(vLead[1]*vMuon[1])+(vLead[2]*vMuon[2]))/(vLead.Mag()*vMuon.Mag()); //note this is the cos(opening angle)   
  En = std::sqrt(std::pow(NEUTRON_MASS,2) + vmiss.Mag2()); //energy of struck nucleon   

  TVector3 vProton;
  if(add_protons){
    vProton.SetXYZ(vLead[0]+vRec[0],vLead[1]+vRec[1],vLead[2]+vRec[2]);
  }else{
    vProton.SetXYZ(vLead[0],vLead[1],vLead[2]);
  }
  delta_pT = (vMuon + vProton).Perp();
  delta_phiT = std::acos( (-vMuon.X()*vProton.X() - vMuon.Y()*vProton.Y()) / (vMuon.XYvector().Mod() * vProton.XYvector().Mod()));
  TVector2 delta_pT_vec = (vMuon + vProton).XYvector();
  delta_alphaT = std::acos( (-vMuon.X()*delta_pT_vec.X()- vMuon.Y()*delta_pT_vec.Y()) / (vMuon.XYvector().Mod() * delta_pT_vec.Mod()) );

  //TLorentzVector betacm(vmiss[0] + vRec[0] + vBeam[0],vmiss[1] + vRec[1] + vBeam[1], vmiss[2] + vRec[2]+ vBeam[2], En + ERec + Eneutrino); //beta for CM              
  TLorentzVector betacm(vRec[0]+vLead[0]+vMuon[0],vRec[1]+vLead[1]+vMuon[1],vRec[2]+vLead[2]+vMuon[2],ERec+ELead+EMuon); 
  TVector3 boost = betacm.BoostVector(); //the boost vector                                                          
  lead.Boost(-boost); //boost leading proton                                                                         
  rec.Boost(-boost); //boost recoil proton                                                                           
  muon.Boost(-boost); 
  
  std::cout<<"Value of the added vectors x : "<<lead[0]+rec[0]+muon[0]<<std::endl;
  std::cout<<"Value of the added vectors y : "<<lead[1]+rec[1]+muon[1]<<std::endl;
  std::cout<<"Value of the added vectors z : "<<lead[2]+rec[2]+muon[2]<<std::endl;
  

  //Sanity Check
  //TLorentzVector betacm(vRec[0]+vLead[0],vRec[1]+vLead[1],vRec[2]+vLead[2],ERec+ELead);
  //TVector3 boost = betacm.BoostVector(); //the boost vector                                                                                                                                                                                                          
  //lead.Boost(-boost); //boost leading proton                                                                                                                               rec.Boost(-boost); //boost recoil proton                                                                                                                                                                                                                          
  //std::cout<<"Value of the added vectors x : "<<lead[0]+rec[0]<<std::endl;
  //std::cout<<"Value of the added vectors y : "<<lead[1]+rec[1]<<std::endl;
  //std::cout<<"Value of the added vectors z : "<<lead[2]+rec[2]<<std::endl;
  

  cos_gamma_cm = cos(lead.Angle(rec.Vect())); //uses Lorentz Vectors                                                                                                  

  //Struck nucleon Momentum:                                                                                                                                             
  TVector3 vector_sum(vMuon[0] + vLead[0] + vRec[0], vMuon[1] + vLead[1] + vRec[1], vMuon[2] + vLead[2] + vRec[2]);
  TVector3 p_struck_nuc_vector(vector_sum[0], vector_sum[1], 0);
  p_struck_nuc = p_struck_nuc_vector.Mag();
  pz_tot = vLead[2] + vRec[2];

  //Some more specific plots
  h_opening_angle_protons->Fill(open_angle,wgt);
  h_opening_angle_mu_leading->Fill(open_angle_mu,wgt);
  h_delta_PT->Fill(delta_pT,wgt);
  h_delta_alphaT->Fill(delta_alphaT*180/3.14,wgt);
  h_delta_phiT->Fill(delta_phiT*180/3.14,wgt);
  h_cos_gamma_cm->Fill(cos_gamma_cm,wgt);
  h_mom_struck_nuc->Fill(p_struck_nuc,wgt);
  h_tot_pz->Fill(pz_tot,wgt);

} //end of Fill Particles
*/
void histogram_funcs::Write_Histograms(bool truth){ 
  if(truth == true){
    for(int i=0; i< num2d; i++){
      h_correlation_overlay[i]->Write();
    }
    for(int i = 0; i < h_list.size(); i++){
      h_list[i]->Write();
    }
  } else if (truth == false){
    for(int i = 0; i < h_list.size(); i++){
      h_list[i]->Write();
    }
  }
}
