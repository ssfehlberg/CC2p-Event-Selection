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
  virtual void Define_Histograms(const char* sample); //defines histograms. works for all samples
  virtual void Fill_Histograms(int i, TVector3 reco_nu_vertex, double CosmicIP, double topological_score,double wgt); //only fills the bnb, ext, & dirt. i indicates cut
  //virtual void Fill_Particles(TVector4 muon, TVector4 p1, TVector4 p2, double wgt); //only fills the bnb, ext, & dirt. Muon ID, Leding ID, Recoil ID
  virtual void Fill_Particles(TVector3 vMuon, TVector3 vLead, TVector3 vRec, double wgt);
  virtual void Write_Histograms(); //writes histograms. works for all samples

  //Total Histograms                                                                                                                
  static const int num = 4;
  const char * total[num] = {"npfp","vtx_npfp","ntrack","nshower"};
  TH1D* h_pfp[num]; //bnb, ext, dirt

  //Now to define all the specific channels and their histograms                                                                   
  //Since I am basically a plot factory now, I am going to try and do this the smart way                                           
  ///////////////////////////////////////////////////////////////////////////////////////                                          
  static const int  number=6; //number cuts                                                                        
  const char * point[number] ={"_before_selection","_after_fv","_after_three_pfps","_after_track_cut","_after_connection_cut","_after_pid"}; //this defines histograms after each cut    
  TH1D* h_vtx_x[number]; //reco x: bnb, ext, dirt                                                                                   
  TH1D* h_vtx_y[number]; //reco y: bnb, ext, dirt                                                                                  
  TH1D* h_vtx_z[number]; //reco z: bnb, ext, dirt   
  TH1D* h_topological_score[number]; //bnb,ext,dirt
  TH1D* h_cosmic_impact_parameter[number];

  //Track related variables
  static const int num_track = 4;
  const char* variable[num_track] = {"_track_score","_track_vertex_distance","_track_length","_track_pid"};
  TH1D* h_track[num_track]; //bnb, ext, dirt
  int num_bins_track[num_track] = {30,10,50,50};
  double xlim_low_track[num_track] = {0.0,0.0,0.0,-1.0};
  double xlim_high_track[num_track] = {1.0,10.0,50.0,1.0};

  //All the single particle plots
  static const int num_var = 4;
  const char* var[num_var] = {"_mom","_E","_theta","_phi"};
  int num_bins[num_var] = {50,50,30,10}; //50 first bin, 10 last bin
  double xlim_low[num_var] = {0,0,-1.5,-3.15}; //0.2 normally first -1.5
  double xlim_high_recoil[num_var] = {0.8,0.35,1.5,3.15};
  double xlim_high_leading[num_var] = {1.2,0.6,1.5,3.15}; //1.5 normally in first, 1.2
  double xlim_high_muon[num_var]={2.5,1,1.5,3.15}; //2.5 first, 1.5 third
  const char* xlabel[num_var] ={"P [GeV/c]","E [GeV]","cos(#theta)","#phi [Rad]"};

  TH1D* h_muon[num_var]; //bnb,ext,dirt
  TH1D* h_recoil[num_var];
  TH1D* h_leading[num_var];
  TH1D* h_opening_angle_protons;
  TH1D* h_opening_angle_mu_leading;
  TH1D* h_delta_PT;
  TH1D* h_delta_alphaT;
  TH1D* h_delta_phiT;
  TH1D* h_cos_gamma_cm;
  TH1D* h_mom_struck_nuc;
  TH1D* h_tot_pz;
  TH1D* h_tot_E;
  TH1D* h_tot_E_minus_beam;
  TH1D* h_E_neutrino;
  TH1D* h_opening_angle_mu_both;
  TH1D* h_PT_squared;

  vector<TH1*> h_list; //list of all the 1D histograms

  variables variables; //variables class

}; //end of class

//DEFINE THE HISTOGRAMS FOR BNB, EXT AND DIRT                                                                                                                                                                                                       
////////////////////////////////////////////////////   
void histogram_funcs::Define_Histograms(const char* sample){

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
    h_cosmic_impact_parameter[i] = new TH1D(Form("h_cosmic_impact_parameter%s_%s",point[i],sample),Form("h_cosmic_impact_parameter%s_%s; Cosmic Impact Distance (cm); No. Events",point[i],sample),20,0,200);
    h_topological_score[i] = new TH1D(Form("h_topological_score%s_%s",point[i],sample),Form("h_topological_score%s_%s; Topological Score; No. Events",point[i],sample),50,0.0,1.0); //30 for Wouter, 50 for steven

    h_list.push_back(h_cosmic_impact_parameter[i]);
    h_list.push_back(h_topological_score[i]);
    h_list.push_back(h_vtx_x[i]);
    h_list.push_back(h_vtx_y[i]);
    h_list.push_back(h_vtx_z[i]);
  }

  for(int j=0; j < num_track; j++){
    h_track[j] = new TH1D(Form("h_track%s",variable[j]),Form("h_track%s",variable[j]),num_bins_track[j],xlim_low_track[j],xlim_high_track[j]);
    h_list.push_back(h_track[j]);
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
  h_delta_PT = new TH1D(Form("h_delta_PT_%s",sample),Form("h_deltaPT_%s;#delta P_{T} [GeV/c];Counts",sample),15,0,1); //normally 10 bins
  h_delta_alphaT = new TH1D(Form("h_delta_alphaT_%s",sample),Form("h_delta_alphaT_%s; #delta #alpha_{T} [Deg.];Counts",sample),10,0,180); //0,180 
  h_delta_phiT = new TH1D(Form("h_delta_phiT_%s",sample),Form("h_delta_phiT_%s; #delta #phi_{T} [Deg.];Counts",sample),10,0,180); //0,180     
  h_cos_gamma_cm = new TH1D(Form("h_cos_gamma_cm_%s",sample),Form("h_cos_gamma_cm_%s;cos(#gamma_{COM});Counts",sample),30,-1.5,1.5);
  h_mom_struck_nuc = new TH1D(Form("h_mom_struck_nuc_%s",sample),Form("h_mom_struck_nuc_%s; P_{Init}; Counts", sample),30, 0, 1);
  h_tot_pz = new TH1D(Form("h_tot_pz_%s",sample),Form("h_tot_pz_%s; P_{Z}^{Total}; Counts",sample), 20, 0, 2);
  h_tot_E = new TH1D(Form("h_tot_E_%s",sample),Form("h_tot_E_%s; Total Energy; Counts;",sample),50,0,2.5);
  h_tot_E_minus_beam = new TH1D(Form("h_tot_E_minus_beam_%s",sample),Form("h_tot_E_minus_beam_%s; Total Energy Remaining (MeV/c); Counts;",sample),100,-100,0);
  h_E_neutrino = new TH1D(Form("h_E_neutrino_%s",sample),Form("h_E_neutrino_%s; Total Energy; Counts;",sample),50,0,2.5);
  h_opening_angle_mu_both = new TH1D(Form("h_opening_angle_mu_both_%s",sample),Form("h_opening_angle_mu_both_%s; Opening Angle btwn Muon and Total Proton Momentum; Counts",sample),30,-1.5,1.5);
  h_PT_squared = new TH1D(Form("h_PT_squared_%s",sample),Form("h_PT_squared_%s; P_{T}^{2}; Counts", sample),50,0,5);

  h_list.push_back(h_PT_squared);
  h_list.push_back(h_opening_angle_mu_both);
  h_list.push_back(h_E_neutrino);
  h_list.push_back(h_tot_E);
  h_list.push_back(h_tot_E_minus_beam);
  h_list.push_back(h_cos_gamma_cm);
  h_list.push_back(h_opening_angle_protons);
  h_list.push_back(h_opening_angle_mu_leading);
  h_list.push_back(h_delta_PT);
  h_list.push_back(h_delta_alphaT);
  h_list.push_back(h_delta_phiT);
  h_list.push_back(h_mom_struck_nuc);
  h_list.push_back(h_tot_pz);

  //make sure to handle the weights correctly
  for (int i = 0; i < h_list.size(); i++){
    h_list[i]->Sumw2();
  }
  //for(int i = 0; i < h_list_2D.size(); i++){
  // h_list_2D[i]->Sumw2();
  // }

} //end of define histograms

//Fills the vertex, cosmicIP, and topo score. Applies appropriate weight.
void histogram_funcs::Fill_Histograms(int i, TVector3 reco_nu_vertex,double CosmicIP, double topological_score, double wgt){ // which cut, reco vertex, wgt to apply
  h_vtx_x[i]->Fill(reco_nu_vertex[0],wgt);
  h_vtx_y[i]->Fill(reco_nu_vertex[1],wgt);
  h_vtx_z[i]->Fill(reco_nu_vertex[2],wgt);
  h_topological_score[i]->Fill(topological_score,wgt);
  h_cosmic_impact_parameter[i]->Fill(CosmicIP,wgt);
}

//Fills more specific plots
void histogram_funcs::Fill_Particles(TVector3 vMuon, TVector3 vLead, TVector3 vRec, double wgt){
  
  // Run the Calculate Variables function inside of variables.h. Returns following:
  // 1) vector: momenta(muon_mom,lead_mom,rec_mom);
  // 2) vector: Energies(KE_muon, TotE_muon, KE_Lead, TotE_Lead, KE_Rec, TotE_Rec);    
  // 3) vector: detector_angles(muon_theta,muon_phi,lead_theta,lead_phi,recoil_theta,recoil_phi);
  // 4) vector: opening_angles(opening_angle_protons_lab,opening_angle_protons_mu_leading,opening_angle_protons_mu_both);  // 5) double: opening_angle_protons_COM 
  // 6) vector: STVS(delta_pT,delta_alphaT,delta_phiT); 
  // 7) double: calculated_nu_E

  variables.Calculate_Variables(vMuon,vLead,vRec,add_protons);

  //first index indicates which variable is being filled: mom, KE energy, theta, phi                                  
  h_muon[0]->Fill(variables.momenta[0],wgt);
  h_leading[0]->Fill(variables.momenta[1],wgt);
  h_recoil[0]->Fill(variables.momenta[2],wgt);

  h_muon[1]->Fill(variables.Energies[0],wgt);
  h_leading[1]->Fill(variables.Energies[2],wgt);
  h_recoil[1]->Fill(variables.Energies[4],wgt);

  h_muon[2]->Fill(variables.detector_angles[0],wgt);
  h_leading[2]->Fill(variables.detector_angles[2],wgt);
  h_recoil[2]->Fill(variables.detector_angles[4],wgt);

  h_muon[3]->Fill(variables.detector_angles[1],wgt);
  h_leading[3]->Fill(variables.detector_angles[3],wgt);
  h_recoil[3]->Fill(variables.detector_angles[5],wgt);
  	
  //Beam Stuff
  double E_tot = (variables.Energies[0] + MASS_MUON) + variables.Energies[2] + variables.Energies[4];
  TVector3 PT_miss(vMuon[0]+vLead[0]+vRec[0],vMuon[1]+vRec[1]+vLead[1],0);
  double Eneutrino = variables.calculated_nu_E;//(EMuon+MASS_MUON) + ELead + ERec +((PT_miss.Mag2())/(2.0*35.37)) + 0.0304;
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

  h_opening_angle_protons->Fill(variables.opening_angles[0],wgt);
  h_cos_gamma_cm->Fill(variables.opening_angle_protons_COM,wgt);
  h_opening_angle_mu_leading->Fill(variables.opening_angles[1],wgt);
  h_opening_angle_mu_both->Fill(variables.opening_angles[2],wgt);
  h_delta_PT->Fill(variables.stvs[0],wgt);
  h_delta_alphaT->Fill(variables.stvs[1],wgt);
  h_delta_phiT->Fill(variables.stvs[2],wgt);
  h_tot_E->Fill(E_tot,wgt);
  h_tot_E_minus_beam->Fill(E_tot_minus_beam,wgt);
  h_E_neutrino->Fill(Eneutrino,wgt);
  h_PT_squared->Fill(PT_miss.Mag2(),wgt);
  h_mom_struck_nuc->Fill(p_struck_nuc,wgt);
  h_tot_pz->Fill(pz_tot,wgt);

}

//Function to write Histograms.
/////////////////////////////////
void histogram_funcs::Write_Histograms(){ 
  for(int i = 0; i < h_list.size(); i++){
    h_list[i]->Write();
  }
}
