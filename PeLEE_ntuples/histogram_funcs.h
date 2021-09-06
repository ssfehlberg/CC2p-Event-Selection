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
  virtual void Fill_Particles(TVector3 vMuon, TLorentzVector muon, TVector3 vLead, TLorentzVector lead, TVector3 vRec, TLorentzVector rec, double wgt);
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
  static const int track_cut = 3;
  const char* which_track_cut[track_cut] ={"_after_3_pfps","_after_track_score","_after_distance_cut"};
  TH1D* h_track[num_track][track_cut]; //bnb, ext, dirt
  int num_bins_track[num_track] = {30,10,50,50};
  double xlim_low_track[num_track] = {0.0,0.0,0.0,-1.0};
  double xlim_high_track[num_track] = {1.0,10.0,50.0,1.0};

  //All the single particle plots
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
  TH1D* h_recoil[num_var];
  TH1D* h_leading[num_var];
  TH1D* h_opening_angle_protons_lab; //opening angle between the protons in lab frame
  TH1D* h_opening_angle_protons_com; //opening angle between the protons in com frame
  TH1D* h_opening_angle_mu_leading; //opening angle between the muon and lead proton
  TH1D* h_opening_angle_mu_both; //opening angle between both protons and the muon
  TH1D* h_delta_PT; //stvs delta PT
  TH1D* h_delta_alphaT; //stvs delta alphaT
  TH1D* h_delta_phiT; //stvs delta phiT
  TH1D* h_pn; //neutron momentum estimate
  TH1D* h_nu_E; //neutrino energy estimate
  TH1D* h_mom_struck_nuc; //estimate for neutron momentum: NOT GOOD
  TH1D* h_tot_pz; //total pz of the system
  TH1D* h_tot_E; //total energy of the system
  TH1D* h_tot_E_minus_beam; //the total energy minus the beam
  TH1D* h_PT_squared; //the PT squared quantitiy in Raquel's note

  vector<TH1*> h_list; //list of all the 1D histograms

  variables variables; //variables class

  //Other parameters:                                                                                                                                                      
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
void histogram_funcs::Define_Histograms(const char* sample){

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

} //end of define histograms

//Fills the vertex, cosmicIP, and topo score. Applies appropriate weight.
void histogram_funcs::Fill_Histograms(int i, TVector3 reco_nu_vertex,double CosmicIP, double topological_score, double wgt){ // which cut, reco vertex, wgt to apply
  h_vtx_x[i]->Fill(reco_nu_vertex[0],wgt);
  h_vtx_y[i]->Fill(reco_nu_vertex[1],wgt);
  h_vtx_z[i]->Fill(reco_nu_vertex[2],wgt);
  h_topological_score[i]->Fill(topological_score,wgt);
  h_cosmic_impact_parameter[i]->Fill(CosmicIP,wgt);
}

void histogram_funcs::Fill_Particles(TVector3 vMuon, TLorentzVector muon, TVector3 vLead, TLorentzVector lead, TVector3 vRec, TLorentzVector rec, double wgt){

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

  //first index indicates which variable is being filled: mom, KE energy, theta, phi                                  
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

  //All the physics variables
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

  //have to make sure to clear the variables before you leave
  //////////////////////////////////////////////////////////
  variables.momenta.clear();
  variables.detector_angles.clear();
  variables.opening_angles.clear();
  variables.stvs.clear();
  variables.Energies.clear();

  /*
  //first index indicates which variable is being filled: mom, energy, theta, phi                                  
  h_muon[0]->Fill(vMuon.Mag(),wgt);
  h_leading[0]->Fill(vLead.Mag(),wgt);
  h_recoil[0]->Fill(vRec.Mag(),wgt);

  h_muon[1]->Fill(std::sqrt(vMuon.Mag2() + std::pow(MASS_MUON,2))-MASS_MUON,wgt);
  h_leading[1]->Fill(std::sqrt( vLead.Mag2() + std::pow(MASS_PROTON,2))-MASS_PROTON,wgt);
  h_recoil[1]->Fill(std::sqrt( vRec.Mag2() + std::pow(MASS_PROTON,2))-MASS_PROTON,wgt);

  h_muon[2]->Fill(cos(vMuon.Theta()),wgt);
  h_leading[2]->Fill(cos(vLead.Theta()),wgt);
  h_recoil[2]->Fill(cos(vRec.Theta()),wgt);

  h_muon[3]->Fill(vMuon.Phi(),wgt);
  h_leading[3]->Fill(vLead.Phi(),wgt);
  h_recoil[3]->Fill(vRec.Phi(),wgt);

  double EMuon = muon[3];
  double ELead = lead[3];
  double ERec = rec[3];
  double E_tot = (EMuon + MASS_MUON) + ELead + ERec;

  //Beam Stuff
  TVector3 PT_miss(vMuon[0]+vLead[0]+vRec[0],vMuon[1]+vRec[1]+vLead[1],0);
  double Eneutrino = (EMuon+MASS_MUON) + ELead + ERec +((PT_miss.Mag2())/(2.0*35.37)) + 0.0304;
  TVector3 vBeam(0.,0.,Eneutrino); // z-direction is defined along the neutrino direction                            
  TVector3 vq = vBeam - vMuon; // Momentum transfer                                                                  
  TVector3 vmiss = vLead - vq; // Missing momentum        
  double E_tot_minus_beam = (E_tot - Eneutrino) * 1000;
  TVector3 vProton;
  if(add_protons){
    vProton.SetXYZ(vLead[0]+vRec[0],vLead[1]+vRec[1],vLead[2]+vRec[2]);
  }else{
    vProton.SetXYZ(vLead[0],vLead[1],vLead[2]);
  }

  std::cout<<"Value of PT_miss magnitude: "<<PT_miss.Mag()<<std::endl;
  std::cout<<"Value of PT_miss magnitude2: "<<PT_miss.Mag2()<<std::endl;
  std::cout<<"Value of PT_miss magnitude2 divided by 2*35.37: "<<(PT_miss.Mag2())/(2.0*35.37)<<std::endl;
  std::cout<<"Value of Eneutrino: "<<Eneutrino<<std::endl;
  std::cout<<"Value of E_tot_minus_beam: "<<E_tot_minus_beam<<std::endl;

  open_angle = ((vLead[0]*vRec[0])+(vLead[1]*vRec[1])+(vLead[2]*vRec[2]))/(vLead.Mag()*vRec.Mag()); //vLead.Angle(vRec);  //note this is the cos(opening angle) opening angle between the protons                             
  open_angle_mu = ((vLead[0]*vMuon[0])+(vLead[1]*vMuon[1])+(vLead[2]*vMuon[2]))/(vLead.Mag()*vMuon.Mag()); //note this is the cos(opening angle) of the angle between the leading proton and the muon   
  open_angle_mu_proton = ((vProton[0]*vMuon[0])+(vProton[1]*vMuon[1])+(vProton[2]*vMuon[2]))/(vProton.Mag()*vMuon.Mag()); //cos(opening angle) between the total proton momentum vector and the muon
  En = std::sqrt(std::pow(MASS_NEUTRON,2) + vmiss.Mag2()); //energy of struck nucleon   
  delta_pT = (vMuon + vProton).Perp(); //perp takes the magnitude as well. stv delta pt
  delta_phiT = std::acos( (-vMuon.X()*vProton.X() - vMuon.Y()*vProton.Y()) / (vMuon.XYvector().Mod() * vProton.XYvector().Mod())); //stv delta phi t
  TVector2 delta_pT_vec = (vMuon + vProton).XYvector(); 
  delta_alphaT = std::acos( (-vMuon.X()*delta_pT_vec.X()- vMuon.Y()*delta_pT_vec.Y()) / (vMuon.XYvector().Mod() * delta_pT_vec.Mod()) ); //stv delta alpha T
  //TLorentzVector betacm(vmiss[0] + vRec[0] + vBeam[0],vmiss[1] + vRec[1] + vBeam[1], vmiss[2] + vRec[2]+ vBeam[2], En + ERec + Eneutrino); //beta for CM              
  TLorentzVector betacm(vRec[0]+vLead[0]+vMuon[0],vRec[1]+vLead[1]+vMuon[1],vRec[2]+vLead[2]+vMuon[2],ERec+ELead+EMuon+MASS_MUON+2*MASS_PROTON); 
  TVector3 boost = betacm.BoostVector(); //the boost vector                                                          
  lead.Boost(-boost); //boost leading proton                                                                         
  rec.Boost(-boost); //boost recoil proton                                                                           
  muon.Boost(-boost); 

  std::cout<<"Value of the added vectors x : "<<lead[0]+rec[0]+muon[0]<<std::endl;
  std::cout<<"Value of the added vectors y : "<<lead[1]+rec[1]+muon[1]<<std::endl;
  std::cout<<"Value of the added vectors z : "<<lead[2]+rec[2]+muon[2]<<std::endl;

  //Sanity Check
  //TLorentzVector betacm(vRec[0]+vLead[0],vRec[1]+vLead[1],vRec[2]+vLead[2],ERec+ELead);
  //TVector3 boost = betacm.BoostVector(); //the boost vector                                                                                                                                                                                           //lead.Boost(-boost); //boost leading proton                                                                                                                               
  //rec.Boost(-boost); //boost recoil proton                                   
  
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
  h_opening_angle_protons_lab->Fill(open_angle,wgt);
  h_opening_angle_protons_com->Fill(cos_gamma_cm,wgt);
  h_opening_angle_mu_leading->Fill(open_angle_mu,wgt);
  h_opening_angle_mu_both->Fill(open_angle_mu_proton,wgt);
  h_delta_PT->Fill(delta_pT,wgt);
  h_delta_alphaT->Fill(delta_alphaT*180/3.14,wgt);
  h_delta_phiT->Fill(delta_phiT*180/3.14,wgt);
  h_nu_E->Fill(Eneutrino,wgt);
  h_mom_struck_nuc->Fill(p_struck_nuc,wgt);
  h_tot_pz->Fill(pz_tot,wgt);
  h_tot_E->Fill(E_tot,wgt);
  h_tot_E_minus_beam->Fill(E_tot_minus_beam,wgt);
  h_PT_squared->Fill(PT_miss.Mag2(),wgt);
  */

}

//Function to write Histograms.
/////////////////////////////////
void histogram_funcs::Write_Histograms(){ 
  for(int i = 0; i < h_list.size(); i++){
    h_list[i]->Write();
  }
}
