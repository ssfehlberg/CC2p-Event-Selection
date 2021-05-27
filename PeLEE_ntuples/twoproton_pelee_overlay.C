#define twoproton_pelee_overlay_cxx
#include "twoproton_pelee_overlay.h"
#include <chrono>
using namespace Constants;
using namespace std::chrono;

void twoproton_pelee_overlay::Loop()
{
  auto start = high_resolution_clock::now(); 

  //Define objects of classes
  ////////////////////////////
  helper_funcs cuts; //helper_funcs.h   

  //Making a new Root File that will contain all the histograms that we will want to plot and files with good RSEs:
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Which_Run();
  TFile *tfile = new TFile(Form("root_files/%s/histograms_pelee_overlay_wgt.root",directory),"RECREATE"); //wgt means applying cenntral value MC value
  ofstream myfile;//File that will contain RSE of good events                                          
  ofstream cc2p; //File that will contain good cc2p events                                                                          
  myfile.open(Form("lists/%s/files_filtered_wgt.list",directory));
  cc2p.open(Form("lists/%s/files_filtered_wgt_cc2p.list",directory));
  myfile<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;
  cc2p<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;

  //Define all the histograms I am going to fill and the mc_wgt                                
  /////////////////////////////////////////////////////////////
  Define_Histograms();
  double mc_wgt; //mc cv weight
  int ohshit = 0;
  int ohshit0 = 0;
  int ohshit1 = 0;
  int ohshit2 = 0;
  int ohshit3 = 0;
  int ohshit_denom = 0;
  int ohshit_num = 0;

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  std::cout<<"Total Number of Entries: "<<nentries<<std::endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;

    std::cout<<"-----------------------------------"<<std::endl;
    std::cout<<"BEGINNING TO PROCESS RUN: " <<run << "  SUBRUN: "<< sub << "  EVENT: " << evt <<std::endl;
    std::cout<<"-----------------------------------"<<std::endl;

    //Checking how many nue's & neutrino slices we have
    /////////////////////////////////
    if(nu_pdg == 12) nue++;
    if(nslice == 0){
      neutrinos_0++;
    }else if(nslice == 1){
      neutrinos_1++;
    }else{
      neutrinos_else++;
    }

    //Definig the MC Weight  & filling the damn mc values
    ///////////////////////////////////////////////////////////                                                  
    if(std::isfinite(weightTune) && weightTune <= 100.) {
      mc_wgt = weightSplineTimesTune;
    } else {
      mc_wgt = 1 * weightSpline;
    }

    int mc_n_threshold_muon = 0;
    int mc_n_threshold_proton = 0;
    int mc_n_threshold_pion0 = 0;
    int mc_n_threshold_pionpm = 0;

    for ( size_t p = 0u; p < mc_pdg->size(); ++p ) {
      int pdg = mc_pdg->at( p );
      float energy = mc_E->at( p );
      if ( pdg == 13) {
	double mom = cuts.real_sqrt( std::pow(energy, 2) - std::pow(MASS_MUON, 2) );
	if(_debug) std::cout<<"Value of the Muon Momentum: "<<mom<<std::endl;
	if ( mom > MUON_MOM_CUT ) {
	  mc_n_threshold_muon++;
	}
      } else if ( pdg == 2212 ) {
	double mom = cuts.real_sqrt( std::pow(energy, 2) - std::pow(MASS_PROTON, 2) );
	if(_debug) std::cout<<"Value of the Proton Momentum: "<<mom<<std::endl;
	if ( mom > PROTON_MOM_CUT_LOW && mom < PROTON_MOM_CUT_HIGH) {
	  mc_n_threshold_proton++;
	}
      } else if ( pdg == 111 ) {
	double mom = cuts.real_sqrt( std::pow(energy, 2) - std::pow(MASS_PION0, 2) );
	if(_debug) std::cout<<"Value of the Pion0 Momentum: "<<mom<<std::endl;
	if ( mom > PION0_MOM_CUT) {
	  mc_n_threshold_pion0++;
	}
      } else if ( std::abs(pdg) == 211 ) {
	double mom = cuts.real_sqrt( std::pow(energy, 2) - std::pow(MASS_PIONPM, 2) );
	if(_debug) std::cout<<"Value of the PionPM Momentum: "<<mom<<std::endl;
	if ( mom > CHARGED_PI_MOM_CUT ) {
	  mc_n_threshold_pionpm++;
	}
      }
    }
    
    //MC_Definitions(mc_n_threshold_muon,mc_n_threshold_proton,mc_n_threshold_pion0,mc_n_threshold_pionpm);
    //MC_Definitions();
    if(_debug) std::cout<<"Value of POT_WGT: "<<pot_wgt<<std::endl;
    if(_debug) std::cout<<"Value of MC_Wgt: "<<mc_wgt<<std::endl;
    if(_debug) std::cout<<"Value of Muon Mom Cut: "<<MUON_MOM_CUT<<std::endl;
    if(_debug) std::cout<<"Value of Proton Low Cut: "<<PROTON_MOM_CUT_LOW<<std::endl;
    if(_debug) std::cout<<"Value of Proton High Cut: "<<PROTON_MOM_CUT_HIGH<<std::endl;
    if(_debug) std::cout<<"Number of threshold muons: "<<mc_n_threshold_muon<<std::endl;
    if(_debug) std::cout<<"Number of threshold protons: "<<mc_n_threshold_proton<<std::endl;
    if(_debug) std::cout<<"Number of threshold pion0: "<<mc_n_threshold_pion0<<std::endl;
    if(_debug) std::cout<<"Number of threshold pionpm: "<<mc_n_threshold_pionpm<<std::endl;

   //making sure that the mc_pdg and npfps are the same length cause fuck me
    //////////////////////////////////////////////////////////////////////////
    std::vector<int> testVector;
    std::vector<double> mc_mom_vector;
    for(size_t i = 0u; i < mc_pdg->size(); i++){
      testVector.push_back(mc_pdg->at(i));
      TVector3 mom_temp(mc_px->at(i),mc_py->at(i),mc_pz->at(i));
      mc_mom_vector.push_back(mom_temp.Mag());
    }
    for(int i = 0; i < n_pfps; i++){
      if(i >= mc_pdg->size()){
	testVector.resize(i+1);
	mc_mom_vector.resize(i+1);
      }
    }

    //Filling histograms before any selection is made
    ////////////////////////////////////////////////
    cuts.Overlay_In_FV(10,10,10,10,10,10,reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z);  // to fill the fv bool
    Fill_Histograms_Mine(0, pot_wgt*mc_wgt, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0 ,mc_n_threshold_pionpm,cuts.fv); 
    Fill_Histograms_Raquel(0, pot_wgt*mc_wgt,cuts.fv);
    
    //Filling the Denominator of the Efficiency
    std::vector<int> mc_protons_id; //the mc ids of the protons
    std::vector<double> mc_proton_mom; //the mc momentum of the protons
    std::vector<std::pair<double,int>> zipped; //pair of the id and the momentum to help me identify the leading and recoil proton
    int leading_id_denom;
    int recoil_id_denom;

    if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_proton == 2 && mc_n_threshold_muon == 1 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && cuts.fv == true){

      //muon stuff for filling the xsec efficiency
      bool true_contained_start;
      bool true_contained_end;
      TVector3 vmuon;
      for(size_t j=0u; j < mc_pdg->size(); j++){
        int pdg = mc_pdg->at(j);
	TVector3 true_mom_vector(mc_px->at(j),mc_py->at(j),mc_pz->at(j));

	//fill the muon
	if (pdg == 13){
	  vmuon.SetXYZ(mc_px->at(j),mc_py->at(j),mc_pz->at(j));
	  if(mc_pdg->size() == n_pfps){
	    true_contained_start = cuts.In_FV(10,10,10,10,10,10,mc_vx->at(j),mc_vy->at(j),mc_vz->at(j)); //I don't have the right size mc_vx....sooo                                                                                                                  
	    true_contained_end = cuts.In_FV(0,0,0,0,0,0,mc_endx->at(j),mc_endy->at(j),mc_endz->at(j));//i don't have the right size mc_vx......                                                                                                                       
	  }
	}
      
	//Now to identify the protons
	if (pdg == 2212){
	  mc_proton_mom.push_back(true_mom_vector.Mag());
	  mc_protons_id.push_back(j);
	  zipped.push_back(std::make_pair(true_mom_vector.Mag(),j));
	}
      } //end loop over pfps
      
      TVector3 leading;
      TVector3 recoil;

      if(mc_proton_mom.size() != 2){
	ohshit_denom++;
      }

      //casse we actually want
      if (mc_proton_mom.size() == 2){
	std::sort(zipped.begin(), zipped.end(), greater());
	for(int j=0; j < mc_protons_id.size(); j++){
	  mc_proton_mom[j] = zipped[j].first;
	  mc_protons_id[j] = zipped[j].second;
	}
	leading_id_denom = mc_protons_id[0]; //leading proton id
	recoil_id_denom = mc_protons_id[1]; //recoil proton id
	leading.SetXYZ(mc_px->at(leading_id_denom),mc_py->at(leading_id_denom),mc_pz->at(leading_id_denom));
	recoil.SetXYZ(mc_px->at(recoil_id_denom),mc_py->at(recoil_id_denom),mc_pz->at(recoil_id_denom));

	std::cout<<"Leading ID FOR DENOM: "<<leading_id_denom<<std::endl;
	std::cout<<"Recoil ID FOR DENOM: "<<recoil_id_denom<<std::endl;
	std::cout<<"Leading Momentum in DENOM: "<<leading.Mag()<<std::endl;
	std::cout<<"Recoil Momentum in DENOM: "<<recoil.Mag()<<std::endl;

	for(size_t j=0u; j < mc_pdg->size(); j++){
	  TVector3 backtracker_mom_vector(mc_px->at(j),mc_py->at(j),mc_pz->at(j)); //mc momentum of particular mc particle
	  bool contained_start = cuts.In_FV(10,10,10,10,10,10,mc_vx->at(j),mc_vy->at(j),mc_vz->at(j)); //is the particle start contained in FV
	  bool contained_end = cuts.In_FV(0,0,0,0,0,0,mc_endx->at(j),mc_endy->at(j),mc_endz->at(j)); //is the particle end contained in the detector
	  Fill_Efficiency_Thresholds(true, mc_pdg->at(j), contained_start, contained_end, leading, recoil, backtracker_mom_vector,mc_wgt*pot_wgt); //filling the efficiency to test the thresholds
	} //end of for loop over the mc particles
      } //end of loop over 2
      Fill_Efficiency_XSec(true,true_contained_start,true_contained_end,vmuon,leading,recoil,1.0);//mc_wgt*pot_wgt); 
    } //end of if same size
    
    mc_protons_id.clear();
    mc_proton_mom.clear();
    zipped.clear();

    //Okay. This Selection requires the following things:
    // 1) The reconstructed neutrino vertex is inside the FV 
    // 2) There are exactly 3 PFP's in the event
    // 3) The 3 PFP's are track like objects i.e they all have a track score > 0.8
    // 4) The 3 PFP's are within 4 cm of the Vertex
    // 5) PID

    //1) Check that the event is in the FV
    //////////////////////////////////////
    if(cuts.In_FV(10,10,10,10,10,10,reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z) == false) continue; //10 cm border except from backend of detector
    fvcntr++;
    Fill_Histograms_Mine(1, pot_wgt*mc_wgt, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0,mc_n_threshold_pionpm,cuts.fv); 
    Fill_Histograms_Raquel(1, pot_wgt*mc_wgt,cuts.fv);

    //2) There are exactly 3 PFP's in the Event 
    ///////////////////////////////////////////////////////
    if(n_pfps != 3) continue;
    threepfps++;
    Fill_Histograms_Mine(2, pot_wgt*mc_wgt, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0,mc_n_threshold_pionpm,cuts.fv);
    Fill_Histograms_Raquel(2, pot_wgt*mc_wgt, cuts.fv);
    for(int i = 0; i < n_pfps; i++){
      int track_pdg = testVector.at(i);
      bool contained_start = cuts.In_FV(10,10,10,10,10,10,trk_start_x_v->at(i),trk_start_y_v->at(i),trk_start_z_v->at(i));
      bool contained_end = cuts.In_FV(10,10,10,10,10,10,trk_end_x_v->at(i),trk_end_y_v->at(i),trk_end_z_v->at(i));
      Fill_Track_Plots(0,i,track_pdg,contained_start,contained_end,pot_wgt*mc_wgt);
    }

    //Require that there are exactly 3 tracks whose vertex distance attachment is less than 4 cm
    /////////////////////////////////////////////////////////////////////////////////////////////
    int y = 0;
    int y1 = 0;
    int muons = 0;
    int protons = 0;
    for(int i = 0; i < n_pfps; i ++){                                                                                       
      float track_score = trk_score_v->at(i);                                                                         
      float track_distance = trk_distance_v->at(i);
      float track_pid = trk_llr_pid_score_v->at(i);
      if(track_score >= TRACK_SCORE_CUT){
	y++; 
      }
      if(track_distance <= TRACK_DIST_CUT){
	y1++;
      }                  
      if(track_pid > PID_CUT){
	muons++;
      }
      if(track_pid < PID_CUT){
	protons++;                                             
      }                                                                                                       
    }                                                                                                     

    //Three pfps with track score above 0.8
    if(y != 3) continue; 
    threetrkcntr++;
    Fill_Histograms_Mine(3, pot_wgt*mc_wgt, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0,mc_n_threshold_pionpm,cuts.fv);
    Fill_Histograms_Raquel(3, pot_wgt*mc_wgt, cuts.fv);
    for(int i = 0; i < n_pfps; i++){
      int track_pdg = testVector.at(i);
      bool contained_start = cuts.In_FV(10,10,10,10,10,10,trk_start_x_v->at(i),trk_start_y_v->at(i),trk_start_z_v->at(i));
      bool contained_end = cuts.In_FV(10,10,10,10,10,10,trk_end_x_v->at(i),trk_end_y_v->at(i),trk_end_z_v->at(i));
      Fill_Track_Plots(1,i,track_pdg,contained_start,contained_end,pot_wgt*mc_wgt);
    }

    //Three Tracks connected to the vertex
    if(y1 != 3) continue;
    threetrk_connected++;
    Fill_Histograms_Mine(4, pot_wgt*mc_wgt, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0,mc_n_threshold_pionpm,cuts.fv);
    Fill_Histograms_Raquel(4, pot_wgt*mc_wgt, cuts.fv);
    for(int i = 0; i < n_pfps; i++){
      int track_pdg = testVector.at(i);
      bool contained_start = cuts.In_FV(10,10,10,10,10,10,trk_start_x_v->at(i),trk_start_y_v->at(i),trk_start_z_v->at(i));
      bool contained_end = cuts.In_FV(10,10,10,10,10,10,trk_end_x_v->at(i),trk_end_y_v->at(i),trk_end_z_v->at(i));	
      Fill_Track_Plots(2,i,track_pdg,contained_start,contained_end,pot_wgt*mc_wgt);
    }
    
    //5) PID: One track with PID > 0.6 and 2 tracks with PID < 0.6
    //////////////////////////////////////////////////////////////
    if(muons != 1 && protons != 2) continue;
    pid++;
    Fill_Histograms_Mine(5, pot_wgt*mc_wgt, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0,mc_n_threshold_pionpm,cuts.fv);
    Fill_Histograms_Raquel(5, pot_wgt*mc_wgt, cuts.fv);

    //Identifying the Partciles
    ///////////////////////////
    int muon_id;
    int leading_proton_id;
    int recoil_proton_id;
    std::vector<int> proton_id_vector;
    for(int i=0; i < trk_pfp_id_v->size(); i++){
      int trk_id = trk_pfp_id_v->at(i);
      double trk_pid = trk_llr_pid_score_v->at(i);	

      if(trk_pid > 1 || trk_pid < -1) continue;
      if(trk_pid > PID_CUT) {
	muon_id = trk_id - 1;
      }

      if(trk_pid < PID_CUT){
	proton_id_vector.push_back(trk_id);
	std::cout<<"trk_id inside of proton PID loop: "<<trk_id<<std::endl;
	float energy = trk_energy_proton_v->at(trk_id - 1);
	float mom = (std::sqrt(std::pow(energy + MASS_PROTON,2) - std::pow(MASS_PROTON,2)));
	std::cout<<"Momentum: "<<mom<<std::endl;
      }
    }

    if(proton_id_vector[0] == 0 || proton_id_vector[1] ==0) continue; 
    uhoh++;

    float mom0 = std::sqrt(std::pow(trk_energy_proton_v->at(proton_id_vector[0]-1) + MASS_PROTON,2) - std::pow(MASS_PROTON,2));
    float mom1 = std::sqrt(std::pow(trk_energy_proton_v->at(proton_id_vector[1]-1) + MASS_PROTON,2) - std::pow(MASS_PROTON,2));

    std::cout<<"Value of mom0: "<<mom0<<std::endl;
    std::cout<<"Value of mom1: "<<mom1<<std::endl;

    if (abs(mom0) > abs(mom1)){
      leading_proton_id = proton_id_vector[0] - 1; //you have to do the -1 cause of course the id's are indexed at one like fucking losers
      recoil_proton_id = proton_id_vector[1] - 1;
    }else{
      leading_proton_id = proton_id_vector[1] - 1;
      recoil_proton_id = proton_id_vector[0] - 1;
    }

    //Finally. Let's define some stuff then fill some variables!
    ////////////////////////////////////////////////////////////    
    
    //Muon
    TVector3 vMuon(1,1,1);
    bool muon_start_contained = cuts.In_FV(10,10,10,10,10,10,trk_sce_start_x_v->at(muon_id),trk_sce_start_y_v->at(muon_id),trk_sce_start_z_v->at(muon_id));
    bool muon_end_contained = cuts.In_FV(0,0,0,0,0,0,trk_sce_end_x_v->at(muon_id),trk_sce_end_y_v->at(muon_id),trk_sce_end_z_v->at(muon_id));
    double EMuon = 0;
    if(muon_start_contained == true && muon_end_contained == true){
      EMuon = std::sqrt(std::pow(trk_range_muon_mom_v->at(muon_id),2)+std::pow(MASS_MUON,2)) - MASS_MUON;
      vMuon.SetMag(trk_range_muon_mom_v->at(muon_id));
      contained++;
    } else if (muon_start_contained == true && muon_end_contained == false){
      EMuon = std::sqrt(std::pow(trk_mcs_muon_mom_v->at(muon_id),2)+std::pow(MASS_MUON,2)) - MASS_MUON;
      vMuon.SetMag(trk_mcs_muon_mom_v->at(muon_id));
      uncontained++;
    }
    vMuon.SetTheta(trk_theta_v->at(muon_id));
    vMuon.SetPhi(trk_phi_v->at(muon_id));
    TLorentzVector muon(vMuon[0],vMuon[1],vMuon[2],EMuon);  

    //Leading Proton
    TVector3 vLead(1,1,1);
    float ELead = trk_energy_proton_v->at(leading_proton_id);
    vLead.SetMag(std::sqrt(std::pow(ELead + MASS_PROTON,2) - std::pow(MASS_PROTON,2)));
    vLead.SetTheta(trk_theta_v->at(leading_proton_id));
    vLead.SetPhi(trk_phi_v->at(leading_proton_id));
    TLorentzVector lead(vLead[0],vLead[1],vLead[2],ELead); 

    //Recoil Proton
    TVector3 vRec(1,1,1);
    float ERec = trk_energy_proton_v->at(recoil_proton_id);
    vRec.SetMag(std::sqrt(std::pow(ERec + MASS_PROTON,2) - std::pow(MASS_PROTON,2)));
    vRec.SetTheta(trk_theta_v->at(recoil_proton_id));
    vRec.SetPhi(trk_phi_v->at(recoil_proton_id));
    TLorentzVector rec(vRec[0],vRec[1],vRec[2],ERec); 

    std::cout<<"Reco Leading Proton ID: "<<leading_proton_id<<std::endl;
    std::cout<<"Reco Recoil Proton ID: "<<recoil_proton_id<<std::endl;
    std::cout<<"Value of Lead Momentum Right Before Filling: "<<vLead.Mag()<<std::endl;
    std::cout<<"Value of Recoil Momentum Right Before Filling: "<<vRec.Mag()<<std::endl;

    Fill_Histograms_Particles(mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0, mc_n_threshold_pionpm, cuts.fv, vMuon, muon, vLead, lead, vRec, rec, 1.0);//mc_wgt*pot_wgt);
    Fill_Histograms_Particles_Raquel(vMuon, muon, vLead, lead, vRec, rec, mc_wgt*pot_wgt, cuts.fv);

    //Make sure to fill the efficiency stuff and the migration matrices
    std::vector<int> reco_protons_id; //the mc ids of the protons                                                                                                                                                                                                            
    std::vector<double> reco_proton_mom; //the mc momentum of the protons                                                                                                                                                                                                    
    std::vector<std::pair<double,int>> zipped_reco; //pair of the id and the momentum to help me identify the leading and recoil proton     
    int leading_id_num;
    int recoil_id_num;

    if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_proton == 2 && mc_n_threshold_muon == 1 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && cuts.fv == true){
      //ohshit++;

      //muon stuff only for filling the matrices and the xsec efficiency
      bool true_contained_start;
      bool true_contained_end;
      TVector3 vmuon;
      
      for(int j=0; j < n_pfps; j++){
	int pdg = backtracked_pdg->at(j);
	int tid = backtracked_tid->at(j);
	std::cout<<"Value of i in NUM: "<<j<<std::endl;
	std::cout<<"Value of PDG in NUM: "<<pdg<<std::endl;
	std::cout<<"Value of TID in NUM: "<<tid<<std::endl;
	TVector3 backtracker_mom_vector(backtracked_px->at(j),backtracked_py->at(j),backtracked_pz->at(j));
	std::cout<<"value of momentum before any particles are id'd: "<<backtracker_mom_vector.Mag()<<std::endl;

	//muon
	if(std::abs(pdg) == 13){
	  if(backtracker_mom_vector.Mag() > MUON_MOM_CUT) {
	    vmuon.SetXYZ(backtracked_px->at(j),backtracked_py->at(j),backtracked_pz->at(j));
	    if(mc_pdg->size() == n_pfps){                                                                                                                                                                                                                   
	      true_contained_start = cuts.In_FV(10,10,10,10,10,10,mc_vx->at(j),mc_vy->at(j),mc_vz->at(j)); //I don't have the right size mc_vx....sooo                                                                                         
	      true_contained_end = cuts.In_FV(0,0,0,0,0,0,mc_endx->at(j),mc_endy->at(j),mc_endz->at(j));//i don't have the right size mc_vx......                                                                                              
	    }
	  }
	} //end of muon loop

	//proton loop used to find the leading and recoil ids
	if(std::abs(pdg) == 2212){
	  if(backtracker_mom_vector.Mag() > PROTON_MOM_CUT_LOW && backtracker_mom_vector.Mag() < PROTON_MOM_CUT_HIGH){	  
	    reco_proton_mom.push_back(backtracker_mom_vector.Mag());
	    std::cout<<"Value of Backtracker Momentum in side the Loop: "<<backtracker_mom_vector.Mag()<<std::endl;

	    reco_protons_id.push_back(j);
	    zipped_reco.push_back(std::make_pair(backtracker_mom_vector.Mag(),j));   
	  } //end of mom threshold
	} //end of proton loop
      } //end of loop over npfps

      TVector3 leading;
      TVector3 recoil;

      std::cout<<"Size of the proton momentum vector in the Num: "<<reco_proton_mom.size()<<std::endl;
      for(int i =0; i < reco_proton_mom.size(); i++){
	std::cout<<"Value of Rec_Proton_Mom Vector at "<<i<<" : "<<reco_proton_mom[i]<<std::endl;
      }

      TVector3 test(-9999.,-9999.,-9999.);
      std::cout<<"Just a Test: "<<test.Mag()<<std::endl;

      if(reco_proton_mom.size() == 0){ //no protons exsist
	ohshit0++;
	leading.SetXYZ(-9999.,-9999.,-9999.);
	recoil.SetXYZ(-9999.,-9999.,-9999.);

      }else if (reco_proton_mom.size() == 1){ //only one proton exists
	ohshit1++;
	leading_id_num = reco_protons_id[0];
	std::cout<<"Leading ID in 1 case :"<<leading_id_num<<std::endl;
	leading.SetXYZ(backtracked_px->at(leading_id_num),backtracked_py->at(leading_id_num),backtracked_pz->at(leading_id_num));
	recoil.SetXYZ(-9999.,-9999.,-9999.);

      }else if(reco_proton_mom.size() > 2){ //more than 2 protons
	ohshit3++;
	  
      }else if (reco_proton_mom.size() >= 2){ //what we want...sort of
	ohshit2++;

	std::sort(zipped_reco.begin(), zipped_reco.end(), greater());
	for(int j=0; j < reco_protons_id.size(); j++){
	  reco_proton_mom[j] = zipped_reco[j].first;
	  reco_protons_id[j] = zipped_reco[j].second;
	}

	leading_id_num = reco_protons_id[0];
	recoil_id_num = reco_protons_id[1];
	leading.SetXYZ(backtracked_px->at(leading_id_num),backtracked_py->at(leading_id_num),backtracked_pz->at(leading_id_num));
	recoil.SetXYZ(backtracked_px->at(recoil_id_num),backtracked_py->at(recoil_id_num),backtracked_pz->at(recoil_id_num));

	std::cout<<"Leading ID in 2>= case :"<<leading_id_num<<std::endl;
	std::cout<<"Recoil ID in 2>= case :"<<recoil_id_num<<std::endl;
	std::cout<<"Value of Leading in 2>= case: "<<leading.Mag()<<std::endl;
	std::cout<<"Value of Recoil in 2>= case: "<<recoil.Mag()<<std::endl;

	for(int j=0; j < n_pfps; j++){
	  TVector3 backtracker_mom_vector(backtracked_px->at(j),backtracked_py->at(j),backtracked_pz->at(j));
	  bool contained_start = cuts.In_FV(10,10,10,10,10,10,mc_vx->at(j),mc_vy->at(j),mc_vz->at(j));
          bool contained_end = cuts.In_FV(0,0,0,0,0,0,mc_endx->at(j),mc_endy->at(j),mc_endz->at(j));
	  Fill_Efficiency_Thresholds(false, backtracked_pdg->at(j), contained_start, contained_end, leading, recoil, backtracker_mom_vector, mc_wgt*pot_wgt);//testing the weight garbage mc_wgt);
	} //end of for
      }//end of if mc_proton_mom_size > 2

      Fill_Matrices(vMuon,vmuon,vLead,leading,vRec,recoil,muon_start_contained,true_contained_start,muon_end_contained,true_contained_end, 1.0);//mc_wgt*pot_wgt);//testing the weight garbage mc_wgt*pot_wgt);                                                                
      Fill_Efficiency_XSec(false,true_contained_start,true_contained_end,vmuon,leading,recoil,1.0);//mc_wgt*pot_wgt);                                       

    }//end of if specific event
    
    //Make sure to clean up before you finish

    if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_proton == 2 && mc_n_threshold_muon == 1 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && cuts.fv == true){
      /*      if(leading_proton_id != leading_id_num || recoil_proton_id != recoil_id_num){
	ohshit++;
	std::cout<<"RIGHT BEFORE CLEANUP"<<std::endl;
	std::cout<<"Leading Denom: "<<leading_id_denom<<std::endl;
	std::cout<<"Leading Reco: "<<leading_proton_id<<std::endl;
	std::cout<<"Leading Num: "<<leading_id_num<<std::endl;
	
	std::cout<<"Recoil Denom: "<<recoil_id_denom<<std::endl;
	std::cout<<"Recoil Reco: "<<recoil_proton_id<<std::endl;
	std::cout<<"Recoil Num: "<<recoil_id_num<<std::endl;
	std::cout<<"YOU ARE NOW CLEARED TO FUCK OFF"<<std::endl;
      }
      */
      if(leading_proton_id != leading_id_num || recoil_proton_id != recoil_id_num) continue;
      
    }

    reco_protons_id.clear();
    reco_proton_mom.clear();
    zipped_reco.clear();
    proton_id_vector.clear();
    testVector.clear();
    events_remaining++;

  } //end of Loop over events

  //Before we finish, we need to make the efficiency and purity plots:
  ///////////////////////////////////////////////////////////////////
  std::vector<int> cut_values = {static_cast<int>(nentries),fvcntr,threepfps,threetrkcntr, threetrk_connected, pid};
  for(int i = 0; i < cut_values.size(); i++){
    double eff = double(cc2p0pi[i]) / double(cc2p0pi[0]);
    double purity = double(cc2p0pi[i]) / double(cut_values[i]);
    std::cout<<"Value of Efficinecy After Cut "<<i<<": "<<eff<<std::endl;
    std::cout<<"Value of Purity After Cut "<<i<<": "<<purity<<std::endl;
    eff_graph->SetPoint(i,i+1,eff);
    pur_graph->SetPoint(i,i+1,purity);
  }

  std::cout<<"-----MODULE SUMMARY-----"<<std::endl;
  std::cout << "[ANALYZER] Initial Number of Events: "<<nentries<<" Fraction of Total: "<<float(100.*float(nentries)/float(nentries))<<"%"<<std::endl;
  std::cout << "[ANALYZER] Number of Events with Vertex in FV: "<<fvcntr<<" Fraction of Total: "<<float(100.*float(fvcntr)/float(nentries))<<"%"<<std::endl;
  std::cout << "[ANALYZER] Number of Events with 3 PFPs: "<<threepfps<<" Fraction of Total: "<<float(100.*float(threepfps)/float(nentries))<<"%"<<std::endl;
  std::cout << "[ANALYZER] Number of Events with 3 Tracks: "<<threetrkcntr<<" Fraction of Total: "<<float(100.*float(threetrkcntr)/float(nentries))<<"%"<<std::endl;
  std::cout << "[ANALYZER] Number of Events with 3 Tracks Connected to Vertex: "<<threetrk_connected<<" Fraction of Total: "<<float(100.*float(threetrk_connected)/float(nentries))<<"%"<<std::endl; 
  std::cout << "[ANALYZER] Number of Events with 1 Muon and 2 Protons: "<<pid<<" Fraction of Total: "<<float(100.*float(pid)/float(nentries))<<"%"<<std::endl;
  std::cout << "[ANALYZER] Sanity Check of the Total Number of Events Remaining: "<<events_remaining<<std::endl;
  std::cout <<"-----CLOSING TIME. YOU DON'T HAVE TO GO HOME, BUT YOU CAN'T STAY HERE-----"<<std::endl;
   
  std::cout<<"-----MC GENERATED SUMMARY: RAQUEL-----"<<std::endl;
  std::cout << "[MC_RAQUEL] Initial Number of Events: "<<nentries<<std::endl;
  std::cout << "[MC_RAQUEL] Number of CCQEL Events: "<<qel[0]<<" Fraction of the Total: "<<float(100.*(float(qel[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC_RAQUEL] Number of CCRES Events: "<<res[0]<<" Fraction of the Total: "<<float(100.*(float(res[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC_RAQUEL] Number of CCMEC Events: "<<mec[0]<<" Fraction of the Total: "<<float(100.*(float(mec[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC_RAQUEL] Number of CCCOH Events: "<<coh[0]<<" Fraction of the Total: "<<float(100.*(float(coh[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC_RAQUEL] Number of CCDIS Events: "<<dis[0]<<" Fraction of the Total: "<<float(100.*(float(dis[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC_RAQUEL] Number of CCNue Events: "<<ccnue_raquel[0]<<" Fraction of the Total: "<<float(100.*(float(ccnue_raquel[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC_RAQUEL] Number of OUTFV Events: "<<outfv_raquel[0]<<" Fraction of the Total: "<<float(100.*(float(outfv_raquel[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC_RAQUEL] Number of NC Events: "<<nc_raquel[0]<<" Fraction of the Total: "<<float(100.*(float(nc_raquel[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC_RAQUEL] Number of Else Events: "<<other_raquel[0]<<" Fraction of the Total: "<<float(100.*(float(other_raquel[0])/float(nentries)))<<"%"<<std::endl;
  std::cout <<"-----MR. SAXOBEAT-----"<<std::endl;
  
  std::cout<<"-----MC GENERATED SUMMARY-----"<<std::endl;
  std::cout << "[MC] Initial Number of Events: "<<nentries<<std::endl;
  std::cout << "[MC] Number of CCOpOpi Events: "<<cc0p0pi[0]<<" Fraction of the Total: "<<float(100.*(float(cc0p0pi[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC] Number of CC1p0pi Events: "<<cc1p0pi[0]<<" Fraction of the Total: "<<float(100.*(float(cc1p0pi[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC] Number of CC2p0pi Events: "<<cc2p0pi[0]<<" Fraction of the Total: "<<float(100.*(float(cc2p0pi[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC] Number of CCNp0pi Events: "<<ccNp0pi[0]<<" Fraction of the Total: "<<float(100.*(float(ccNp0pi[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC] Number of CCNp1pi Events: "<<ccNp1pi[0]<<" Fraction of the Total: "<<float(100.*(float(ccNp1pi[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC] Number of CCNpNpi Events: "<<ccNpNpi[0]<<" Fraction of the Total: "<<float(100.*(float(ccNpNpi[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC] Number of CCNue Events: "<<ccnue[0]<<" Fraction of the Total: "<<float(100.*(float(ccnue[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC] Number of OUTFV Events: "<<outfv[0]<<" Fraction of the Total: "<<float(100.*(float(outfv[0])/float(nentries)))<<"%"<<std::endl; 
  std::cout << "[MC] Number of NC Events: "<<nc[0]<<" Fraction of the Total: "<<float(100.*(float(nc[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC] Number of Else Events: "<<other[0]<<" Fraction of the Total: "<<float(100.*(float(other[0])/float(nentries)))<<"%"<<std::endl;
  std::cout <<"-----OPEN UP MY EAGER EYES! CAUSE I'M MR. BRIGHTSIDE-----"<<std::endl;

  std::cout<<"-----MC RECO'D SUMMARY: RAQUEL-----"<<std::endl;
  std::cout << "[MC_RECO_RAQUEL] Initial Number of Events That were Reconstructed: "<<events_remaining<<std::endl;
  std::cout << "[MC_RECO_RAQUEL] Number of CCQEL Events: "<<qel[number-1]<<" Fraction of the Total: "<<float(100.*(float(qel[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO_RAQUEL] Number of CCRES Events: "<<res[number-1]<<" Fraction of the Total: "<<float(100.*(float(res[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO_RAQUEL] Number of CCMEC Events: "<<mec[number-1]<<" Fraction of the Total: "<<float(100.*(float(mec[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO_RAQUEL] Number of CCCOH Events: "<<coh[number-1]<<" Fraction of the Total: "<<float(100.*(float(coh[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO_RAQUEL] Number of CCDIS Events: "<<dis[number-1]<<" Fraction of the Total: "<<float(100.*(float(dis[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO_RAQUEL] Number of CCNue Events: "<<ccnue_raquel[number-1]<<" Fraction of the Total: "<<float(100.*(float(ccnue_raquel[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO_RAQUEL] Number of OUTFV Events: "<<outfv_raquel[number-1]<<" Fraction of the Total: "<<float(100.*(float(outfv_raquel[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO_RAQUEL] Number of NC Events: "<<nc_raquel[number-1]<<" Fraction of the Total: "<<float(100.*(float(nc_raquel[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO_RAQUEL] Number of Else Events: "<<other_raquel[number-1]<<" Fraction of the Total: "<<float(100.*(float(other_raquel[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout <<"-----ONE FOR THE DAGGER, AND ONE FOR THE ONE YOU BELIEVE!!-----"<<std::endl;

  std::cout<<"-----MC RECO'D SUMMARY-----"<<std::endl;
  std::cout << "[MC_RECO] Initial Number of Events That were Reconstructed: "<<events_remaining<<std::endl;
  std::cout << "[MC_RECO] Number of CCOpOpi Events: "<<cc0p0pi[number-1]<<" Fraction of the Total: "<<float(100.*(float(cc0p0pi[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO] Number of CC1p0pi Events: "<<cc1p0pi[number-1]<<" Fraction of the Total: "<<float(100.*(float(cc1p0pi[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO] Number of CC2p0pi Events: "<<cc2p0pi[number-1]<<" Fraction of the Total: "<<float(100.*(float(cc2p0pi[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO] Number of CCNp0pi Events: "<<ccNp0pi[number-1]<<" Fraction of the Total: "<<float(100.*(float(ccNp0pi[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO] Number of CCNp1pi Events: "<<ccNp1pi[number-1]<<" Fraction of the Total: "<<float(100.*(float(ccNp1pi[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO] Number of CCNpNpi Events: "<<ccNpNpi[number-1]<<" Fraction of the Total: "<<float(100.*(float(ccNpNpi[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO] Number of CCNue Events: "<<ccnue[number-1]<<" Fraction of the Total: "<<float(100.*(float(ccnue[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO] Number of OUTFV Events: "<<outfv[number-1]<<" Fraction of the Total: "<<float(100.*(float(outfv[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO] Number of NC Events: "<<nc[number-1]<<" Fraction of the Total: "<<float(100.*(float(nc[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout << "[MC_RECO] Number of Else Events: "<<other[number-1]<<" Fraction of the Total: "<<float(100.*(float(other[number-1])/float(events_remaining)))<<"%"<<std::endl;
  std::cout <<"-----NOTHING REALLY MATTERS. ANYONE CAN SEE. NOTHING REALLY MATTERS. NOTHING REALLY MATTERS TO ME-----"<<std::endl;

  std::cout<<"Neutrinos 0: "<<neutrinos_0<<std::endl;
  std::cout<<"Neutrinos 1: "<<neutrinos_1<<std::endl;
  std::cout<<"Neutrinos Else: "<<neutrinos_else<<std::endl;
  std::cout<<"1mu2p"<<res_count[0]<<std::endl;
  std::cout<<"1mu1p1pi"<<res_count[1]<<std::endl;
  std::cout<<"1muNp"<<res_count[2]<<std::endl;
  //std::cout<<"else"<<res_count[number-1]<<std::endl;
  std::cout<<"Number of Nue: "<<nue<<std::endl;

  std::cout<<"Other Else: "<<other_else<<std::endl;
  std::cout<<"Neutron: "<<neutron<<std::endl;
  std::cout<<"Neutrino: "<<neutrino<<std::endl;
  std::cout<<"Zeros: "<<zeros<<std::endl;

  std::cout<<"Total Protons: "<<total_protons<<std::endl;
  std::cout<<"Contained Protons: "<<contain<<std::endl;
  std::cout<<"Uncontained Protons: "<<uncontain<<std::endl;
  std::cout<<"UhOh: "<<uhoh<<std::endl;

  std::cout<<"cc2p0pi 0: "<<cc2p0pi[0]<<std::endl;
  std::cout<<"cc2p0pi 1: "<<cc2p0pi[1]<<std::endl;
  std::cout<<"cc2p0pi 2: "<<cc2p0pi[2]<<std::endl;

  std::cout<<"Contained: "<<contained<<std::endl;
  std::cout<<"Uncontained: "<<uncontained<<std::endl;
  std::cout<<"Denom Contained: "<<denom_contained<<std::endl;
  std::cout<<"Denom Uncontained: "<<denom_uncontained<<std::endl;
  std::cout<<"Num Contained: "<<num_contained<<std::endl;
  std::cout<<"Num Uncontained: "<<num_uncontained<<std::endl;

  std::cout<<"OHSHI: "<<ohshit<<std::endl;
  std::cout<<"OHshit0: "<<ohshit0<<std::endl;
  std::cout<<"OHshit1: "<<ohshit1<<std::endl;
  std::cout<<"OHshit2: "<<ohshit2<<std::endl;
  std::cout<<"OHshit2: "<<ohshit3<<std::endl;
  std::cout<<"Ohshit_denom: "<<ohshit_denom<<std::endl;
  std::cout<<"OHSHIT_num: "<<ohshit_num<<std::endl;

  //Don't forget to write all of your histograms before you leave!                                                                       
  ///////////////////////////////////////////////////////////////                                                                 
  tfile->cd();
  Write_Histograms(); //function that writes all our histograms                                                      
  eff_graph->Write("eff_graph");
  pur_graph->Write("pur_graph");
  tfile->Close(); //write the root file that contains our histograms                                                    
  myfile.close(); //Write the file that contains the RSE of good events                                                 
  cc2p.close(); //Write the file that contains the RSE of good 1mu2p events
  //ccNp0pi_file.close();

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<minutes>(stop - start); 
  std::cout<<"Program Run Time: "<<duration.count()<<std::endl;
  
} //end of progrm
