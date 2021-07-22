#define mc_efficiency_cxx
#include "mc_efficiency.h"
#include <chrono>
using namespace Constants;
using namespace std::chrono;

void mc_efficiency::Loop()
{
  auto start = high_resolution_clock::now(); 

  //Define objects of classes
  ////////////////////////////
  Cuts cuts; //cuts.h   
  helper_funcs help; //helper_funcs.h

  //Making a new Root File that will contain all the histograms that we will want to plot and files with good RSEs:
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Which_Run();
  TFile *tfile = new TFile(Form("root_files/%s/histograms_mc_eff.root",directory),"RECREATE"); //wgt means applying cenntral value MC value
  ofstream cc2p_reco;//File that will contain RSE of good events                                          
  ofstream cc2p_true; //File that will contain good cc2p events                                                                          
  cc2p_reco.open(Form("lists/%s/files_filtered_eff_reco.list",directory));
  cc2p_true.open(Form("lists/%s/files_filtered_eff_true.list",directory));
  cc2p_reco<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;
  cc2p_true<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;

  //Define all the histograms I am going to fill and the mc_wgt                                
  /////////////////////////////////////////////////////////////
  Define_Histograms();


  double mc_wgt; //mc cv weight
  int denom = 0;
  int god_dammit = 0;
  int ohshit_num = 0;
  int num = 0;
  int total_muon = 0;
  int flip_muon = 0;
  int total_lead = 0;
  int flip_lead = 0;
  int total_recoil = 0;
  int flip_recoil = 0;

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
      if ( std::abs(pdg) == 13) {
	double mom = help.real_sqrt( std::pow(energy, 2) - std::pow(MASS_MUON, 2) );
	if(_debug) std::cout<<"Value of the Muon Momentum: "<<mom<<std::endl;
	if ( mom > MUON_MOM_CUT_LOW && mom < MUON_MOM_CUT_HIGH ) {
	  mc_n_threshold_muon++;
	}
      } else if (std::abs(pdg) == 2212 ) {
	double mom = help.real_sqrt( std::pow(energy, 2) - std::pow(MASS_PROTON, 2) );
	if(_debug) std::cout<<"Value of the Proton Momentum: "<<mom<<std::endl;
	if ( mom > PROTON_MOM_CUT_LOW && mom < PROTON_MOM_CUT_HIGH) {
	  mc_n_threshold_proton++;
	}
      } else if ( pdg == 111 ) {
	double mom = help.real_sqrt( std::pow(energy, 2) - std::pow(MASS_PION0, 2) );
	if(_debug) std::cout<<"Value of the Pion0 Momentum: "<<mom<<std::endl;
	if ( mom > PION0_MOM_CUT) {
	  mc_n_threshold_pion0++;
	}
      } else if (std::abs(pdg) == 211 ) {
	double mom = help.real_sqrt( std::pow(energy, 2) - std::pow(MASS_PIONPM, 2) );
	if(_debug) std::cout<<"Value of the PionPM Momentum: "<<mom<<std::endl;
	if ( mom > CHARGED_PI_MOM_CUT ) {
	  mc_n_threshold_pionpm++;
	}
      }
    }
    
    if(_debug) std::cout<<"Value of POT_WGT: "<<pot_wgt<<std::endl;
    if(_debug) std::cout<<"Value of MC_Wgt: "<<mc_wgt<<std::endl;
    if(_debug) std::cout<<"Value of Muon Low Mom Cut: "<<MUON_MOM_CUT_LOW<<std::endl;
    if(_debug) std::cout<<"Value of Muon High Mom Cut: "<<MUON_MOM_CUT_HIGH<<std::endl;
    if(_debug) std::cout<<"Value of Proton Low Mom Cut: "<<PROTON_MOM_CUT_LOW<<std::endl;
    if(_debug) std::cout<<"Value of Proton High Mom Cut: "<<PROTON_MOM_CUT_HIGH<<std::endl;
    if(_debug) std::cout<<"Number of threshold muons: "<<mc_n_threshold_muon<<std::endl;
    if(_debug) std::cout<<"Number of threshold protons: "<<mc_n_threshold_proton<<std::endl;
    if(_debug) std::cout<<"Number of threshold pion0: "<<mc_n_threshold_pion0<<std::endl;
    if(_debug) std::cout<<"Number of threshold pionpm: "<<mc_n_threshold_pionpm<<std::endl;

    //Okay. This Selection requires the following things
    // 1) The reconstructed neutrino vertex is inside the FV                                                                             
    // 2) There are exactly 3 PFP's in the event
    // 3) The 3 PFP's are track like objects i.e they all have a track score > 0.8
    // 4) The 3 PFP's are within 4 cm of the Vertex
    // 5) PID  
    
    //Is this event a true 1mu2p?
    /////////////////////////////////
    cuts.True_CC2p(true_nu_vtx_sce_x,true_nu_vtx_sce_y,true_nu_vtx_sce_z, ccnc, nu_pdg, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0, mc_n_threshold_pionpm);
    if(cuts.true_cc2p == true){
      cc2p_true << run << " " << sub << " " << evt << " " ;
      cc2p_true << endl;
    }

    //Reconstructed FV Stuff                                                                                                             
    ///////////////////////////////// 
    cuts.In_FV(10,10,10,10,10,10,reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z);

    //now we have to fill some variables real quick
    //////////////////////////////////////////////
    int tracks_w_good_score = 0;
    int tracks_w_good_distance = 0;
    int muons = 0;
    int protons = 0;
    for(int i = 0; i < n_pfps; i ++){
      float track_score = trk_score_v->at(i);
      float track_distance = trk_distance_v->at(i);
      float track_pid = trk_llr_pid_score_v->at(i);
      if(track_score >= TRACK_SCORE_CUT){
        tracks_w_good_score++;
      }
      if(track_distance <= TRACK_DIST_CUT){
        tracks_w_good_distance++;
      }
      if(track_pid >= PID_CUT && track_pid < 1 && track_pid > -1){
        muons++;
      }
      if(track_pid < PID_CUT && track_pid < 1 && track_pid > -1){
        protons++;
      }
    }

    //Are there 3 PFPs, that are tracks, and within 4cm of the vertex. Also 1 mu and 2 protons?
    ///////////////////////////////////////////
    cuts.Event_Selection(n_pfps, tracks_w_good_score, tracks_w_good_distance, muons, protons);

    //This is the hard part. Finding the ID's of the particles and filling the variables
    /////////////////////////////////////////////////////////////////////////////////
    int muon_id;
    bool muon_start_contained;
    bool muon_end_contained;
    TVector3 vMuon(1,1,1);
    int leading_proton_id;
    TVector3 vLead(1,1,1);
    int recoil_proton_id;
    TVector3 vRec(1,1,1);
    std::vector<int> proton_id_vector;
    TVector3 nu_vtx_reco;
    TVector3 nu_vtx_true;

    if(cuts.good_pid == true){

      for(int i=0; i < trk_pfp_id_v->size(); i++){
	int trk_id = trk_pfp_id_v->at(i);
	double trk_pid = trk_llr_pid_score_v->at(i);	
	if(trk_pid >= PID_CUT && trk_pid < 1 && trk_pid > -1) { //muon
	  muon_id = trk_id - 1;
	}
	if(trk_pid < PID_CUT && trk_pid < 1 && trk_pid > -1.0){ //proton
	  proton_id_vector.push_back(trk_id);
	}
      }

      if(proton_id_vector.size() != 2){
	std::cout<<"FUCK IT ALL!"<<std::endl;
      }

      if(proton_id_vector.size() == 2){
	float mom0 = std::sqrt(std::pow(trk_energy_proton_v->at(proton_id_vector[0]-1) + MASS_PROTON,2) - std::pow(MASS_PROTON,2));
	float mom1 = std::sqrt(std::pow(trk_energy_proton_v->at(proton_id_vector[1]-1) + MASS_PROTON,2) - std::pow(MASS_PROTON,2));

	if (abs(mom0) > abs(mom1)){
	  leading_proton_id = proton_id_vector[0] - 1; //you have to do the -1 cause of course the id's are indexed at one like fucking losers
	  recoil_proton_id = proton_id_vector[1] - 1;
	}else{
	  leading_proton_id = proton_id_vector[1] - 1;
	  recoil_proton_id = proton_id_vector[0] - 1;
	}
      }

      //Neutrino Vertex Stuff for Flipping the Protons
      nu_vtx_true.SetXYZ(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z);
      nu_vtx_reco.SetXYZ(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z); 
      if(_debug) std::cout<<"Location of True Neutrino Vertex: ("<<nu_vtx_true[0]<<","<<nu_vtx_true[2]<<","<<nu_vtx_true[2]<<")"<<std::endl;
      if(_debug) std::cout<<"Location of Reco Neutrino Vertex: ("<<nu_vtx_reco[0]<<","<<nu_vtx_reco[2]<<","<<nu_vtx_reco[2]<<")"<<std::endl;

      ///Muon///
    float  muon_track_start_distance_reco = trk_distance_v->at(muon_id); //distance from start to vertex: reconstructed                                                                        
    TVector3 muon_track_end_reco(trk_sce_end_x_v->at(muon_id),trk_sce_end_y_v->at(muon_id),trk_sce_end_z_v->at(muon_id)); //leading track end reco                                            
    muon_track_end_reco -= nu_vtx_reco;
    double muon_track_end_distance_reco = muon_track_end_reco.Mag(); //distance from end to vertex: reconstructed  
      muon_start_contained = cuts.In_FV(10,10,10,10,10,10,trk_sce_start_x_v->at(muon_id),trk_sce_start_y_v->at(muon_id),trk_sce_start_z_v->at(muon_id));
      muon_end_contained = cuts.In_FV(0,0,0,0,0,0,trk_sce_end_x_v->at(muon_id),trk_sce_end_y_v->at(muon_id),trk_sce_end_z_v->at(muon_id));
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
      total_muon++;

    if(_debug) std::cout<<"Muon 4 Vector: ("<<muon[0]<<","<<muon[1]<<","<<muon[2]<<","<<muon[3]<<")"<<std::endl;
    if(_debug) std::cout<<"[Reconstructed] Muon Start Distance: "<<muon_track_start_distance_reco<<" Muon End Distance: "<<muon_track_end_distance_reco<<std::endl;
    if(_debug) std::cout<<"Muon PID Value: "<<trk_llr_pid_score_v->at(muon_id)<<std::endl;

    //flip momentum if this occurs
    if(muon_track_start_distance_reco > muon_track_end_distance_reco){
      vMuon *= (-1.0); //three vector
      muon.SetPxPyPzE(vMuon[0],vMuon[1],vMuon[2],EMuon); //four vector
      flip_muon++;
    }
    if(_debug) std::cout<<"After Flipping: Muon 4 Vector: ("<<muon[0]<<","<<muon[1]<<","<<muon[2]<<","<<muon[3]<<")"<<std::endl;

      ///Leading proton///
      float lead_track_start_distance_reco = trk_distance_v->at(leading_proton_id); //distance from start to vertex: reconstructed 
      TVector3 lead_track_end_reco(trk_sce_end_x_v->at(leading_proton_id),trk_sce_end_y_v->at(leading_proton_id),trk_sce_end_z_v->at(leading_proton_id)); //leading track end reco
      lead_track_end_reco -= nu_vtx_reco;
      double lead_track_end_distance_reco = lead_track_end_reco.Mag(); //distance from end to vertex: reconstructed                                                                                                                                                                 
      double ELead = trk_energy_proton_v->at(leading_proton_id);
      vLead.SetMag(std::sqrt(std::pow(ELead + MASS_PROTON,2) - std::pow(MASS_PROTON,2)));
      vLead.SetTheta(trk_theta_v->at(leading_proton_id));
      vLead.SetPhi(trk_phi_v->at(leading_proton_id));
      TLorentzVector lead(vLead[0],vLead[1],vLead[2],ELead); 
      total_lead++;

      if(_debug) std::cout<<"Leading Proton 4 Vector: ("<<lead[0]<<","<<lead[1]<<","<<lead[2]<<","<<lead[3]<<")"<<std::endl;                                                                                                             
      if(_debug) std::cout<<"[Reconstructed] Leading Start Distance: "<<lead_track_start_distance_reco<<" Leading End Distance: "<<lead_track_end_distance_reco<<std::endl;                                                       
      if(_debug) std::cout<<"Leading PID Value: "<<trk_llr_pid_score_v->at(leading_proton_id)<<std::endl;                                                                                                                                                                                  
      //flip momentum if this occurs                                                                                                                                                                                             
      if(lead_track_start_distance_reco > lead_track_end_distance_reco){
	vLead *= (-1.0); //three vector                                                                                                                                                                                           
	lead.SetPxPyPzE(vLead[0],vLead[1],vLead[2],ELead); //four vector                                                                                                                                                          
	flip_lead++;
      }
      if(_debug) std::cout<<"After Flipping: Leading Proton 4 Vector: ("<<lead[0]<<","<<lead[1]<<","<<lead[2]<<","<<lead[3]<<")"<<std::endl;

      ///Recoil proton///
      float recoil_track_start_distance_reco = trk_distance_v->at(recoil_proton_id); //distance from start to vertex: reconstructed                                                                                                
      TVector3 recoil_track_end_reco(trk_sce_end_x_v->at(recoil_proton_id),trk_sce_end_y_v->at(recoil_proton_id),trk_sce_end_z_v->at(recoil_proton_id));
      recoil_track_end_reco -= nu_vtx_reco;
      double recoil_track_end_distance_reco = recoil_track_end_reco.Mag(); //distance from end to vertex: reconstructed                                                                                                                                                             
      double ERec = trk_energy_proton_v->at(recoil_proton_id);
      vRec.SetMag(std::sqrt(std::pow(ERec + MASS_PROTON,2) - std::pow(MASS_PROTON,2)));
      vRec.SetTheta(trk_theta_v->at(recoil_proton_id));
      vRec.SetPhi(trk_phi_v->at(recoil_proton_id));
      TLorentzVector rec(vRec[0],vRec[1],vRec[2],ERec); 
      total_recoil++;

      if(_debug) std::cout<<"Recoil Proton 4 Vector: ("<<rec[0]<<","<<rec[1]<<","<<rec[2]<<","<<rec[3]<<")"<<std::endl;
      if(_debug) std::cout<<"[Reconstructed] Recoil Start Distance: "<<recoil_track_start_distance_reco<<" Recoil End Distance: "<<recoil_track_end_distance_reco<<std::endl;
      if(_debug) std::cout<<"Recoil PID Value: "<<trk_llr_pid_score_v->at(recoil_proton_id)<<std::endl;

      //flip the momentum if this occurs
      if(recoil_track_start_distance_reco > recoil_track_end_distance_reco){
	vRec *= (-1.0); //three vector                                                                                                                                                                                             
	rec.SetPxPyPzE(vRec[0],vRec[1],vRec[2],ERec); //four vector 
	flip_recoil++;
      }
      if(_debug) std::cout<<"After Flipping: Recoil Proton 4 Vector: ("<<rec[0]<<","<<rec[1]<<","<<rec[2]<<","<<rec[3]<<")"<<std::endl;
    } //end of good pid event

    //We had to add another cut: The muon and protons must have reconstructed momentum within the thresholdss defined                              
    // this is the only way to get the closure test to work                                                                                                                                                                            
    /////////////////////////////////////////////////////                                                                                                                                                                          
    cuts.Reco_Momentum(vMuon,vLead,vRec);

    //Tell me if this is a reco event?
    ////////////////////////////////
    cuts.Reco_Event();
                                                                                              
    //Fill The Counters
    /////////////////////////////////                                                                                               
    cuts.Fill_Counters(true, true_nu_vtx_sce_x,true_nu_vtx_sce_y,true_nu_vtx_sce_z, ccnc, nu_pdg,mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0, mc_n_threshold_pionpm, interaction);

    //Filling the Denominator of the Efficiency
    ////////////////////////////////////////////
    TVector3 vmuon_denom;
    bool true_contained_start;
    bool true_contained_end;
    std::vector<int> mc_protons_id; //the mc ids of the protons
    std::vector<double> mc_proton_mom; //the mc momentum of the protons
    std::vector<std::pair<double,int>> zipped; //pair of the id and the momentum to help me identify the leading and recoil proton
    int leading_id_denom;
    int recoil_id_denom;
    TVector3 vleading_denom;
    TVector3 vrecoil_denom;

    if(cuts.true_cc2p == true){
      denom++;
      
      for(int j=0u; j < mc_pdg->size(); j++){
	int pdg = mc_pdg->at(j);
	TVector3 true_mom_vector(mc_px->at(j),mc_py->at(j),mc_pz->at(j));
	if(std::abs(pdg) == 13){
	  if(true_mom_vector.Mag() > MUON_MOM_CUT_LOW && true_mom_vector.Mag() < MUON_MOM_CUT_HIGH){
	    true_contained_start = cuts.In_FV(10,10,10,10,10,10,mc_vx->at(j),mc_vy->at(j),mc_vz->at(j));                                                                                                         
	    true_contained_end =  cuts.In_FV(0,0,0,0,0,0,mc_endx->at(j),mc_endy->at(j),mc_endz->at(j));   
	    vmuon_denom.SetXYZ(mc_px->at(j),mc_py->at(j),mc_pz->at(j));
	    }
	}
	  
	if(std::abs(pdg) == 2212){
	  if(true_mom_vector.Mag() > PROTON_MOM_CUT_LOW && true_mom_vector.Mag() < PROTON_MOM_CUT_HIGH){
	    mc_proton_mom.push_back(true_mom_vector.Mag());
	    mc_protons_id.push_back(j);
	    zipped.push_back(std::make_pair(true_mom_vector.Mag(),j));
	    }
	}
      } //end loop over pfps   

      if(mc_proton_mom.size() != 2){
	std::cout<<"GOD Dammit"<<std::endl;
	god_dammit++;
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
	vleading_denom.SetXYZ(mc_px->at(leading_id_denom),mc_py->at(leading_id_denom),mc_pz->at(leading_id_denom));
	vrecoil_denom.SetXYZ(mc_px->at(recoil_id_denom),mc_py->at(recoil_id_denom),mc_pz->at(recoil_id_denom));
	
	for(size_t j=0u; j < mc_pdg->size(); j++){
	  TVector3 backtracker_mom_vector(mc_px->at(j),mc_py->at(j),mc_pz->at(j)); //mc momentum of particular mc particle
	  bool contained_start = cuts.In_FV(10,10,10,10,10,10,mc_vx->at(j),mc_vy->at(j),mc_vz->at(j)); //is the particle start contained in FV
	  bool contained_end = cuts.In_FV(10,10,10,10,10,10,mc_endx->at(j),mc_endy->at(j),mc_endz->at(j));
	  Fill_Efficiency_Thresholds(true, mc_pdg->at(j), contained_start, contained_end, vleading_denom, vrecoil_denom, backtracker_mom_vector,mc_wgt*pot_wgt); //filling the efficiency to test the thresholds
	} //end of for loop over the mc particles
      } //end loop over 2
      Fill_Efficiency_XSec(true,true_contained_start,true_contained_end,vmuon_denom,vleading_denom,vrecoil_denom,mc_wgt*pot_wgt);
    } //end of loop over true cc2p event

    //make sure to clear shit 
    mc_protons_id.clear();
    mc_proton_mom.clear();
    zipped.clear();

    //Now to Fill the Numerator of the Efficiency
    /////////////////////////////////////////////////////////////////////
    bool true_contained_start_num;
    bool true_contained_end_num;
    int muon_id_num;
    int leading_id_num;
    int recoil_id_num;
    TVector3 vmuon_num;
    TVector3 vleading_num;
    TVector3 vrecoil_num;
    
    if(cuts.true_cc2p == true && cuts.reconstructed_event == true){
      num++;
      for(int j=0; j < mc_pdg->size(); j++){
        int pdg = mc_pdg->at(j);
        TVector3 true_mom_vector(mc_px->at(j),mc_py->at(j),mc_pz->at(j));
        
	//fill the muon                                                                                                                                                                                                           
	if(std::abs(pdg) == 13){
          if(true_mom_vector.Mag() >  MUON_MOM_CUT_LOW && true_mom_vector.Mag() < MUON_MOM_CUT_HIGH){
            vmuon_num.SetXYZ(mc_px->at(j),mc_py->at(j),mc_pz->at(j));
	    muon_id_num = j;
	    true_contained_start = cuts.In_FV(10,10,10,10,10,10,mc_vx->at(j),mc_vy->at(j),mc_vz->at(j)); //I don't have the right size mc_vx....sooo                                                                              
	    true_contained_end = cuts.In_FV(0,0,0,0,0,0,mc_endx->at(j),mc_endy->at(j),mc_endz->at(j));//i don't have the right size mc_vx......                                                                                   
	  }
        }
        
	//Now to identify the protons                                                                                                                                                                                             
	if(std::abs(pdg) == 2212){
          if(true_mom_vector.Mag() > PROTON_MOM_CUT_LOW && true_mom_vector.Mag() < PROTON_MOM_CUT_HIGH){
            mc_proton_mom.push_back(true_mom_vector.Mag());
            mc_protons_id.push_back(j);
            zipped.push_back(std::make_pair(true_mom_vector.Mag(),j));
	    }
        }
      } //end loop over pfps   

      if(mc_proton_mom.size() != 2){
        ohshit_num++;
      }

      //casse we actually want                                                                                                                                                                                                    
      if (mc_proton_mom.size() == 2){
	std::sort(zipped.begin(), zipped.end(), greater());
        for(int j=0; j < mc_protons_id.size(); j++){
          mc_proton_mom[j] = zipped[j].first;
          mc_protons_id[j] = zipped[j].second;
        }
        leading_id_num = mc_protons_id[0]; //leading proton id                                                                                                                                                                    
	recoil_id_num = mc_protons_id[1]; //recoil proton id                                                                                                                                                                      
	vleading_num.SetXYZ(mc_px->at(leading_id_num),mc_py->at(leading_id_num),mc_pz->at(leading_id_num));
        vrecoil_num.SetXYZ(mc_px->at(recoil_id_num),mc_py->at(recoil_id_num),mc_pz->at(recoil_id_num));

	for(size_t j=0u; j < mc_pdg->size(); j++){
          TVector3 backtracker_mom_vector(mc_px->at(j),mc_py->at(j),mc_pz->at(j)); //mc momentum of particular mc particle                                                                                                        
	  bool contained_start = cuts.In_FV(10,10,10,10,10,10,mc_vx->at(j),mc_vy->at(j),mc_vz->at(j)); //is the particle start contained in FV                                                                                    
	  bool contained_end = cuts.In_FV(0,0,0,0,0,0,mc_endx->at(j),mc_endy->at(j),mc_endz->at(j)); //is the particle end contained in the detector                                                                              
	  Fill_Efficiency_Thresholds(false, mc_pdg->at(j), contained_start, contained_end, vleading_num, vrecoil_num, backtracker_mom_vector,mc_wgt*pot_wgt); //filling the efficiency to test the thresholds                       
	} //end of for loop over the mc particles                                                                                                                                                                                 
      } //end of loop over 2                                                                                                                                                                                                      
      Fill_Matrices(vMuon,vmuon_num,vLead,vleading_num,vRec,vrecoil_num,muon_start_contained,true_contained_start,muon_end_contained,true_contained_end, mc_wgt*pot_wgt);
      Fill_Efficiency_XSec(false,true_contained_start,true_contained_end,vmuon_num,vleading_num,vrecoil_num,mc_wgt*pot_wgt);   
    }//end of good reco event

    //Make sure to clean up before you finish
    mc_protons_id.clear();
    mc_proton_mom.clear();
    zipped.clear();
    proton_id_vector.clear();

  } //end of Loop over events

  //Before we finish, we need to make the efficiency and purity plots:
  ///////////////////////////////////////////////////////////////////
  std::vector<int> cut_values = {static_cast<int>(nentries),fvcntr,threepfps,threetrkcntr, threetrk_connected, pid, reco_muon_mom_cut, reco_lead_mom_cut, reco_recoil_mom_cut};
  for(int i = 0; i < cut_values.size(); i++){
    double eff = double(cuts.cc2p0pi[i]) / double(cuts.cc2p0pi[0]);
    double purity = double(cuts.cc2p0pi[i]) / double(cut_values[i]);
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
  std::cout << "[ANALYZER] Number of Events with Reco. Muon Momentum above 0.1 GeV and below 2.5 GeV: "<<reco_muon_mom_cut<<" Fraction of Total: "<<float(100.*float(reco_muon_mom_cut)/float(nentries))<<"%"<<std::endl;
  std::cout << "[ANALYZER] Number of Events with Reco. Lead Momentum above 0.1 GeV and below 2.5 GeV: "<<reco_lead_mom_cut<<" Fraction of Total: "<<float(100.*float(reco_lead_mom_cut)/float(nentries))<<"%"<<std::endl;
  std::cout << "[ANALYZER] Number of Events with Reco. Recoil Momentum above 0.1 GeV and below 2.5 GeV: "<<reco_recoil_mom_cut<<" Fraction of Total: "<<float(100.*float(reco_recoil_mom_cut)/float(nentries))<<"%"<<std::endl;
  std::cout <<"-----CLOSING TIME. YOU DON'T HAVE TO GO HOME, BUT YOU CAN'T STAY HERE-----"<<std::endl;
   
  std::cout << "[MC] Initial Number of Events: "<<nentries<<std::endl;
  std::cout << "[MC] Number of Generated CC2p0pi Events: "<<cuts.cc2p0pi[0]<<" Fraction of the Total: "<<float(100.*(float(cuts.cc2p0pi[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC] Number of Denominator Entries Events: "<<denom<<" Fraction of the Total: "<<float(100.*(float(denom)/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC] Number of Generated and Reconstructed CC2p0pi Events: "<<cuts.cc2p0pi[6]<<" Fraction of the Total: "<<float(100.*(float(cuts.cc2p0pi[6])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC] Number of Numerator Entries Events: "<<num<<" Fraction of the Total: "<<float(100.*(float(num)/float(nentries)))<<"%"<<std::endl;

  std::cout<<"God dammit: "<<god_dammit<<std::endl;
  std::cout<<"Ohshitnum: "<<ohshit_num<<std::endl;

  std::cout<<"Total Protons: "<<total_protons<<std::endl;
  std::cout<<"Contained Protons: "<<contain<<std::endl;
  std::cout<<"Uncontained Protons: "<<uncontain<<std::endl;

  std::cout<<"Contained: "<<contained<<std::endl;
  std::cout<<"Uncontained: "<<uncontained<<std::endl;
  std::cout<<"Denom Contained: "<<denom_contained<<std::endl;
  std::cout<<"Denom Uncontained: "<<denom_uncontained<<std::endl;
  std::cout<<"Num Contained: "<<num_contained<<std::endl;
  std::cout<<"Num Uncontained: "<<num_uncontained<<std::endl;

  std::cout<<"Total Number of Muons: "<<total_muon<<std::endl;
  std::cout<<"Flip Muon: "<<flip_muon<<std::endl;
  std::cout<<"Total Number of Lead Protons: "<<total_lead<<std::endl;
  std::cout<<"Flip Lead: "<<flip_lead<<std::endl;
  std::cout<<"Total Number of Recoil Protons: "<<total_recoil<<std::endl;
  std::cout<<"Flip Recoil: "<<flip_recoil<<std::endl;

  //Don't forget to write all of your histograms before you leave!                                                                       
  ///////////////////////////////////////////////////////////////                                                                 
  tfile->cd();
  Write_Histograms(); //function that writes all our histograms                                                      
  eff_graph->Write("eff_graph");
  pur_graph->Write("pur_graph");
  tfile->Close(); //write the root file that contains our histograms                                                    
  cc2p_true.close(); //Write the file that contains the RSE of good events                                                 
  cc2p_reco.close(); //Write the file that contains the RSE of good 1mu2p events

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<minutes>(stop - start); 
  std::cout<<"Program Run Time: "<<duration.count()<<std::endl;
  
} //end of progrm