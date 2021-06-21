#define efficiency_cxx
#include "efficiency.h"
#include <chrono>
using namespace Constants;
using namespace std::chrono;

void efficiency::Loop()
{
  auto start = high_resolution_clock::now(); 

  //Define objects of classes
  ////////////////////////////
  helper_funcs cuts; //helper_funcs.h   

  //Making a new Root File that will contain all the histograms that we will want to plot and files with good RSEs:
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Which_Run();
  TFile *tfile = new TFile(Form("root_files/%s/histograms_eff.root",directory),"RECREATE"); //wgt means applying cenntral value MC value
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
  int num = 0;
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
	double mom = cuts.real_sqrt( std::pow(energy, 2) - std::pow(MASS_MUON, 2) );
	if(_debug) std::cout<<"Value of the Muon Momentum: "<<mom<<std::endl;
	if ( mom > MUON_MOM_CUT_LOW && mom < MUON_MOM_CUT_HIGH ) {
	  mc_n_threshold_muon++;
	}
      } else if (std::abs(pdg) == 2212 ) {
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
      } else if (std::abs(pdg) == 211 ) {
	double mom = cuts.real_sqrt( std::pow(energy, 2) - std::pow(MASS_PIONPM, 2) );
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
    bool true_cc2p;
    cuts.Overlay_In_FV(10,10,10,10,10,10,true_nu_vtx_sce_x,true_nu_vtx_sce_y,true_nu_vtx_sce_z); //overlay FV requirment
    if(ccnc == 0 && nu_pdg == 14 && mc_n_threshold_proton == 2 && mc_n_threshold_muon == 1 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && cuts.fv == true){
      true_cc2p = true;
      cc2p0pi[0]++;
      cc2p_true << run << " " << sub << " " << evt << " " ;
      cc2p_true << endl;
    }else{
      true_cc2p = false;
    }

    //Reconstructed FV Stuff                                                                                                                                                                                                              ///////////////////////////////// 
    bool reco_in_FV = cuts.In_FV(10,10,10,10,10,10,reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z); 

    //Are there 3 pfps?
    /////////////////////
    bool three_pfps;
    if(n_pfps == 3){
      three_pfps = true;
    }else{
      three_pfps = false;
    }

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

    //Are there 3 tracks with a good track score?
    ///////////////////////////////////////////
    bool three_tracks;
    if(tracks_w_good_score == 3){
      three_tracks = true;
    }else{
      three_tracks = false;
    }
    
    //Are there 3 tracks within 4cm of reco vertex?
    ///////////////////////////////////////////////
    bool three_tracks_4cm_vertex;
    if(tracks_w_good_distance == 3){
      three_tracks_4cm_vertex = true;
    }else{
      three_tracks_4cm_vertex = false;
    }

    //are there exactly 1 muon and 2 protons
    ///////////////////////////////////////
    bool reco_cc2p_event;
    if(muons == 1 && protons == 2){
      reco_cc2p_event = true;
    }else{
      reco_cc2p_event = false;
    }

    //Filling a final bool that just indicates that the event is a good reco event:                                                                                                                                                       ///////////////////////////////////////////////////////////////////////////////                                                                                                                                                  
    bool reconstructed_event;
    if(reco_in_FV == true && three_pfps == true && three_tracks == true && three_tracks_4cm_vertex == true && reco_cc2p_event == true){
      reconstructed_event = true;
    }else{
      reconstructed_event = false;
    }

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

    if(reconstructed_event == true){

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

      ///Muon///
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
      
      //Neutrino Vertex Stuff for Flipping the Protons
      nu_vtx_true.SetXYZ(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z);
      nu_vtx_reco.SetXYZ(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z); 
      if(_debug) std::cout<<"Location of True Neutrino Vertex: ("<<nu_vtx_true[0]<<","<<nu_vtx_true[2]<<","<<nu_vtx_true[2]<<")"<<std::endl;
      if(_debug) std::cout<<"Location of Reco Neutrino Vertex: ("<<nu_vtx_reco[0]<<","<<nu_vtx_reco[2]<<","<<nu_vtx_reco[2]<<")"<<std::endl;

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
    } //end of good reco event

    //We had to add another cut: The muon and protons must have reconstructed momentum within the thresholdss defined                                                                                                                
    // this is the only way to get the closure test to work                                                                                                                                                                         
    /////////////////////////////////////////////////////                                                                                                                                                                            
    bool good_muon_mom;
    if(vMuon.Mag() > MUON_MOM_CUT_LOW && vMuon.Mag() < MUON_MOM_CUT_HIGH){
      good_muon_mom = true;
    }else{
      good_muon_mom = false;
    }

    //Make sure to fill the counters                                                                                                                                                                                                      /////////////////////////////////                                                                                                                                                                                                 

    //fv
    if(reco_in_FV == true){
      fvcntr++;
      if(true_cc2p == true){
	cc2p0pi[1]++;
      }
    }

    //three pfps
    if(reco_in_FV && three_pfps == true){
    threepfps++;
    if(true_cc2p == true){
      cc2p0pi[2]++;
    }
    }

    //three tracks
    if(reco_in_FV==true && three_pfps == true && three_tracks == true){
      threetrkcntr++;
      if(true_cc2p == true){
        cc2p0pi[3]++;
      }
    }

    //three connected tracks
    if(reco_in_FV == true && three_pfps == true && three_tracks == true && three_tracks_4cm_vertex == true){
      threetrk_connected++;
      if(true_cc2p == true){
        cc2p0pi[4]++;
      }
    }
    
    //pid
    if(reco_in_FV == true && three_pfps == true && three_tracks == true && three_tracks_4cm_vertex == true && reco_cc2p_event == true){
      if(true_cc2p == true){
        cc2p0pi[5]++;
      }
      pid++;
    }

    //muon_mom_cut
    if(reconstructed_event == true && good_muon_mom == true){
      reco_muon_mom_cut++;
      cc2p_reco << run << " " << sub << " " << evt << " " ;
      cc2p_reco << endl;
      if(true_cc2p == true){
        cc2p0pi[6]++;
      } 
    }

    //Filling the Denominator of the Efficiency
    ////////////////////////////////////////////
    TVector3 vmuon_denom;
    bool true_contained_start;
    bool true_contained_end;
    std::vector<int> mc_protons_id; //the mc ids of the protons
    std::vector<double> mc_proton_mom; //the mc momentum of the protons
    std::vector<std::pair<double,int>> zipped; //pair of the id and the momentum to help me identify the leading and recoil proton
    TVector3 vleading_denom;
    TVector3 vrecoil_denom;

    if(true_cc2p == true){
      denom++;
      
      //if it is a reco event, then we can used the ids we have
      if(reconstructed_event == true){
	true_contained_start = cuts.In_FV(10,10,10,10,10,10,backtracked_start_x->at(muon_id),backtracked_start_y->at(muon_id),backtracked_start_z->at(muon_id));
	//i have to do this to fill the true_contained_end because of course the backtracked_end is not saved >:-(
	for(int j=0u; j < mc_pdg->size(); j++){
          int pdg = mc_pdg->at(j);
          if(std::abs(pdg) == 13){
	    true_contained_end =  cuts.In_FV(0,0,0,0,0,0,mc_endx->at(j),mc_endy->at(j),mc_endz->at(j));
	  }
	}
	vmuon_denom.SetXYZ(backtracked_px->at(muon_id),backtracked_py->at(muon_id),backtracked_pz->at(muon_id));
	vleading_denom.SetXYZ(backtracked_px->at(leading_proton_id),backtracked_py->at(leading_proton_id),backtracked_pz->at(leading_proton_id));
	vrecoil_denom.SetXYZ(backtracked_px->at(recoil_proton_id),backtracked_py->at(recoil_proton_id),backtracked_pz->at(recoil_proton_id));
	
	//if it is not a reco event, we have to find the ids ourselves
      } else{
	for(int j=0u; j < mc_pdg->size(); j++){
	  int pdg = mc_pdg->at(j);
	  TVector3 true_mom_vector(mc_px->at(j),mc_py->at(j),mc_pz->at(j));
	  if(std::abs(pdg) == 13){
	    true_contained_start = cuts.In_FV(10,10,10,10,10,10,mc_vx->at(j),mc_vy->at(j),mc_vz->at(j));                                                                                                         
	    true_contained_end =  cuts.In_FV(0,0,0,0,0,0,mc_endx->at(j),mc_endy->at(j),mc_endz->at(j));   
	    vmuon_denom.SetXYZ(mc_px->at(j),mc_py->at(j),mc_pz->at(j));
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
	  int leading_id_denom = mc_protons_id[0]; //leading proton id                                                                                                                                                               
	  int recoil_id_denom = mc_protons_id[1]; //recoil proton id                                                                                                                                                                 
	  vleading_denom.SetXYZ(mc_px->at(leading_id_denom),mc_py->at(leading_id_denom),mc_pz->at(leading_id_denom));
	  vrecoil_denom.SetXYZ(mc_px->at(recoil_id_denom),mc_py->at(recoil_id_denom),mc_pz->at(recoil_id_denom));
	} //if mc_proton_mom_size == 2

      } //if not a reco event

      for(size_t j=0u; j < mc_pdg->size(); j++){
	TVector3 backtracker_mom_vector(mc_px->at(j),mc_py->at(j),mc_pz->at(j)); //mc momentum of particular mc particle
	bool contained_start = cuts.In_FV(10,10,10,10,10,10,mc_vx->at(j),mc_vy->at(j),mc_vz->at(j)); //is the particle start contained in FV
	bool contained_end = cuts.In_FV(10,10,10,10,10,10,mc_endx->at(j),mc_endy->at(j),mc_endz->at(j));
	Fill_Efficiency_Thresholds(true, mc_pdg->at(j), contained_start, contained_end, vleading_denom, vrecoil_denom, backtracker_mom_vector,mc_wgt*pot_wgt); //filling the efficiency to test the thresholds
      } //end of for loop over the mc particles
      Fill_Efficiency_XSec(true,true_contained_start,true_contained_end,vmuon_denom,vleading_denom,vrecoil_denom,mc_wgt*pot_wgt);
    } //end of loop over true cc2p event
  
    //Now to Fill the Numerator of the Efficiency
    /////////////////////////////////////////////////////////////////////
    bool true_contained_start_num;
    bool true_contained_end_num;
    TVector3 vmuon_num;
    TVector3 vleading_num;
    TVector3 vrecoil_num;
    
    if(reconstructed_event == true && good_muon_mom == true && true_cc2p == true){
      num++;
      true_contained_start_num = cuts.In_FV(10,10,10,10,10,10,backtracked_start_x->at(muon_id),backtracked_start_y->at(muon_id),backtracked_start_z->at(muon_id));
      //i have to do this to fill the true_contained_end because of course the backtracked_end is not saved >:-(                                                                                                                     
      for(int j=0u; j < mc_pdg->size(); j++){
	int pdg = mc_pdg->at(j);
	if(std::abs(pdg) == 13){
	  true_contained_end_num =  cuts.In_FV(0,0,0,0,0,0,mc_endx->at(j),mc_endy->at(j),mc_endz->at(j));
	}
      }
      vmuon_num.SetXYZ(backtracked_px->at(muon_id),backtracked_py->at(muon_id),backtracked_pz->at(muon_id));
      vleading_num.SetXYZ(backtracked_px->at(leading_proton_id),backtracked_py->at(leading_proton_id),backtracked_pz->at(leading_proton_id));
      vrecoil_num.SetXYZ(backtracked_px->at(recoil_proton_id),backtracked_py->at(recoil_proton_id),backtracked_pz->at(recoil_proton_id));

      for(size_t j=0u; j < mc_pdg->size(); j++){
        TVector3 backtracker_mom_vector(mc_px->at(j),mc_py->at(j),mc_pz->at(j)); //mc momentum of particular mc particle                                                                                  
	bool contained_start = cuts.In_FV(10,10,10,10,10,10,mc_vx->at(j),mc_vy->at(j),mc_vz->at(j)); //is the particle start contained in FV    
	bool contained_end = cuts.In_FV(0,0,0,0,0,0,mc_endx->at(j),mc_endy->at(j),mc_endz->at(j)); //is the particle end contained in the detector                                                                                 
	Fill_Efficiency_Thresholds(true, mc_pdg->at(j), contained_start, contained_end, vleading_num, vrecoil_num, backtracker_mom_vector,mc_wgt*pot_wgt); //filling the efficiency to test the thresholds                                                             
      } //end of for loop over the mc particles
      Fill_Matrices(vMuon,vmuon_num,vLead,vleading_num,vRec,vrecoil_num,muon_start_contained,true_contained_start_num,muon_end_contained,true_contained_end_num, mc_wgt*pot_wgt);
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
  std::vector<int> cut_values = {static_cast<int>(nentries),fvcntr,threepfps,threetrkcntr, threetrk_connected, pid, reco_muon_mom_cut};
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
  std::cout << "[ANALYZER] Number of Events with Reco. Muon Momentum above 0.1 GeV and below 2.5 GeV: "<<reco_muon_mom_cut<<" Fraction of Total: "<<float(100.*float(reco_muon_mom_cut)/float(nentries))<<"%"<<std::endl;
  std::cout <<"-----CLOSING TIME. YOU DON'T HAVE TO GO HOME, BUT YOU CAN'T STAY HERE-----"<<std::endl;
   
  std::cout << "[MC] Initial Number of Events: "<<nentries<<std::endl;
  std::cout << "[MC] Number of Generated CC2p0pi Events: "<<cc2p0pi[0]<<" Fraction of the Total: "<<float(100.*(float(cc2p0pi[0])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC] Number of Denominator Entries Events: "<<denom<<" Fraction of the Total: "<<float(100.*(float(denom)/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC] Number of Generated anbd Reconstructed CC2p0pi Events: "<<cc2p0pi[6]<<" Fraction of the Total: "<<float(100.*(float(cc2p0pi[6])/float(nentries)))<<"%"<<std::endl;
  std::cout << "[MC] Number of Numerator Entries Events: "<<num<<" Fraction of the Total: "<<float(100.*(float(num)/float(nentries)))<<"%"<<std::endl;

  std::cout<<"God dammit: "<<god_dammit<<std::endl;

  std::cout<<"Total Protons: "<<total_protons<<std::endl;
  std::cout<<"Contained Protons: "<<contain<<std::endl;
  std::cout<<"Uncontained Protons: "<<uncontain<<std::endl;

  std::cout<<"Contained: "<<contained<<std::endl;
  std::cout<<"Uncontained: "<<uncontained<<std::endl;
  std::cout<<"Denom Contained: "<<denom_contained<<std::endl;
  std::cout<<"Denom Uncontained: "<<denom_uncontained<<std::endl;
  std::cout<<"Num Contained: "<<num_contained<<std::endl;
  std::cout<<"Num Uncontained: "<<num_uncontained<<std::endl;

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
