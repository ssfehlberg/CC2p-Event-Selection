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
  helper_funcs cuts; //helper_funcs.h   

  //Making a new Root File that will contain all the histograms that we will want to plot and files with good RSEs:
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TFile* tfile = new TFile(Form("/uboone/data/users/sfehlber/Systematics/%s/histograms_%s_mc_eff.root",directory,sample),"RECREATE");

  //Define all the histograms I am going to fill and the mc_wgt                                
  /////////////////////////////////////////////////////////////
  Define_Histograms();

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

    //Getting the Number of Threshold Particles:
    ///////////////////////////////////////////////
    int mc_n_threshold_muon = 0;
    int mc_n_threshold_proton = 0;
    int mc_n_threshold_pion0 = 0;
    int mc_n_threshold_pionpm = 0;

    for ( size_t p = 0u; p < mc_pdg->size(); ++p ) {
      int pdg = mc_pdg->at( p );
      float energy = mc_E->at( p );
      if ( std::abs(pdg) == 13) {
	double mom = real_sqrt( std::pow(energy, 2) - std::pow(MASS_MUON, 2) );
	if(_debug) std::cout<<"Value of the Muon Momentum: "<<mom<<std::endl;
	if ( mom > MUON_MOM_CUT_LOW && mom < MUON_MOM_CUT_HIGH){
	  mc_n_threshold_muon++;
	}
      } else if (std::abs(pdg) == 2212 ) {
	double mom = real_sqrt( std::pow(energy, 2) - std::pow(MASS_PROTON, 2) );
	if(_debug) std::cout<<"Value of the Proton Momentum: "<<mom<<std::endl;
	if ( mom > PROTON_MOM_CUT_LOW && mom < PROTON_MOM_CUT_HIGH){
	  mc_n_threshold_proton++;
	}
      } else if ( pdg == 111 ) {
	  mc_n_threshold_pion0++;
      
      } else if (std::abs(pdg) == 211 ) {
	double mom = real_sqrt( std::pow(energy, 2) - std::pow(MASS_PIONPM, 2) );
	if(_debug) std::cout<<"Value of the PionPM Momentum: "<<mom<<std::endl;
	if ( mom > CHARGED_PI_MOM_CUT ) {
	  mc_n_threshold_pionpm++;
	}
      }
    }
    
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
    }else{
      true_cc2p = false;
    }

    //Reconstructed FV Stuff                                                                                                           
    ///////////////////////////////// 
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
    TVector3 nu_vtx_reco(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z);
    for(int i = 0; i < n_pfps; i ++){
      float track_score = trk_score_v->at(i);
      float track_distance = trk_distance_v->at(i);
      float track_pid = trk_llr_pid_score_v->at(i);
      TVector3 track_end(trk_sce_end_x_v->at(i),trk_sce_end_y_v->at(i),trk_sce_end_z_v->at(i)); //leading track end reco           
      track_end -= nu_vtx_reco;
      double track_end_distance = track_end.Mag();
      if(track_score >= TRACK_SCORE_CUT){
        tracks_w_good_score++;
      }
      if(track_distance <= TRACK_DIST_CUT || track_end_distance <= TRACK_DIST_CUT){
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

    //Filling a final bool that just indicates that the event is a good reco event:                                                                                                                                                       
    ///////////////////////////////////////////////////////////////////////////////                                                                                                                                                  
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
    TVector3 nu_vtx_true;

    bool muon_contained_bool;
    bool lead_contained_bool;
    bool recoil_contained_bool;

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

    TVector3 nu_vtx_true(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z);
    std::cout<<"Location of True Neutrino Vertex: ("<<nu_vtx_true[0]<<","<<nu_vtx_true[2]<<","<<nu_vtx_true[2]<<")"<<std::endl;
    std::cout<<"Location of Reco Neutrino Vertex: ("<<nu_vtx_reco[0]<<","<<nu_vtx_reco[2]<<","<<nu_vtx_reco[2]<<")"<<std::endl;

    //Muon
    //////////
    
    //First check is the muon if fully contained
    bool muon_start_contained = cuts.In_FV(10,10,10,10,10,10,trk_sce_start_x_v->at(muon_id),trk_sce_start_y_v->at(muon_id),trk_sce_start_z_v->at(muon_id)); //is the muon start within the FV?
    bool muon_end_contained = cuts.In_FV(0,0,0,0,0,0,trk_sce_end_x_v->at(muon_id),trk_sce_end_y_v->at(muon_id),trk_sce_end_z_v->at(muon_id)); //is the muon end within the detector?
    if(muon_start_contained == true && muon_end_contained == true){
      muon_contained_bool = true;    
    }else{
      muon_contained_bool = false;
    }

    //now to define the momentum
    vMuon.SetMag(trk_range_muon_mom_v->at(muon_id));
    double EMuon = std::sqrt(std::pow(trk_range_muon_mom_v->at(muon_id),2)+std::pow(MASS_MUON,2)) - MASS_MUON;
    vMuon.SetTheta(trk_theta_v->at(muon_id));
    vMuon.SetPhi(trk_phi_v->at(muon_id));
    TLorentzVector muon(vMuon[0],vMuon[1],vMuon[2],EMuon);

    //flip momentum if the track start and end are flipped
    float  muon_track_start_distance_reco = trk_distance_v->at(muon_id); //distance from start to vertex: reconstructed                                                                                                                                                   
    TVector3 muon_track_end_reco(trk_sce_end_x_v->at(muon_id),trk_sce_end_y_v->at(muon_id),trk_sce_end_z_v->at(muon_id)); //leading track end reco                                                                                                                        
    muon_track_end_reco -= nu_vtx_reco;
    double muon_track_end_distance_reco = muon_track_end_reco.Mag(); //distance from end to vertex: reconstructed   
    if(_debug) std::cout<<"Muon 4 Vector: ("<<muon[0]<<","<<muon[1]<<","<<muon[2]<<","<<muon[3]<<")"<<std::endl;
    if(_debug) std::cout<<"[Reconstructed] Muon Start Distance: "<<muon_track_start_distance_reco<<" Muon End Distance: "<<muon_track_end_distance_reco<<std::endl;
    if(_debug) std::cout<<"Muon PID Value: "<<trk_llr_pid_score_v->at(muon_id)<<std::endl;

    if(muon_track_start_distance_reco > muon_track_end_distance_reco){
      vMuon *= (-1.0); //three vector
      muon.SetPxPyPzE(vMuon[0],vMuon[1],vMuon[2],EMuon); //four vector
    }
    if(_debug) std::cout<<"After Flipping: Muon 4 Vector: ("<<muon[0]<<","<<muon[1]<<","<<muon[2]<<","<<muon[3]<<")"<<std::endl;

    //Leading Proton
    /////////////////////
    
    //first check is lead proton is fully contained
    bool lead_start_contained = cuts.In_FV(10,10,10,10,10,10,trk_sce_start_x_v->at(leading_proton_id),trk_sce_start_y_v->at(leading_proton_id),trk_sce_start_z_v->at(leading_proton_id)); //start of the lead proton within the FV
    bool lead_end_contained = cuts.In_FV(0,0,0,0,0,0,trk_sce_end_x_v->at(leading_proton_id),trk_sce_end_y_v->at(leading_proton_id),trk_sce_end_z_v->at(leading_proton_id)); //is end of the lead proton within the detector
    if(lead_start_contained == true && lead_end_contained == true){
      lead_contained_bool = true;
    }else{
      lead_contained_bool= false;
    }

    //now define the momentum
    float ELead = trk_energy_proton_v->at(leading_proton_id);
    vLead.SetMag(std::sqrt(std::pow(ELead + MASS_PROTON,2) - std::pow(MASS_PROTON,2)));
    vLead.SetTheta(trk_theta_v->at(leading_proton_id));
    vLead.SetPhi(trk_phi_v->at(leading_proton_id));
    TLorentzVector lead(vLead[0],vLead[1],vLead[2],ELead);

    //flip momentum if the track start and end are flipped
    float lead_track_start_distance_reco = trk_distance_v->at(leading_proton_id); //distance from start to vertex: reconstructed 
    TVector3 lead_track_end_reco(trk_sce_end_x_v->at(leading_proton_id),trk_sce_end_y_v->at(leading_proton_id),trk_sce_end_z_v->at(leading_proton_id)); //leading track end reco
    lead_track_end_reco -= nu_vtx_reco;
    double lead_track_end_distance_reco = lead_track_end_reco.Mag(); //distance from end to vertex: reconstructed
    if(_debug) std::cout<<"Leading Proton 4 Vector: ("<<lead[0]<<","<<lead[1]<<","<<lead[2]<<","<<lead[3]<<")"<<std::endl;
    if(_debug) std::cout<<"[Reconstructed] Leading Start Distance: "<<lead_track_start_distance_reco<<" Leading End Distance: "<<lead_track_end_distance_reco<<std::endl;
    if(_debug) std::cout<<"Leading PID Value: "<<trk_llr_pid_score_v->at(leading_proton_id)<<std::endl;

    if(lead_track_start_distance_reco > lead_track_end_distance_reco){
      vLead *= (-1.0); //three vector
      lead.SetPxPyPzE(vLead[0],vLead[1],vLead[2],ELead); //four vector
    }
    if(_debug) std::cout<<"After Flipping: Leading Proton 4 Vector: ("<<lead[0]<<","<<lead[1]<<","<<lead[2]<<","<<lead[3]<<")"<<std::endl;

    //Recoil Proton
    ////////////////////////

    //first check if recoil proton is fully contained
    bool recoil_start_contained = cuts.In_FV(10,10,10,10,10,10,trk_sce_start_x_v->at(recoil_proton_id),trk_sce_start_y_v->at(recoil_proton_id),trk_sce_start_z_v->at(recoil_proton_id)); //start of the recoil proton within the FV                              
    bool recoil_end_contained = cuts.In_FV(0,0,0,0,0,0,trk_sce_end_x_v->at(recoil_proton_id),trk_sce_end_y_v->at(recoil_proton_id),trk_sce_end_z_v->at(recoil_proton_id)); //is end of the recoil proton within the detector                                    
    if(recoil_start_contained == true && recoil_end_contained == true){
      recoil_contained_bool = true;
    }else{
      recoil_contained_bool= false;
    }

    //now define the momentum
    float ERec = trk_energy_proton_v->at(recoil_proton_id);
    vRec.SetMag(std::sqrt(std::pow(ERec + MASS_PROTON,2) - std::pow(MASS_PROTON,2)));
    vRec.SetTheta(trk_theta_v->at(recoil_proton_id));
    vRec.SetPhi(trk_phi_v->at(recoil_proton_id));
    TLorentzVector rec(vRec[0],vRec[1],vRec[2],ERec);

    //flip the momentum if the track start and end are flipped
    float recoil_track_start_distance_reco = trk_distance_v->at(recoil_proton_id); //distance from start to vertex: reconstructed                                                                                                                                         
    TVector3 recoil_track_end_reco(trk_sce_end_x_v->at(recoil_proton_id),trk_sce_end_y_v->at(recoil_proton_id),trk_sce_end_z_v->at(recoil_proton_id));
    recoil_track_end_reco -= nu_vtx_reco;
    double recoil_track_end_distance_reco = recoil_track_end_reco.Mag(); //distance from end to vertex: reconstructed    
    if(_debug) std::cout<<"Recoil Proton 4 Vector: ("<<rec[0]<<","<<rec[1]<<","<<rec[2]<<","<<rec[3]<<")"<<std::endl;
    if(_debug) std::cout<<"[Reconstructed] Recoil Start Distance: "<<recoil_track_start_distance_reco<<" Recoil End Distance: "<<recoil_track_end_distance_reco<<std::endl;
    if(_debug) std::cout<<"Recoil PID Value: "<<trk_llr_pid_score_v->at(recoil_proton_id)<<std::endl;

    if(recoil_track_start_distance_reco > recoil_track_end_distance_reco){
      vRec *= (-1.0); //three vector                                                                                                         
      rec.SetPxPyPzE(vRec[0],vRec[1],vRec[2],ERec); //four vector 
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

    bool good_lead_mom;
    if(vLead.Mag() > PROTON_MOM_CUT_LOW && vLead.Mag() < PROTON_MOM_CUT_HIGH){
      good_lead_mom = true;
    }else{
      good_lead_mom = false;
    }

    bool good_recoil_mom;
    if(vRec.Mag() > PROTON_MOM_CUT_LOW && vRec.Mag() < PROTON_MOM_CUT_HIGH){
      good_recoil_mom = true;
    }else{
      good_recoil_mom = false;
    }

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

    if(true_cc2p == true){
      
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
	
      } //end loop over 2

      //Now to get the correct weight from the map and to fill the map
      //If we are looking at the detector variations, then the weight is simply just the POT Weight
      //We can also just directly fill the histograms
      ///////////////////////////////////////////////////////////////////
      double map_weight;
      double event_weight;
      
      //DetVar
      if(response == 0){
	event_weight = apply_cv_correction_weights(pot_weight, 1.0);
	Fill_Efficiency_XSec(0,true,true_contained_start,true_contained_end,vmuon_denom,vleading_denom,vrecoil_denom,event_weight);   
	
      } else {
	auto it = weights->find(Form("%s",directory));
	//std::cout<<"it"<<(*it).first<<std::endl;
	vector <double> inVect = (*it).second; //vector that contins all the values of the systematic.  
	double number_of_universes = inVect.size();//the size of the vector indicates if it is a multisim or a unisim. 500-1000 is a multisim. 1-2 is a unisim 
	
	for (unsigned j=0; j < number_of_universes; j++){                                                                         
	  map_weight = inVect[j];
	  event_weight = apply_cv_correction_weights(pot_weight, map_weight);
	  Fill_Efficiency_XSec(j,true,true_contained_start,true_contained_end,vmuon_denom,vleading_denom,vrecoil_denom,event_weight);   
	} //end of loop over the number of universes                                                                               
      } //end of the else loop
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
    
    if(true_cc2p == true && reconstructed_event == true && muon_contained_bool == true && lead_contained_bool == true && recoil_contained_bool == true && good_muon_mom == true && good_lead_mom == true && good_recoil_mom == true){
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
      } //end of loop over 2                                                                                                                                                                                                      
      //Now to get the correct weight from the map and to fill the map
      //If we are looking at the detector variations, then the weight is simply just the POT Weight
      //We can also just directly fill the histograms
      ///////////////////////////////////////////////////////////////////
      double map_weight;
      double event_weight;
      
      //DetVar
      if(response == 0){
	event_weight = apply_cv_correction_weights(pot_weight, 1.0);
	Fill_Matrices(0,vMuon,vmuon_num,vLead,vleading_num,vRec,vrecoil_num,muon_start_contained,true_contained_start,muon_end_contained,true_contained_end, event_weight);
	Fill_Efficiency_XSec(0,false,true_contained_start,true_contained_end,vmuon_num,vleading_num,vrecoil_num,event_weight);   
	
      } else {
	auto it = weights->find(Form("%s",directory));
	//std::cout<<"it"<<(*it).first<<std::endl;
	vector <double> inVect = (*it).second; //vector that contins all the values of the systematic.  
	double number_of_universes = inVect.size();//the size of the vector indicates if it is a multisim or a unisim. 500-1000 is a multisim. 1-2 is a unisim 
	
	for (unsigned j=0; j < number_of_universes; j++){                                                                         
	  map_weight = inVect[j];
	  event_weight = apply_cv_correction_weights(pot_weight, map_weight);
	  Fill_Matrices(j,vMuon,vmuon_num,vLead,vleading_num,vRec,vrecoil_num,muon_start_contained,true_contained_start,muon_end_contained,true_contained_end, event_weight);
	  Fill_Efficiency_XSec(j,false,true_contained_start,true_contained_end,vmuon_num,vleading_num,vrecoil_num,event_weight);   
	} //end of loop over the number of universes                                                                               
      } //end of the else loop
    } //end of good reco event
    
    //Make sure to clean up before you finish
    mc_protons_id.clear();
    mc_proton_mom.clear();
    zipped.clear();
    proton_id_vector.clear();

  } //end of Loop over events


  //Don't forget to write all of your histograms before you leave!                                                                   
  ///////////////////////////////////////////////////////////////                                                                 
  tfile->cd();
  Write_Histograms(); //function that writes all our histograms                                                      
  tfile->Close(); //write the root file that contains our histograms                                                    
  
} //end of progrm
