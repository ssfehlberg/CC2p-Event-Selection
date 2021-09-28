#define twoproton_pelee_cxx
#include "twoproton_pelee.h"

void twoproton_pelee::Loop()
{

  //Define objects of classes
  ////////////////////////////
  Selection cuts; //cuts.h: contains all the cuts I will apply

  //Making a new Root File that will contain all the histograms that we will want to plot and files with good RSEs:
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TFile* tfile;
  if(use_xsec_binning == true){
     tfile = new TFile(Form("root_files/%s/histograms_pelee_xsec_%s.root",directory,sample),"RECREATE");
  } else {
    tfile = new TFile(Form("root_files/%s/histograms_pelee_%s.root",directory,sample),"RECREATE");
  }

  //Open files to contain the RSE of Good events
  ////////////////////////////////////////////////
  ofstream myfile;//File that will contain RSE of good events                                          
  myfile.open(Form("lists/%s/%s_selected_events.csv",directory,sample));
  myfile<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;

  ofstream cc2p; //File that will contain good cc2p events
  cc2p.open(Form("lists/%s/%s_selected_cc2p_events.csv",directory,sample));
  cc2p<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;

  //Define all the histograms I am going to Fill
  /////////////////////////////////////////////////////////////
  Define_Histograms();

  //Make sure to have the event weight double ready
  ///////////////////////////////////////////////
  double event_weight; //event weight
  double mc_wgt; //MC weight. Dirt Only

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

    //Need to set the event weight for samples
    ////////////////////////////////////////////////
    if(_debug) std::cout<<"GET THE MC WEIGHT"<< std::endl;

    if(std::isfinite(weightTune) && weightTune <= 100.) {
      mc_wgt = weightSplineTimesTune;
    } else {
      mc_wgt = 1 * weightSpline;
    }

    event_weight = pot_weight;// * mc_wgt;

    if(_debug) std::cout<<"Value of the MC Weight: "<< mc_wgt<< std::endl;
    if(_debug) std::cout<<"Value of the POT Weight: "<<pot_weight << std::endl;
    if(_debug) std::cout<<"Value of Event Weight: "<< event_weight << std::endl;

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

    //Fill the Overlay FV Bool
    //////////////////////////
    cuts.Overlay_In_FV(10,10,10,10,10,10,true_nu_vtx_sce_x,true_nu_vtx_sce_y,true_nu_vtx_sce_z); //overlay FV requirment  

    //Checking how many nue's & neutrino slices we have
    /////////////////////////////////
    if(_debug) std::cout<<"CHECKING THE NUMBER OF NEUTRINO SLICES"<< std::endl;

    if(nu_pdg == 12) nue++;
    if(nslice == 0){
      neutrinos_0++;
    }else if(nslice == 1){
      neutrinos_1++;
    }else{
      neutrinos_else++;
    }

    //Create a TVector3 of the reco and true vertex
    /////////////////////////////////////////////////
    TVector3 reco_nu_vtx(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z);
    TVector3 true_nu_vtx(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z);
    TVector3 true_nu_vtx_sce(true_nu_vtx_sce_x,true_nu_vtx_sce_y,true_nu_vtx_sce_z);

    if(_debug) std::cout<<"Location of the Reconstructed Vertex: "<<reco_nu_vtx[0] << ", "<< reco_nu_vtx[1] << ", " <<reco_nu_vtx[2] << std::endl;
    if(_debug) std::cout<<"Location of the True Vertex: "<<true_nu_vtx[0] << ", "<< true_nu_vtx[1] << ", " <<true_nu_vtx[2] << std::endl;
    if(_debug) std::cout<<"Location of the True SCE Vertex: "<<true_nu_vtx_sce[0] << ", "<< true_nu_vtx_sce[1] << ", " <<true_nu_vtx_sce[2] << std::endl;

    //Now to apply the selection. We require the following:
    // 1) Event must be in the fiducial volume. We define this to be 10 cm from any TPC edge
    // 2) There must be exactly 3 PFPs in the event (npfps == 3), 3 tracks with track score above 0.8, 3 tracks with distance less than 4cm to the vertex, and events with exactly 1 muon and 2 protons (pid)
    //3) The particles must have good momentum. 
    ///////////////////////////////////////////////

    //1) Event must be in the FV, defined to be 10cm from any TPC edge 
    if(_debug) std::cout<<"APPLY FV CUT"<< std::endl;
    cuts.In_FV(10,10,10,10,10,10,reco_nu_vtx[0],reco_nu_vtx[1],reco_nu_vtx[2]); //returns fv bool
    if(cuts.reco_fv != true) continue;    
    fvcntr++;

    //2)  There must be exactly 3 PFPs in the event (npfps == 3), 3 tracks with track score above 0.8, 3 tracks with distance less than 4cm to the vertex, and events with exactly 1 muon and 2 protons (pid)
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

    //Here is the function that calculates this.
    cuts.Event_Selection(n_pfps, tracks_w_good_score, tracks_w_good_distance, muons, protons); //returns bool for the pfp selection + pid
    
    //3 PFPs
    ///////////////////////////////////
    if(_debug) std::cout<<"APPLY 3 PFPs"<< std::endl;
    if(cuts.three_pfps == false) continue;
    threepfps++;

    //3 Tracks
    ////////////////////////////////
    if(_debug) std::cout<<"APPLY 3 Tracks"<< std::endl;
    if(cuts.three_tracks == false) continue;
    threetrkcntr++;

    //3 Connected Tracks
    ////////////////////////
    if(_debug) std::cout<<"APPLY 3 4cm TRACKS"<< std::endl;
    if(cuts.three_tracks_4cm_vertex == false) continue;
    threetrk_connected++;

    //PID
    /////////////////////
    if(_debug) std::cout<<"FILLING HISTOGRAMS AFTER PID"<< std::endl;
    if(cuts.good_pid == false) continue;
    pid++;

    /////////////////////////////////////////////////////////////
    //Okay. Now we have to identify the leading and recoil proton
    // Determine if these are contained and define there momentum vectors
    ////////////////////////////////////////////////////////////
    if(_debug) std::cout<<"NOW TO IDENTIFY THE PARTICLES AND FILL THEIR MOMENTUM"<< std::endl;
    int muon_id;
    int leading_proton_id;
    int recoil_proton_id;
    std::vector<int> proton_id_vector;
    for(int i=0; i < trk_pfp_id_v->size(); i++){
      int trk_id = trk_pfp_id_v->at(i);
      double trk_pid = trk_llr_pid_score_v->at(i);	
      if(trk_pid >= PID_CUT && trk_pid < 1 && trk_pid > -1.0) {
	muon_id = trk_id - 1;
      }
      if(trk_pid < PID_CUT && trk_pid < 1 && trk_pid > -1.0){
	proton_id_vector.push_back(trk_id);
      }
    }

    float mom0 = trk_energy_proton_v->at(proton_id_vector[0]-1);
    float mom1 = trk_energy_proton_v->at(proton_id_vector[1]-1);

    if (abs(mom0) > abs(mom1)){
      leading_proton_id = proton_id_vector[0] - 1; //you have to do the -1 cause of course the id's are indexed at one like fucking losers
      recoil_proton_id = proton_id_vector[1] - 1;
    }else{
      leading_proton_id = proton_id_vector[1] - 1;
      recoil_proton_id = proton_id_vector[0] - 1;
    }

    if(_debug) std::cout<<"Muon ID: "<< muon_id << std::endl;
    if(_debug) std::cout<<"Leading Proton ID: "<< leading_proton_id << std::endl;
    if(_debug) std::cout<<"Recoil Proton ID: "<< recoil_proton_id << std::endl;

    //Now to get the three momentum vectors
    /////////////////////////////////////////

    //Muon
    //////////
    if(_debug) std::cout<<"---MUON---"<< std::endl;
    //First check is the muon if fully contained
    bool muon_start_contained = cuts.In_FV(10,10,10,10,10,10,trk_sce_start_x_v->at(muon_id),trk_sce_start_y_v->at(muon_id),trk_sce_start_z_v->at(muon_id)); //is the muon start within the FV?
    bool muon_end_contained = cuts.In_FV(0,0,0,0,0,0,trk_sce_end_x_v->at(muon_id),trk_sce_end_y_v->at(muon_id),trk_sce_end_z_v->at(muon_id)); //is the muon end within the detector?
    if(muon_start_contained == false) continue;
    muon_contained[0]++;
    if(muon_start_contained == true && muon_end_contained == false)continue;
    muon_contained[1]++;
    if(muon_start_contained == false && muon_end_contained == false) continue;
    muon_contained[2]++;

    //now to define the momentum
    TVector3 vMuon(1,1,1);
    vMuon.SetMag(trk_range_muon_mom_v->at(muon_id));
    double EMuon = std::sqrt(std::pow(trk_range_muon_mom_v->at(muon_id),2)+std::pow(MASS_MUON,2)) - MASS_MUON;
    vMuon.SetTheta(trk_theta_v->at(muon_id));
    vMuon.SetPhi(trk_phi_v->at(muon_id));
    TLorentzVector muon(vMuon[0],vMuon[1],vMuon[2],EMuon);
    total_muon++;

    //flip momentum if the track start and end are flipped
    float  muon_track_start_distance_reco = trk_distance_v->at(muon_id); //distance from start to vertex: reconstructed         
    TVector3 muon_track_end_reco(trk_sce_end_x_v->at(muon_id),trk_sce_end_y_v->at(muon_id),trk_sce_end_z_v->at(muon_id)); //leading track end reco                                                                                                               
    muon_track_end_reco -= reco_nu_vtx;
    double muon_track_end_distance_reco = muon_track_end_reco.Mag(); //distance from end to vertex: reconstructed   
    if(_debug) std::cout<<"Muon 4 Vector: ("<<muon[0]<<","<<muon[1]<<","<<muon[2]<<","<<muon[3]<<")"<<std::endl;
    if(_debug) std::cout<<"[Reconstructed] Muon Start Distance: "<<muon_track_start_distance_reco<<" Muon End Distance: "<<muon_track_end_distance_reco<<std::endl;
    if(_debug) std::cout<<"Muon PID Value: "<<trk_llr_pid_score_v->at(muon_id)<<std::endl;

    if(muon_track_start_distance_reco > muon_track_end_distance_reco){
      vMuon *= (-1.0); //three vector
      muon.SetPxPyPzE(vMuon[0],vMuon[1],vMuon[2],EMuon); //four vector
      flip_muon++;
    }
    if(_debug) std::cout<<"After Flipping: Muon 4 Vector: ("<<muon[0]<<","<<muon[1]<<","<<muon[2]<<","<<muon[3]<<")"<<std::endl;

    //Leading Proton
    /////////////////////
    if(_debug) std::cout<<"---LEADING---"<< std::endl;

    //first check is lead proton is fully contained
    bool lead_start_contained = cuts.In_FV(10,10,10,10,10,10,trk_sce_start_x_v->at(leading_proton_id),trk_sce_start_y_v->at(leading_proton_id),trk_sce_start_z_v->at(leading_proton_id)); //start of the lead proton within the FV
    bool lead_end_contained = cuts.In_FV(0,0,0,0,0,0,trk_sce_end_x_v->at(leading_proton_id),trk_sce_end_y_v->at(leading_proton_id),trk_sce_end_z_v->at(leading_proton_id)); //is end of the lead proton within the detector
    if(lead_start_contained == false) continue;
    lead_contained[0]++;
    if(lead_start_contained == true && lead_end_contained == false)continue;
    lead_contained[1]++;
    if(lead_start_contained == false && lead_end_contained == false) continue;
    lead_contained[2]++;

    //now define the momentum
    TVector3 vLead(1,1,1);
    float ELead = trk_energy_proton_v->at(leading_proton_id);
    vLead.SetMag(std::sqrt(std::pow(ELead + MASS_PROTON,2) - std::pow(MASS_PROTON,2)));
    vLead.SetTheta(trk_theta_v->at(leading_proton_id));
    vLead.SetPhi(trk_phi_v->at(leading_proton_id));
    TLorentzVector lead(vLead[0],vLead[1],vLead[2],ELead);
    total_lead++;

    //flip momentum if the track start and end are flipped
    float lead_track_start_distance_reco = trk_distance_v->at(leading_proton_id); //distance from start to vertex: reconstructed 
    TVector3 lead_track_end_reco(trk_sce_end_x_v->at(leading_proton_id),trk_sce_end_y_v->at(leading_proton_id),trk_sce_end_z_v->at(leading_proton_id)); //leading track end reco
    lead_track_end_reco -= reco_nu_vtx;
    double lead_track_end_distance_reco = lead_track_end_reco.Mag(); //distance from end to vertex: reconstructed
    if(_debug) std::cout<<"Leading Proton 4 Vector: ("<<lead[0]<<","<<lead[1]<<","<<lead[2]<<","<<lead[3]<<")"<<std::endl;
    if(_debug) std::cout<<"[Reconstructed] Leading Start Distance: "<<lead_track_start_distance_reco<<" Leading End Distance: "<<lead_track_end_distance_reco<<std::endl;
    if(_debug) std::cout<<"Leading PID Value: "<<trk_llr_pid_score_v->at(leading_proton_id)<<std::endl;

    if(lead_track_start_distance_reco > lead_track_end_distance_reco){
      vLead *= (-1.0); //three vector
      lead.SetPxPyPzE(vLead[0],vLead[1],vLead[2],ELead); //four vector
      flip_lead++;
    }
    if(_debug) std::cout<<"After Flipping: Leading Proton 4 Vector: ("<<lead[0]<<","<<lead[1]<<","<<lead[2]<<","<<lead[3]<<")"<<std::endl;

    //Recoil Proton
    ////////////////////////
    if(_debug) std::cout<<"---RECOIL---"<< std::endl;

    //first check if recoil proton is fully contained
    bool recoil_start_contained = cuts.In_FV(10,10,10,10,10,10,trk_sce_start_x_v->at(recoil_proton_id),trk_sce_start_y_v->at(recoil_proton_id),trk_sce_start_z_v->at(recoil_proton_id)); //start of the recoil proton within the FV                              
    bool recoil_end_contained = cuts.In_FV(0,0,0,0,0,0,trk_sce_end_x_v->at(recoil_proton_id),trk_sce_end_y_v->at(recoil_proton_id),trk_sce_end_z_v->at(recoil_proton_id)); //is end of the recoil proton within the detector                                    
    if(recoil_start_contained == false) continue;
    recoil_contained[0]++;
    if(recoil_start_contained == true && recoil_end_contained == false) continue;
    recoil_contained[1]++;
    if(recoil_start_contained == false && recoil_end_contained == false) continue;
    recoil_contained[2]++;

    //now define the momentum
    TVector3 vRec(1,1,1);
    float ERec = trk_energy_proton_v->at(recoil_proton_id);
    vRec.SetMag(std::sqrt(std::pow(ERec + MASS_PROTON,2) - std::pow(MASS_PROTON,2)));
    vRec.SetTheta(trk_theta_v->at(recoil_proton_id));
    vRec.SetPhi(trk_phi_v->at(recoil_proton_id));
    TLorentzVector rec(vRec[0],vRec[1],vRec[2],ERec);
    total_recoil++;

    //flip the momentum if the track start and end are flipped
    float recoil_track_start_distance_reco = trk_distance_v->at(recoil_proton_id); //distance from start to vertex: reconstructe
    TVector3 recoil_track_end_reco(trk_sce_end_x_v->at(recoil_proton_id),trk_sce_end_y_v->at(recoil_proton_id),trk_sce_end_z_v->at(recoil_proton_id));
    recoil_track_end_reco -= reco_nu_vtx;
    double recoil_track_end_distance_reco = recoil_track_end_reco.Mag(); //distance from end to vertex: reconstructed    
    if(_debug) std::cout<<"Recoil Proton 4 Vector: ("<<rec[0]<<","<<rec[1]<<","<<rec[2]<<","<<rec[3]<<")"<<std::endl;
    if(_debug) std::cout<<"[Reconstructed] Recoil Start Distance: "<<recoil_track_start_distance_reco<<" Recoil End Distance: "<<recoil_track_end_distance_reco<<std::endl;
    if(_debug) std::cout<<"Recoil PID Value: "<<trk_llr_pid_score_v->at(recoil_proton_id)<<std::endl;

    if(recoil_track_start_distance_reco > recoil_track_end_distance_reco){
      vRec *= (-1.0); //three vector                                                                                                         
      rec.SetPxPyPzE(vRec[0],vRec[1],vRec[2],ERec); //four vector 
      flip_recoil++;
    }
    if(_debug) std::cout<<"After Flipping: Recoil Proton 4 Vector: ("<<rec[0]<<","<<rec[1]<<","<<rec[2]<<","<<rec[3]<<")"<<std::endl;

    /////////////////////////////////////////////////////////////////////////////////////
    //3) Require the muon, leading proton, and recoil proton to be within Momentum Limits
    ////////////////////////////////////////////////////////////////////////////////////
    cuts.Reco_Momentum(vMuon, vLead, vRec);
    if(cuts.good_muon_mom == false) continue;
    reco_muon_mom_cut++;
    if(cuts.good_lead_mom == false) continue;
    reco_lead_mom_cut++;
    if(cuts.good_recoil_mom == false) continue;
    reco_recoil_mom_cut++;

    Fill_Histograms_Particles(mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0, mc_n_threshold_pionpm, cuts.true_fv, vMuon, muon, vLead, lead, vRec, rec, event_weight);
    Fill_Histograms_Particles_Raquel(vMuon, muon, vLead, lead, vRec, rec, event_weight, cuts.true_fv);

    //////////////////////////////////////////////////////////////////
    //Now to fill the Reco_Event Boolean, and fill the final counters
    /////////////////////////////////////////////////////////////////
    cuts.Reco_Event(); 

    /////////////////////////////////////////
    //Making sure to save the RSE in the Myfile
    //Also get number of remaining events
    ////////////////////////////////////
    myfile << run << " " << sub << " " << evt << " " ;
    myfile << endl;
    events_remaining++;
    
  } //END OF LOOP OVER EVENTS

  ///////////////////////////////////////////////////////////////
  // Saving the RSE for the Events in the MyFILE CSV
  //Also Making an Output CSV to save the number of events after cuts.
  ////////////////////////////////////////////////////////////////
  std::ofstream csv_file;
  csv_file.open(Form("lists/%s/%s.csv",directory,sample));
  csv_file << "Cut, Number of Events, Fraction of Total (%) \n";
  csv_file << Form("Number of Events to Begin With, %lld,  %f \n", nentries, float(100.*float(nentries)/float(nentries)));
  csv_file << Form("Number of Events with Vertex in FV, %d , %f \n",fvcntr,float(100.*float(fvcntr)/float(nentries)));
  csv_file << Form("Number of Events with 3 PFPs , %d , %f \n",threepfps,float(100.*float(threepfps)/float(nentries)));
  csv_file << Form("Number of Events with 3 Tracks , %d , %f \n",threetrkcntr,float(100.*float(threetrkcntr)/float(nentries)));
  csv_file << Form("Number of Events with 3 Tracks Connected to the Vertex, %d , %f \n",threetrk_connected,float(100.*float(threetrk_connected)/float(nentries)));
  csv_file << Form("Number of Events with 1 Muon and 2 Protons, %d , %f \n",pid,float(100.*float(pid)/float(nentries)));
  csv_file << Form("Number of Events with Reco. Muon Momentum above %f and below %f GeV/c, %d , %f \n",MUON_MOM_CUT_LOW, MUON_MOM_CUT_HIGH,reco_muon_mom_cut,float(100.*float(reco_muon_mom_cut)/float(nentries)));
  csv_file << Form("Number of Events with Reco. Leading Momentum above %f and below %f GeV/c, %d , %f \n",PROTON_MOM_CUT_LOW, PROTON_MOM_CUT_HIGH,reco_lead_mom_cut,float(100.*float(reco_lead_mom_cut)/float(nentries)));
  csv_file << Form("Number of Events with Reco. Recoil Momentum above %f and below %f GeV/c, %d , %f \n",PROTON_MOM_CUT_LOW, PROTON_MOM_CUT_HIGH,reco_recoil_mom_cut,float(100.*float(reco_recoil_mom_cut)/float(nentries)));  
  csv_file << Form("Sanity Check of the  Total Number of Events Remaining, %d \n",events_remaining);
  //csv_file.close();

  if(print_module_summary){
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
    std::cout << "[ANALYZER] Sanity Check of the Total Number of Events Remaining: "<<events_remaining<<std::endl;
    std::cout <<"-----CLOSING TIME. YOU DON'T HAVE TO GO HOME, BUT YOU CAN'T STAY HERE-----"<<std::endl;
  }

  //Create CSV file with the selected event information, MC Breakdowns
  /////////////////////////////////////////////////////////////////////
  std::ofstream overlay_csv_file;
  overlay_csv_file.open(Form("lists/%s/%s.csv",directory,sample));

  overlay_csv_file << "GENERATED EVENTS: RAQUEL \n";
  overlay_csv_file << "Cut, Number of Events, Fraction of Total (%) \n";
  overlay_csv_file << Form("Initial Number of Events, %lld, %f \n",nentries,float(100.*float(nentries)/float(nentries)));
  overlay_csv_file << Form("Number of CCQEL Events, %d, %f \n",qel_0,float(100.*(float(qel_0)/float(nentries))));
  overlay_csv_file << Form("Number of CCRES Events, %d, %f \n",res_0,float(100.*(float(res_0)/float(nentries))));
  overlay_csv_file << Form("Number of CCDIS Events, %d, %f \n",dis_0,float(100.*(float(dis_0)/float(nentries))));
  overlay_csv_file << Form("Number of CCCOH Events, %d, %f \n",coh_0,float(100.*(float(coh_0)/float(nentries))));
  overlay_csv_file << Form("Number of CCNuE Events, %d, %f \n",ccnue_raquel_0,float(100.*(float(ccnue_raquel_0)/float(nentries))));
  overlay_csv_file << Form("Number of NC Events, %d, %f \n",nc_raquel_0,float(100.*(float(nc_raquel_0)/float(nentries))));
  overlay_csv_file << Form("Number of OOFV Events, %d, %f \n",outfv_raquel_0,float(100.*(float(outfv_raquel_0)/float(nentries))));
  overlay_csv_file << Form("Number of Else Events, %d, %f \n",other_raquel_0,float(100.*(float(other_raquel_0)/float(nentries))));
  overlay_csv_file << "\n";
      
  overlay_csv_file << "SELECTED EVENTS: RAQUEL \n";
  overlay_csv_file << "Cut, Number of Events, Fraction of Total (%) \n";
  overlay_csv_file << Form("Number of Selected Events, %d, %f \n",events_remaining,float(100.*float(events_remaining)/float(events_remaining)));
  overlay_csv_file << Form("Number of CCQEL Events, %d, %f \n",qel,float(100.*(float(qel)/float(events_remaining))));
  overlay_csv_file << Form("Number of CCRES Events, %d, %f \n",res,float(100.*(float(res)/float(events_remaining))));
  overlay_csv_file << Form("Number of CCDIS Events, %d, %f \n",dis,float(100.*(float(dis)/float(events_remaining))));
  overlay_csv_file << Form("Number of CCCOH Events, %d, %f \n",coh,float(100.*(float(coh)/float(events_remaining))));
  overlay_csv_file << Form("Number of CCNuE Events, %d, %f \n",ccnue_raquel,float(100.*(float(ccnue_raquel)/float(events_remaining))));
  overlay_csv_file << Form("Number of NC Events, %d, %f \n",nc_raquel,float(100.*(float(nc_raquel)/float(events_remaining))));
  overlay_csv_file << Form("Number of OOFV Events, %d, %f \n",outfv_raquel,float(100.*(float(outfv_raquel)/float(events_remaining))));
  overlay_csv_file << Form("Number of Else Events, %d, %f \n",other_raquel,float(100.*(float(other_raquel)/float(events_remaining))));
  overlay_csv_file << "\n";

  overlay_csv_file << "GENERATED EVENTS: MINE \n";
  overlay_csv_file << "Cut, Number of Events, Fraction of Total (%) \n";
  overlay_csv_file << Form("Initial Number of Events, %lld, %f \n",nentries,float(100.*float(nentries)/float(nentries)));
  overlay_csv_file << Form("Number of CC0p0pi Events, %d, %f \n",cc0p0pi_0,float(100.*float(cc0p0pi_0)/float(nentries)));
  overlay_csv_file << Form("Number of CC1p0pi Events, %d, %f \n",cc1p0pi_0,float(100.*float(cc1p0pi_0)/float(nentries)));
  overlay_csv_file << Form("Number of CC2p0pi Events, %d, %f \n",cc2p0pi_0,float(100.*float(cc2p0pi_0)/float(nentries)));
  overlay_csv_file << Form("Number of CCNp0pi Events, %d, %f \n",ccNp0pi_0,float(100.*float(ccNp0pi_0)/float(nentries)));
  overlay_csv_file << Form("Number of CCNp1pi Events, %d, %f \n",ccNp1pi_0,float(100.*float(ccNp1pi_0)/float(nentries)));
  overlay_csv_file << Form("Number of CCNpNpi Events, %d, %f \n",ccNpNpi_0,float(100.*float(ccNpNpi_0)/float(nentries)));
  overlay_csv_file << Form("Number of CCNue Events, %d, %f \n",ccnue_0, float(100.*float(ccnue_0)/float(nentries)));
  overlay_csv_file << Form("Number of NC Events, %d, %f \n",nc_0, float(100.*float(nc_0)/float(nentries)));
  overlay_csv_file << Form("Number of OOFV Events, %d, %f \n",outfv_0, float(100.*float(outfv_0)/float(nentries)));
  overlay_csv_file << Form("Number of Other Events, %d, %f \n",other_0, float(100.*float(other_0)/float(nentries)));
  overlay_csv_file << "\n";

  overlay_csv_file << "SELECTED EVENTS: MINE \n";
  overlay_csv_file << "Cut, Number of Events, Fraction of Total (%) \n";
  overlay_csv_file << Form("Number of Selected Events, %d, %f \n",events_remaining,float(100.*float(events_remaining)/float(events_remaining)));
  overlay_csv_file << Form("Number of CC0p0pi Events, %d, %f \n",cc0p0pi,float(100.*float(cc0p0pi)/float(events_remaining)));
  overlay_csv_file << Form("Number of CC1p0pi Events, %d, %f \n",cc1p0pi,float(100.*float(cc1p0pi)/float(events_remaining)));
  overlay_csv_file << Form("Number of CC2p0pi Events, %d, %f \n",cc2p0pi,float(100.*float(cc2p0pi)/float(events_remaining)));
  overlay_csv_file << Form("Number of CCNp0pi Events, %d, %f \n",ccNp0pi,float(100.*float(ccNp0pi)/float(events_remaining)));
  overlay_csv_file << Form("Number of CCNp1pi Events, %d, %f \n",ccNp1pi,float(100.*float(ccNp1pi)/float(events_remaining)));
  overlay_csv_file << Form("Number of CCNpNpi Events, %d, %f \n",ccNpNpi,float(100.*float(ccNpNpi)/float(events_remaining)));
  overlay_csv_file << Form("Number of CCNue Events, %d, %f \n",ccnue, float(100.*float(ccnue)/float(events_remaining)));
  overlay_csv_file << Form("Number of NC Events, %d, %f \n",nc, float(100.*float(nc)/float(events_remaining)));
  overlay_csv_file << Form("Number of OOFV Events, %d, %f \n",outfv, float(100.*float(outfv)/float(events_remaining)));
  overlay_csv_file << Form("Number of Other Events, %d, %f \n",other_0, float(100.*float(other_0)/float(events_remaining)));
  overlay_csv_file << "\n";
 
  if(print_module_summary){

    std::cout<<"-----MC GENERATED SUMMARY: RAQUEL-----"<<std::endl;
    std::cout << "[MC_RAQUEL] Initial Number of Events: "<<nentries<<std::endl;
    std::cout << "[MC_RAQUEL] Number of CCQEL Events: "<<qel_0<<" Fraction of the Total: "<<float(100.*(float(qel_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC_RAQUEL] Number of CCRES Events: "<<res_0<<" Fraction of the Total: "<<float(100.*(float(res_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC_RAQUEL] Number of CCMEC Events: "<<mec_0<<" Fraction of the Total: "<<float(100.*(float(mec_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC_RAQUEL] Number of CCCOH Events: "<<coh_0<<" Fraction of the Total: "<<float(100.*(float(coh_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC_RAQUEL] Number of CCDIS Events: "<<dis_0<<" Fraction of the Total: "<<float(100.*(float(dis_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC_RAQUEL] Number of CCNue Events: "<<ccnue_raquel_0<<" Fraction of the Total: "<<float(100.*(float(ccnue_raquel_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC_RAQUEL] Number of OUTFV Events: "<<outfv_raquel_0<<" Fraction of the Total: "<<float(100.*(float(outfv_raquel_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC_RAQUEL] Number of NC Events: "<<nc_raquel_0<<" Fraction of the Total: "<<float(100.*(float(nc_raquel_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC_RAQUEL] Number of Else Events: "<<other_raquel_0<<" Fraction of the Total: "<<float(100.*(float(other_raquel_0)/float(nentries)))<<"%"<<std::endl;
    std::cout <<"-----MR. SAXOBEAT-----"<<std::endl;
    
    std::cout<<"-----MC GENERATED SUMMARY-----"<<std::endl;
    std::cout << "[MC] Initial Number of Events: "<<nentries<<std::endl;
    std::cout << "[MC] Number of CCOpOpi Events: "<<cc0p0pi_0<<" Fraction of the Total: "<<float(100.*(float(cc0p0pi_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC] Number of CC1p0pi Events: "<<cc1p0pi_0<<" Fraction of the Total: "<<float(100.*(float(cc1p0pi_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC] Number of CC2p0pi Events: "<<cc2p0pi_0<<" Fraction of the Total: "<<float(100.*(float(cc2p0pi_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC] Number of CCNp0pi Events: "<<ccNp0pi_0<<" Fraction of the Total: "<<float(100.*(float(ccNp0pi_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC] Number of CCNp1pi Events: "<<ccNp1pi_0<<" Fraction of the Total: "<<float(100.*(float(ccNp1pi_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC] Number of CCNpNpi Events: "<<ccNpNpi_0<<" Fraction of the Total: "<<float(100.*(float(ccNpNpi_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC] Number of CCNue Events: "<<ccnue_0<<" Fraction of the Total: "<<float(100.*(float(ccnue_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC] Number of OUTFV Events: "<<outfv_0<<" Fraction of the Total: "<<float(100.*(float(outfv_0)/float(nentries)))<<"%"<<std::endl; 
    std::cout << "[MC] Number of NC Events: "<<nc_0<<" Fraction of the Total: "<<float(100.*(float(nc_0)/float(nentries)))<<"%"<<std::endl;
    std::cout << "[MC] Number of Else Events: "<<other_0<<" Fraction of the Total: "<<float(100.*(float(other_0)/float(nentries)))<<"%"<<std::endl;
    std::cout <<"-----OPEN UP MY EAGER EYES! CAUSE I'M MR. BRIGHTSIDE-----"<<std::endl;

    std::cout<<"-----MC RECO'D SUMMARY: RAQUEL-----"<<std::endl;
    std::cout << "[MC_RECO_RAQUEL] Initial Number of Events That were Reconstructed: "<<events_remaining<<std::endl;
    std::cout << "[MC_RECO_RAQUEL] Number of CCQEL Events: "<<qel<<" Fraction of the Total: "<<float(100.*(float(qel)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO_RAQUEL] Number of CCRES Events: "<<res<<" Fraction of the Total: "<<float(100.*(float(res)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO_RAQUEL] Number of CCMEC Events: "<<mec<<" Fraction of the Total: "<<float(100.*(float(mec)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO_RAQUEL] Number of CCCOH Events: "<<coh<<" Fraction of the Total: "<<float(100.*(float(coh)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO_RAQUEL] Number of CCDIS Events: "<<dis<<" Fraction of the Total: "<<float(100.*(float(dis)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO_RAQUEL] Number of CCNue Events: "<<ccnue_raquel<<" Fraction of the Total: "<<float(100.*(float(ccnue_raquel)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO_RAQUEL] Number of OUTFV Events: "<<outfv_raquel<<" Fraction of the Total: "<<float(100.*(float(outfv_raquel)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO_RAQUEL] Number of NC Events: "<<nc_raquel<<" Fraction of the Total: "<<float(100.*(float(nc_raquel)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO_RAQUEL] Number of Else Events: "<<other_raquel<<" Fraction of the Total: "<<float(100.*(float(other_raquel)/float(events_remaining)))<<"%"<<std::endl;
    std::cout <<"-----ONE FOR THE DAGGER, AND ONE FOR THE ONE YOU BELIEVE!!-----"<<std::endl;
    
    std::cout<<"-----MC RECO'D SUMMARY-----"<<std::endl;
    std::cout << "[MC_RECO] Initial Number of Events That were Reconstructed: "<<events_remaining<<std::endl;
    std::cout << "[MC_RECO] Number of CCOpOpi Events: "<<cc0p0pi<<" Fraction of the Total: "<<float(100.*(float(cc0p0pi)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO] Number of CC1p0pi Events: "<<cc1p0pi<<" Fraction of the Total: "<<float(100.*(float(cc1p0pi)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO] Number of CC2p0pi Events: "<<cc2p0pi<<" Fraction of the Total: "<<float(100.*(float(cc2p0pi)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO] Number of CCNp0pi Events: "<<ccNp0pi<<" Fraction of the Total: "<<float(100.*(float(ccNp0pi)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO] Number of CCNp1pi Events: "<<ccNp1pi<<" Fraction of the Total: "<<float(100.*(float(ccNp1pi)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO] Number of CCNpNpi Events: "<<ccNpNpi<<" Fraction of the Total: "<<float(100.*(float(ccNpNpi)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO] Number of CCNue Events: "<<ccnue<<" Fraction of the Total: "<<float(100.*(float(ccnue)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO] Number of OUTFV Events: "<<outfv<<" Fraction of the Total: "<<float(100.*(float(outfv)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO] Number of NC Events: "<<nc<<" Fraction of the Total: "<<float(100.*(float(nc)/float(events_remaining)))<<"%"<<std::endl;
    std::cout << "[MC_RECO] Number of Else Events: "<<other<<" Fraction of the Total: "<<float(100.*(float(other)/float(events_remaining)))<<"%"<<std::endl;
    std::cout <<"-----NOTHING REALLY MATTERS. ANYONE CAN SEE. NOTHING REALLY MATTERS. NOTHING REALLY MATTERS TO ME-----"<<std::endl;
  
  } //end of print module summary

  if(_debug){
    std::cout<<"Neutrinos 0: "<<neutrinos_0<<std::endl;
    std::cout<<"Neutrinos 1: "<<neutrinos_1<<std::endl;
    std::cout<<"Neutrinos Else: "<<neutrinos_else<<std::endl;

    std::cout<<"Other Else: "<<other_else<<std::endl;
    std::cout<<"Neutron: "<<neutron<<std::endl;
    std::cout<<"Neutrino: "<<neutrino<<std::endl;
    std::cout<<"Zeros: "<<zeros<<std::endl;

    std::cout<<"Total Protons: "<<total_protons<<std::endl;
    std::cout<<"Contained Protons: "<<contain<<std::endl;
    std::cout<<"Uncontained Protons: "<<uncontain<<std::endl;
    std::cout<<"UhOh: "<<uhoh<<std::endl;
    
    std::cout<<"Contained: "<<contained<<std::endl;
    std::cout<<"Uncontained: "<<uncontained<<std::endl;
    std::cout<<"Denom Contained: "<<denom_contained<<std::endl;
    std::cout<<"Denom Uncontained: "<<denom_uncontained<<std::endl;
    std::cout<<"Num Contained: "<<num_contained<<std::endl;
    std::cout<<"Num Uncontained: "<<num_uncontained<<std::endl;

    std::cout<<"Total Number of Muon: "<<total_muon<<std::endl;
    std::cout<<"Flip Muon: "<<flip_muon<<std::endl;
    std::cout<<"Total Number of Lead Protons: "<<total_lead<<std::endl;
    std::cout<<"Flip Lead: "<<flip_lead<<std::endl;
    std::cout<<"Total Number of Recoil Protons: "<<total_recoil<<std::endl;
    std::cout<<"Flip Recoil: "<<flip_recoil<<std::endl;
  }

  ////////////////////////////////////////////////////////////////
  //Don't forget to write all of your histograms before you leave!                                                                 
  ////////////////////////////////////////////////////////////////                                                              
  tfile->cd();
  Write_Histograms(); //function that writes all our histograms                                                      
  tfile->Close(); //write the root file that contains our histograms                                                    
  myfile.close(); //Write the file that contains the RSE of good events                                                 
  cc2p.close();
  csv_file.close(); //Write the file that contains the module summary
  overlay_csv_file.close();
  
} //end of progrm
