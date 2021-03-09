#define twoproton_pelee_overlay_cxx
#include "twoproton_pelee_overlay.h"
#include "constants.h"
#include <chrono>
using namespace Constants;
using namespace std::chrono;

void twoproton_pelee_overlay::Loop()
{
  auto start = high_resolution_clock::now(); 

  //Define objects of classes
  ////////////////////////////
  histogram_funcs hist; //histogram_funcs.h
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
	if ( mom > MUON_MOM_CUT ) {
	  mc_n_threshold_muon++;
	}
      } else if ( pdg == 2212 ) {
	double mom = cuts.real_sqrt( std::pow(energy, 2) - std::pow(MASS_PROTON, 2) );
	if ( mom > PROTON_MOM_CUT ) 
	  mc_n_threshold_proton++;
      } else if ( pdg == 111 ) {
	double mom = cuts.real_sqrt( std::pow(energy, 2) - std::pow(MASS_PION0, 2) );
	if ( mom > PION0_MOM_CUT) 
	mc_n_threshold_pion0++;
      } else if ( std::abs(pdg) == 211 ) {
	double mom = cuts.real_sqrt( std::pow(energy, 2) - std::pow(MASS_PIONPM, 2) );
	if ( mom > CHARGED_PI_MOM_CUT ) {
	  mc_n_threshold_pionpm++;
	}
      }
    }

    //making sure that the mc_pdg and npfps are the same length cause fuck me
    //////////////////////////////////////////////////////////////////////////
    std::vector<int> testVector;
    std::vector<double> mc_mom_vector;
    for(int i = 0; i < mc_pdg->size(); i++){
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

    //studd for efficiency
    std::vector<double> mc_proton_mom;
    if(ccnc == 0 && abs(nu_pdg) == 14 && nproton == 2 && nmuon == 1 && npion == 0 && npi0 == 0){

      //grabbing proton momentum first
      for(int j=0; j < n_pfps; j++){
	if (std::abs(backtracked_pdg->at(j)) == 2212){
	  TVector3 backtracker_mom_vector(backtracked_px->at(j),backtracked_py->at(j),backtracked_pz->at(j));
	  mc_proton_mom.push_back(backtracker_mom_vector.Mag());
	}
      }
      sort(mc_proton_mom.begin(),mc_proton_mom.end());

      //now for the rest
      for(int j=0; j < n_pfps; j++){
	TVector3 backtracker_mom_vector(backtracked_px->at(j),backtracked_py->at(j),backtracked_pz->at(j));
	if(std::abs(backtracked_pdg->at(j)) == 13){
	  h_mom_threshold_denom[0]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //muon all  
	  bool contained_start = cuts.In_FV(10,10,10,10,10,10,trk_start_x_v->at(j),trk_start_y_v->at(j),trk_start_z_v->at(j));
          bool contained_end = cuts.In_FV(10,10,10,10,10,10,trk_end_x_v->at(j),trk_end_y_v->at(j),trk_end_z_v->at(j));
	  if(contained_start == true && contained_end == true){
	    h_mom_threshold_denom[1]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //muon contained
	  }else if(contained_start == true && contained_end == false){
	    h_mom_threshold_denom[2]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //muon uncontained 
	  } 
	}else if (std::abs(backtracked_pdg->at(j)) == 2212){
	  h_mom_threshold_denom[3]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //proton all
	  if(mc_proton_mom.size() == 2){
	    if(backtracker_mom_vector.Mag() == mc_proton_mom[1]){
	      h_mom_threshold_denom[4]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //proton leading
	    } else if (backtracker_mom_vector.Mag() == mc_proton_mom[0]){
	      h_mom_threshold_denom[5]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //proton recoil 
	    }
	  }
	}else if (backtracked_pdg->at(j) == 211){
	  h_mom_threshold_denom[6]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //pion+
	}else if (backtracked_pdg->at(j) == -211){
	  h_mom_threshold_denom[7]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //pion-
	}else if(std::abs(backtracked_pdg->at(j)) == 111){
	  h_mom_threshold_denom[8]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //pion0  
	} //end of else if
      } //end of for
    } //end of if

    mc_proton_mom.clear();

    //2) There are exactly 3 PFP's in the Event 
    ///////////////////////////////////////////////////////
    if(n_pfps != 3) continue;
    threepfps++;
    Fill_Histograms_Mine(2, pot_wgt*mc_wgt, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0,mc_n_threshold_pionpm,cuts.fv);
    Fill_Histograms_Raquel(2, pot_wgt*mc_wgt, cuts.fv);

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

    //Three Tracks connected to the vertex
    if(y1 != 3) continue;
    threetrk_connected++;
    Fill_Histograms_Mine(4, pot_wgt*mc_wgt, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0,mc_n_threshold_pionpm,cuts.fv);
    Fill_Histograms_Raquel(4, pot_wgt*mc_wgt, cuts.fv);
    for(int i = 0; i < n_pfps; i++){
      int track_pdg = testVector.at(i);
      bool contained_start = cuts.In_FV(10,10,10,10,10,10,trk_start_x_v->at(i),trk_start_y_v->at(i),trk_start_z_v->at(i));
      bool contained_end = cuts.In_FV(10,10,10,10,10,10,trk_end_x_v->at(i),trk_end_y_v->at(i),trk_end_z_v->at(i));	
      Fill_Track_Plots(i,track_pdg,contained_start,contained_end,pot_wgt*mc_wgt);
    }
    
    //5) PID: One track with PID > 0.6 and 2 tracks with PID < 0.6
    //////////////////////////////////////////////////////////////
    if(muons != 1 && protons != 2) continue;
    pid++;
    Fill_Histograms_Mine(5, pot_wgt*mc_wgt, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0,mc_n_threshold_pionpm,cuts.fv);
    Fill_Histograms_Raquel(5, pot_wgt*mc_wgt, cuts.fv);
    
    //filling the numerator of the efficiency
    if(ccnc == 0 && abs(nu_pdg) == 14 && nproton == 2 && nmuon == 1 && npion == 0 && npi0 == 0){

      //grabbing proton momentum first                                                                                                                                                                                               
      for(int j=0; j < n_pfps; j++){
	if (std::abs(backtracked_pdg->at(j)) == 2212){
	  TVector3 backtracker_mom_vector(backtracked_px->at(j),backtracked_py->at(j),backtracked_pz->at(j));
	  mc_proton_mom.push_back(backtracker_mom_vector.Mag());
	}
      }
      sort(mc_proton_mom.begin(),mc_proton_mom.end());

      for(int j=0; j < n_pfps; j++){
	TVector3 backtracker_mom_vector(backtracked_px->at(j),backtracked_py->at(j),backtracked_pz->at(j));
	if(std::abs(backtracked_pdg->at(j)) == 13){
	  h_mom_threshold_num[0]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //muon all 
	  bool contained_start = cuts.In_FV(10,10,10,10,10,10,trk_start_x_v->at(j),trk_start_y_v->at(j),trk_start_z_v->at(j));
	  bool contained_end = cuts.In_FV(10,10,10,10,10,10,trk_end_x_v->at(j),trk_end_y_v->at(j),trk_end_z_v->at(j));
	  if(contained_start == true && contained_end == true){
	    h_mom_threshold_num[1]->Fill(backtracker_mom_vector.Mag(),1); //muon contained
	  }else if(contained_start == true && contained_end == false){
	    h_mom_threshold_num[2]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //muon uncontained 
	  } 
	}else if (std::abs(backtracked_pdg->at(j)) == 2212){
	  h_mom_threshold_num[3]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //proton all
	  if(mc_proton_mom.size() == 2){
	    if(backtracker_mom_vector.Mag() == mc_proton_mom[1]){
	      h_mom_threshold_num[4]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //proton leading                                                                                                                                  
	    } else if (backtracker_mom_vector.Mag() == mc_proton_mom[0]){
	      h_mom_threshold_num[5]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //proton recoil                                                                                                                                   
	    }
	  }
	}else if (backtracked_pdg->at(j) == 211){
	  h_mom_threshold_num[6]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //pion+
	}else if (backtracked_pdg->at(j) == -211){
	  h_mom_threshold_num[7]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //pion-
	}else if(std::abs(backtracked_pdg->at(j)) == -111){
	  h_mom_threshold_num[8]->Fill(backtracker_mom_vector.Mag(),mc_wgt); //pion0  
	} //end of else if
      } //end of for
    } //end of the if

    mc_proton_mom.clear();

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

      if(trk_pid > cuts.PID_CUT) {
	muon_id = trk_id - 1;
      }
      if(trk_pid < cuts.PID_CUT){
	proton_id_vector.push_back(trk_id);
      }
    }

    if(proton_id_vector[0] == 0 || proton_id_vector[1] ==0) continue; 
    uhoh++;

    float mom0 = trk_energy_proton_v->at(proton_id_vector[0]-1);
    float mom1 = trk_energy_proton_v->at(proton_id_vector[1]-1);
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
    bool muon_end_contained = cuts.In_FV(10,10,10,10,10,10,trk_sce_end_x_v->at(muon_id),trk_sce_end_y_v->at(muon_id),trk_sce_end_z_v->at(muon_id));
    double EMuon = 0;
    if(muon_start_contained == true && muon_end_contained == true){
      EMuon = std::sqrt(std::pow(trk_range_muon_mom_v->at(muon_id),2)+std::pow(MASS_MUON,2)) - MASS_MUON;
      vMuon.SetMag(trk_range_muon_mom_v->at(muon_id));
    } else if (muon_start_contained == true && muon_end_contained == false){
      EMuon = std::sqrt(std::pow(trk_mcs_muon_mom_v->at(muon_id),2)+std::pow(MASS_MUON,2)) - MASS_MUON;
      vMuon.SetMag(trk_mcs_muon_mom_v->at(muon_id));
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

    Fill_Histograms_Particles(mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0, mc_n_threshold_pionpm, cuts.fv, vMuon, muon, vLead, lead, vRec, rec, mc_wgt*pot_wgt);
    Fill_Histograms_Particles_Raquel(vMuon, muon, vLead, lead, vRec, rec, mc_wgt*pot_wgt, cuts.fv);

    std::cout<<"After Filling Histograms"<<std::endl;

    //Make sure to clean up before you finish
    proton_id_vector.clear();
    testVector.clear();
    mc_mom_vector.clear();
    events_remaining++;

  } //end of Loop over events

  //Before we finish, we need to make the efficiency and purity plots:
  ///////////////////////////////////////////////////////////////////
  std::vector<int> cut_values = {static_cast<int>(nentries),fvcntr,threepfps,threetrkcntr, threetrk_connected, pid};
  for(int i = 0; i < number; i++){
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
