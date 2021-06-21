#define twoproton_pelee_ext_cxx
#include "twoproton_pelee_ext.h"
#include "histogram_funcs.h"
#include "constants.h"
#include <chrono> 
using namespace Constants;
using namespace std::chrono;

void twoproton_pelee_ext::Loop()
{
  auto start = high_resolution_clock::now();

  //Define objects of classes
  ////////////////////////////
  histogram_funcs hist; //histogram_funcs.h
  helper_funcs cuts; //helper_funcs.h    
  
  //Making a new Root File that will contain all the histograms that we will want to plot and file with good RSEs:                      
  ///////////////////////////////////////////////////////////////////////////////////////                          
  Which_Run();
  TFile *tfile = new TFile(Form("root_files/%s/histograms_pelee_ext.root",directory),"RECREATE"); //root file
  ofstream myfile;//File that will contain RSE of good events                                                                 
  myfile.open(Form("lists/%s/files_filtered_ext.list",directory));
  myfile<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;

  //Define all the histograms I am going to fill
 ////////////////////////////////////////////                                                                                
  hist.Define_Histograms("ext");

  //Random counters
  /////////////////
  int total_lead = 0;
  int flip_lead = 0;
  int total_recoil = 0;
  int flip_recoil = 0;

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

    std::cout<<"-----------------------------------"<<std::endl;
    std::cout<<"BEGINNING TO PROCESS RUN: " <<run << "  SUBRUN: "<< sub << "  EVENT: " << evt <<std::endl;
    std::cout<<"-----------------------------------"<<std::endl;

    //Filling histograms before any selection is made
    ////////////////////////////////////////////////
    hist.Fill_Histograms(0, TVector3(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z), CosmicIP, topological_score,pot_wgt);
    
    //Just casually checking how many neutrino slices we have
    if(nslice == 0){
      neutrinos_0++;
    }else if(nslice == 1){
      neutrinos_1++;
    }else{
      neutrinos_else++;
    }

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
    hist.Fill_Histograms(1, TVector3(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z), CosmicIP, topological_score,pot_wgt);

    //2) There are exactly 3 PFP's in the Event
    ////////////////////////////////////////////
    if(n_pfps != 3) continue;
    threepfps++;
    hist.Fill_Histograms(2, TVector3(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z), CosmicIP, topological_score,pot_wgt);
    for(int i = 0; i < n_pfps; i++){
      hist.h_track[0][0]->Fill(trk_score_v->at(i),pot_wgt);
      hist.h_track[1][0]->Fill(trk_distance_v->at(i),pot_wgt);
      hist.h_track[2][0]->Fill(trk_len_v->at(i),pot_wgt);
      hist.h_track[3][0]->Fill(trk_llr_pid_score_v->at(i),pot_wgt);
    }
    
    //3) Require that the 3 PFP's are tracks. Defined to have a track Score above 0.8
    /////////////////////////////////////////////////////////////////////////////////
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
      if(track_pid >= PID_CUT){
	muons++;
      }
      if(track_pid < PID_CUT){
	protons++;                                             
      }  
    }                                                                                                     

    //3 Pfps
    if(y != 3) continue;
    threetrkcntr++;
    hist.Fill_Histograms(3, TVector3(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z),CosmicIP, topological_score, pot_wgt);
    for(int i = 0; i < n_pfps; i++){
      hist.h_track[0][1]->Fill(trk_score_v->at(i),pot_wgt);
      hist.h_track[1][1]->Fill(trk_distance_v->at(i),pot_wgt);
      hist.h_track[2][1]->Fill(trk_len_v->at(i),pot_wgt);
      hist.h_track[3][1]->Fill(trk_llr_pid_score_v->at(i),pot_wgt);
    }

    //3 PFPs attached to the verted
    if(y1 != 3) continue;//three tracks connected to vertex
    threetrk_connected++;
    hist.Fill_Histograms(4, TVector3(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z),CosmicIP, topological_score, pot_wgt);
    for(int i = 0; i < n_pfps; i++){
      hist.h_track[0][2]->Fill(trk_score_v->at(i),pot_wgt);
      hist.h_track[1][2]->Fill(trk_distance_v->at(i),pot_wgt);
      hist.h_track[2][2]->Fill(trk_len_v->at(i),pot_wgt);
      hist.h_track[3][2]->Fill(trk_llr_pid_score_v->at(i),pot_wgt);
    }
    
    //5) PID: One track with PID > 0.6 and 2 tracks with PID < 0.6
    //////////////////////////////////////////////////////////////
    if(muons != 1 && protons != 2) continue;
    pid++;
    hist.Fill_Histograms(5, TVector3(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z),CosmicIP, topological_score, pot_wgt);

    ///////////////////////////////////////////////////////////////////////////////////////////
    //Huzzah! We are done with the inital selection. Now to make some particle specific plots:
    //////////////////////////////////////////////////////////////////////////////////////////
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
    } else if (muon_start_contained == true && muon_end_contained == false){
      EMuon = std::sqrt(std::pow(trk_mcs_muon_mom_v->at(muon_id),2)+std::pow(MASS_MUON,2)) - MASS_MUON;
      vMuon.SetMag(trk_mcs_muon_mom_v->at(muon_id));
    }
    vMuon.SetTheta(trk_theta_v->at(muon_id));
    vMuon.SetPhi(trk_phi_v->at(muon_id));
    TLorentzVector muon(vMuon[0],vMuon[1],vMuon[2],EMuon);

    //We will need these to flip the protons
    TVector3 nu_vtx_reco(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z);
    std::cout<<"Location of Reco Neutrino Vertex: ("<<nu_vtx_reco[0]<<","<<nu_vtx_reco[2]<<","<<nu_vtx_reco[2]<<")"<<std::endl;

    //Leading Proton
    float lead_track_start_distance_reco = trk_distance_v->at(leading_proton_id); //distance from start to vertex: reconstructed 
    TVector3 lead_track_end_reco(trk_sce_end_x_v->at(leading_proton_id),trk_sce_end_y_v->at(leading_proton_id),trk_sce_end_z_v->at(leading_proton_id)); //leading track end reco
    lead_track_end_reco -= nu_vtx_reco;
    double lead_track_end_distance_reco = lead_track_end_reco.Mag(); //distance from end to vertex: reconstructed                          

    TVector3 vLead(1,1,1);
    float ELead = trk_energy_proton_v->at(leading_proton_id);
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

    //Recoil Proton
    float recoil_track_start_distance_reco = trk_distance_v->at(recoil_proton_id); //distance from start to vertex: reconstructed        
    TVector3 recoil_track_end_reco(trk_sce_end_x_v->at(recoil_proton_id),trk_sce_end_y_v->at(recoil_proton_id),trk_sce_end_z_v->at(recoil_proton_id));
    recoil_track_end_reco -= nu_vtx_reco;
    double recoil_track_end_distance_reco = recoil_track_end_reco.Mag(); //distance from end to vertex: reconstructed                      

    TVector3 vRec(1,1,1);
    float ERec = trk_energy_proton_v->at(recoil_proton_id);
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


    //We had to add another cut: The muon and protons must have reconstructed momentum within the thresholdss defined
    // this is the only way to get the closure test to work
    if(vMuon.Mag() < MUON_MOM_CUT_LOW || vMuon.Mag() > MUON_MOM_CUT_HIGH) continue;
    reco_muon_mom_cut++;

    hist.Fill_Particles(vMuon,muon,vLead,lead,vRec,rec,pot_wgt);

    //Make sure to clean up before you finish
    proton_id_vector.clear();
    events_remaining++;

   } //end of Loop over events

  std::cout<<"-----MODULE SUMMARY-----"<<std::endl;
  std::cout << "[ANALYZER] Initial Number of Events: "<<nentries<<" Fraction of Total: "<<float(100.*float(nentries)/float(nentries))<<"%"<<std::endl;
  std::cout << "[ANALYZER] Number of Events with Vertex in FV: "<<fvcntr<<" Fraction of Total: "<<float(100.*float(fvcntr)/float(nentries))<<"%"<<std::endl;
  std::cout << "[ANALYZER] Number of Events with 3 PFPs: "<<threepfps<<" Fraction of Total: "<<float(100.*float(threepfps)/float(nentries))<<"%"<<std::endl;
  std::cout << "[ANALYZER] Number of Events with 3 Tracks: "<<threetrkcntr<<" Fraction of Total: "<<float(100.*float(threetrkcntr)/float(nentries))<<"%"<<std::endl;
  std::cout << "[ANALYZER] Number of Events with 3 Tracks Connected to Vertex: "<<threetrk_connected<<" Fraction of Total: "<<float(100.*float(threetrk_connected)/float(nentries))<<"%"<<std::endl; 
  std::cout << "[ANALYZER] Number of Events with 1 Muon and 2 Protons: "<<pid<<" Fraction of Total: "<<float(100.*float(pid)/float(nentries))<<"%"<<std::endl;
  std::cout << "[ANALYZER] Number of Events with Reco. Muon Momentum above 0.1 GeV and below 2.5 GeV: "<<reco_muon_mom_cut<<" Fraction of Total: "<<float(100.*float(reco_muon_mom_cut)/float(nentries))<<"%"<<std::endl;
  std::cout << "[ANALYZER] Sanity Check of the Total Number of Events Remaining: "<<events_remaining<<std::endl;
  std::cout <<"-----CLOSING TIME. YOU DON'T HAVE TO GO HOME, BUT YOU CAN'T STAY HERE-----"<<std::endl;

  std::cout<<"Neutrinos 0: "<<neutrinos_0<<std::endl;
  std::cout<<"Neutrinos 1: "<<neutrinos_1<<std::endl;
  std::cout<<"Neutrinos Else: "<<neutrinos_else<<std::endl;

  std::cout<<"Total Number of Lead Protons: "<<total_lead<<std::endl;
  std::cout<<"Flip Lead: "<<flip_lead<<std::endl;
  std::cout<<"Total Number of Recoil Protons: "<<total_recoil<<std::endl;
  std::cout<<"Flip Recoil: "<<flip_recoil<<std::endl;

   //Don't forget to write all of your histograms before you leave!                                                                       
   ///////////////////////////////////////////////////////////////                                                                 
   tfile->cd();
   hist.Write_Histograms(); //function that writes all our histograms                                                              
   tfile->Close(); //write the root file that contains our histograms                                                         
   myfile.close(); //Write the file that contains the RSE of good events                                                     

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<minutes>(stop - start); 
  std::cout<<"Program Run Time: "<<duration.count()<<std::endl;

} //end of progrm
