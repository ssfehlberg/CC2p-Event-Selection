#define twoproton_eff_cxx
#include "twoproton_eff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>

void twoproton_eff::Loop()
{
  //Making a new Root File that will contain all the histograms that we will want to plot:           
  ///////////////////////////////////////////////////////////////////////////////////////
  TFile *tfile = new TFile("histograms_efficiency.root","RECREATE");

  //Define all the histograms I am going to fill
  ////////////////////////////////////////////
  const double wgt = 1; //this is NOT the POT WEIGHT CAUSE EFF! It can however be the MC_CV weight

  //List of all the 1D histograms
  std::vector<TH1D*> h_list;

  static const int number = 6;
  const char * total[number] = {"_nue","_q2","_X","_Y","_muon_mom","_proton_mom"};//"_recoil_mom","_leading_mom"};
  std::vector<float> nbins = {20,20,20,20,50,50};//,50};
  std::vector<float> xlim_low = {0.22,0,0,0,0,0};//,0};
  std::vector<float> xlim_high = {2.2,2,2,1,2,2};//,2};
  TH1D *h_num_overlay[number];
  TH1D * h_denom_overlay[number];
  
  //Efficiency Histograms
  for(int i = 0; i < number; i++){
    h_num_overlay[i] = new TH1D(Form("h_num_overlay%s",total[i]),Form("h_num_overlay%s",total[i]),nbins[i],xlim_low[i],xlim_high[i]);
    h_denom_overlay[i] = new TH1D(Form("h_denom_overlay%s",total[i]),Form("h_denom_overlay%s",total[i]),nbins[i],xlim_low[i],xlim_high[i]); 
    h_list.push_back(h_num_overlay[i]);
    h_list.push_back(h_denom_overlay[i]);
  }

  //Defining all the constans we will use later
  //////////////////////////////
  bool _debug = false; //debug statements

  //Counters
  int fvcntr = 0; //Number of events with reconstructed vertex within the FV
  int truthcntr = 0; //Requiring reco-truth vertex to be less than 5cm
  int isfromnucntr = 0; //how many pfp's are from the neutrino slice
  int has3pfp = 0; //how many events have exactly 3 pfps
  int has0shower = 0;//how many events has 0 showers (i.e. three tracks)
  int threetrkcntr = 0; //Number of events with three tracks    
  int vectorsize3 = 0; //Number of events with 3 tracks whose start is < 5cm from the reco vertex
  int secondtrkgood = 0; //Number of events where the second shortest/longest track is contained
  int shortesttrkgood=0; //Number of events where the shortest track is contained
  int pid_cut0 = 0; //sanity pid cut
  int pid_cut1 = 0; //sannity pid cut
  int proton_cut = 0; //number of events with 2 protons
  int muon_cut = 0; //number of events with 1 muon in final state
  int mom_cut = 0; //Number of events after a momentum threshold cut
  int events_remaining = 0; //sanity check for number of events remaining                                            
  int num = 0;
  int denom = 0;
  int num_cc2p = 0;
  int denom_cc2p = 0;
  int events_2mu = 0; //events with 2 muons
  int events_2other = 0; //events with 2 others
  int n_mom_mu = 0;
  int n_mom_p1 = 0;
  int n_mom_p2 = 0;

  //FV stuff
  float_t FV_edge = 10.0; //How far away we want to be from any TPC edge.  
  float_t xmin = 0.0 + FV_edge; //Just defining the x, y, and z locations of thee FV. 
  float_t xmax = 256.35 - FV_edge;
  float_t ymin = -116.5 + FV_edge;
  float_t ymax = 116.5 - FV_edge;
  float_t zmin = 0.0 + FV_edge;
  float_t zmax = 1036.8 - FV_edge;

  if (fChain == 0) return;
  Long64_t nentries  = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) { //Loop over every entry in the file
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    std::cout<<"-----------------------------------"<<std::endl;
    std::cout<<"BEGINNING TO PROCESS RUN: " <<run << "  SUBRUN: "<< subrun << "  EVENT: " << event <<std::endl;
    std::cout<<"-----------------------------------"<<std::endl;

      if(_debug) std::cout<<"Entry: "<<jentry<<std::endl;
      if(_debug) std::cout<<"-----WHAT KIND OF THINGS ARE IN THE EVENT-----"<<std::endl;
      if(_debug) std::cout<<"[DEBUG] Number of PFParticles in Event:" <<n_pfp_per_event << std::endl;
      if(_debug) std::cout<<"[DEBUG] Number of Tracks in Event:" << n_trk_per_event <<std::endl;
      if(_debug) std::cout<<"[DEBUG] Number of Vertex PFPs:" << vtx_n_pfp <<std::endl;
      if(_debug) std::cout<<"-----------------------------------"<<std::endl;
      
      //First off, we are going to check if the event is within the FV                                                       
      //function returns a boolean for the reconstructed vertex. First number represents how far away from any TPC edge we want to be in cm                                         
      In_FV(10.0, reco_nu_vtxx, reco_nu_vtxy, reco_nu_vtxz);

      ///////////////////////////////////////////////////////////////
      //NOW TO GET STARTED WITH OUR SELECTION
      ////////////////////////////////////////////////////////////////
      std::vector<float> good_trk; //vector used to define good daughters
      std::vector<int> good_trk_id; //id of the good tracks
      std::vector<float> good_trk_length; //vector of the good daughters lengths  
      std::vector<float> good_trk_start_x; //vector of the good tracks start
      std::vector<float> good_trk_end_x; //vector of the good tracks ends
      std::vector<float> good_trk_start_y; //vector of the good tracks start                                                 
      std::vector<float> good_trk_end_y; //vector of the good tracks ends    
      std::vector<float> good_trk_start_z; //vector of the good tracks start                                                 
      std::vector<float> good_trk_end_z; //vector of the good tracks ends    

      //Okay. First thing we ned to do if find events that have a reconstructed vertex within the fiducial volume. 
      //The MicroBooNE FV is is defined as follows: 0 < x < 256: -116 < y < 116: 0 < z < 1056
      ///////////////////////////////////////////////////////
      if(_debug) std::cout<<"-----NOW TO CHECK THAT THE RECONSTRUCTED VERTEX IS IN THE FV-----"<<std::endl;
      if(_debug) std::cout<<"[DEBUG] Location of the Vertex: "<<reco_nu_vtxx<<" "<<reco_nu_vtxy<<" "<<reco_nu_vtxz<<std::endl;
      if(_debug) std::cout<<"-----------------------------------"<<std::endl;

      //Here is where the cut is actually applied
      if((reco_nu_vtxx <= xmin || reco_nu_vtxx >= xmax) || (reco_nu_vtxy <= ymin || reco_nu_vtxy >= ymax) || (reco_nu_vtxz <= zmin || reco_nu_vtxz >= zmax)) continue;
      fvcntr++;

      //Fill the stuff for the denominator:
      ////////////////////////////////////////////
      if(mc_ccnc == 0 && mc_nupdg == 14 && mc_n_proton == 2 && mc_n_pion0 == 0 && mc_n_pionpm == 0 && fv == true){  
	if(_debug) std::cout<<"Denom"<<std::endl;
	h_denom_overlay[0]->Fill(mc_enu,wgt);
	h_denom_overlay[1]->Fill(mc_q2,wgt);
	h_denom_overlay[2]->Fill(mc_X,wgt);
	h_denom_overlay[3]->Fill(mc_Y,wgt);
	for(int i=0; i < mc_g4_pdg->size(); i++){
	  if(mc_g4_pdg->at(i) == 13){
	    for(int j=0; j < mc_g4_mom_muon->size(); j++){
	      h_denom_overlay[4]->Fill(mc_g4_mom_muon->at(j),wgt);
	    }
	  }
	  if(mc_g4_pdg->at(i) == 2212){
	    for(int j=0; j < mc_g4_mom_proton->size(); j++){
	      h_denom_overlay[5]->Fill(mc_g4_mom_proton->at(j),wgt);
	      //h_denom_overlay[6]->Fill(mc_g4_mom_proton->at(i),wgt);
	    }
	  }
	}
	denom++;
	denom_cc2p++;
      }

      //Okay Next: We need to require things to be from the neutrino slice cause otherwise this is going to be a hot mess
      for (int i = 0; i < is_from_nu_slice->size(); i++){

	   if(n_pfp->at(i) == -9999.) continue;
           if(n_trk->at(i) == -9999.) continue;
           if(n_shower->at(i) == -9999.) continue;

	   if(_debug) std::cout<<"-----SANITY CHECK: For Events with Exactly 3 PFP's attached to the Vertex-----" <<std::endl;
	   if(_debug) std::cout<<"[DEBUG] Number of PFParticles in Event: " <<n_pfp_per_event << std::endl;
	   if(_debug) std::cout<<"[DEBUG] Number of Tracks in Event: " << n_trk_per_event <<std::endl;
	   if(_debug) std::cout<<"[DEBUG] Number of Showers in the Event: "<<n_shower_per_event<<std::endl;
	   if(_debug) std::cout<<"[DEBUG] Number of PFParticles at i: "<<n_pfp->at(i)<<std::endl;
	   if(_debug) std::cout<<"[DEBUG] Number of Tracks at i: "<<n_trk->at(i)<<std::endl;
	   if(_debug) std::cout<<"[DEBUG] Number of Showers at i: "<<n_shower->at(i)<<std::endl;
	   if(_debug) std::cout<<"-----------------------------------"<<std::endl;

	   //If a PFP is deemed to come from the neutrino slice we can continue on our merry way                                  
	   ////////////////////////////////////////////////                                                                       
	   //Want to make sure the PFP's are from the nu slice
	   //NOTE: //this is the same thing as if (parentPDG->at(i) != 14 || parentPDG->at(i) !=-14) continue; 
	   if(is_from_nu_slice->at(i) == 0) continue;
	   isfromnucntr++;

	   //this is dumb but it works
	   if(n_pfp->at(i) != 3) continue;
	   has3pfp++;
	   if(n_trk->at(i) != 3) continue;
	   threetrkcntr++;
	   if(n_shower->at(i) != 0) continue;
	   has0shower++;

	   //Final cut! Want to ensure that the 3D location of the start/end of PFP/track is < 5cm from  reco vertex
	   /////////////////////////////////////////////////////
	   //Here i where we define our 3d track start difference  
	   float reco_3d_diff_start = sqrt(pow((reco_nu_vtxx - reco_start_x->at(i)),2) + 
					   pow((reco_nu_vtxy - reco_start_y->at(i)),2) + 
					   pow((reco_nu_vtxz - reco_start_z->at(i)),2)); 

	   float reco_3d_diff_end = sqrt(pow((reco_nu_vtxx - reco_end_x->at(i)),2) +
                                           pow((reco_nu_vtxy - reco_end_y->at(i)),2) +
                                           pow((reco_nu_vtxz - reco_end_z->at(i)),2));

	   if(_debug) std::cout<<"-----------------------------------"<<std::endl;
	   if(_debug) std::cout<<"[DEBUG] The Vertex Resolution From Track Start: "<<reco_3d_diff_start<<std::endl;
	   if(_debug) std::cout<<"[DEBUG] The Vertex Resolution From Track End: "<<reco_3d_diff_end<<std::endl;
	   if(_debug) std::cout<<"-----------------------------------"<<std::endl;

	   //Now we apply the vertex resolution cut
	   if(reco_3d_diff_start < 5.0 || reco_3d_diff_end < 5.0) {
	     good_trk.push_back(float(reco_3d_diff_start));
	     good_trk_id.push_back(int(id_pfp->at(i)));
	     good_trk_length.push_back(float(reco_length->at(i)));
	     good_trk_start_x.push_back(float(reco_start_x->at(i)));
	     good_trk_end_x.push_back(float(reco_end_x->at(i)));			
	     good_trk_start_y.push_back(float(reco_start_y->at(i)));
	     good_trk_end_y.push_back(float(reco_end_y->at(i)));
	     good_trk_start_z.push_back(float(reco_start_z->at(i)));
	     good_trk_end_z.push_back(float(reco_end_z->at(i)));
	   }
	   
	   if(_debug) std::cout<<"-----SANITY CHECK 2: For Events with Exactly 3 PFP's attached to the Vertex-----" <<std::endl;
	   if(_debug) std::cout<<"[DEBUG] Number of PFParticles in Event: " <<n_pfp_per_event << std::endl;
	   if(_debug) std::cout<<"[DEBUG] Number of Tracks in Event: " << n_trk_per_event <<std::endl;
	   if(_debug) std::cout<<"[DEBUG] Number of Showers in the Event: "<<n_shower_per_event<<std::endl;
	   if(_debug) std::cout<<"[DEBUG] Number of PFParticles at i: "<<n_pfp->at(i)<<std::endl;
	   if(_debug) std::cout<<"[DEBUG] Number of Tracks at i: "<<n_trk->at(i)<<std::endl;
	   if(_debug) std::cout<<"[DEBUG] Number of Showers at i: "<<n_shower->at(i)<<std::endl;
	   if(_debug) std::cout<<"-----------------------------------"<<std::endl;

      }// end of if from nu slice

      //Make sure that there are 3 good daughter tracks
      //////////////////////////////////////////////////
      if(_debug) std::cout<<"Number of Good Tracks: : "<<good_trk.size()<<std::endl;
      if(good_trk.size() != 3) continue;
      vectorsize3++;
	 
      //Final Cut in Intial Selection. We have to check containment of the two shortest tracks
      ///////////////////////////////////////////////////////////////////
      int longest_trk_index = max_element(good_trk_length.begin(),good_trk_length.end())-good_trk_length.begin();
      float largest_trk = *max_element(good_trk_length.begin(),good_trk_length.end());
      int shortest_trk_index = min_element(good_trk_length.begin(),good_trk_length.end())-good_trk_length.begin();
      float shortest_trk = *min_element(good_trk_length.begin(),good_trk_length.end());
      
      if(_debug) std::cout<<"The Longest Track's Index: "<<longest_trk_index<<std::endl;
      if(_debug) std::cout<<"The Longest Track in the Vector: "<<largest_trk<<std::endl;
      if(_debug) std::cout<<"The Shortest Track's Index: "<<shortest_trk_index<<std::endl;
      if(_debug) std::cout<<"The Shortest Track in the Vector: "<<shortest_trk<<std::endl;

      //Check containment of the two shortest tracks: Produces a boolean was can use for the actual cut
      Check_Containment(10, shortest_trk_index, longest_trk_index,good_trk,good_trk_start_x,good_trk_start_y,good_trk_start_z,good_trk_end_x,good_trk_end_y,good_trk_end_z);

      //Here is where we make the actual cut
      if(second_trk != true) continue;
      secondtrkgood++;
      if(short_trk != true) continue;
      shortesttrkgood++;

      if(_debug) std::cout<<"[DEBUG] We finished the initial Selection. Yeah! Onto PID"<<std::endl;
	 
      /////////////////////////////
      //Okay......now to apply PID
      /////////////////////////////
      std::vector<float> protons_vector; //vector of the protons
      std::vector<float> muon_vector; //vector of the muon particles 
      std::vector<float> other_vector; //vector of the other particles 
      float muon_id;
      float leading_proton_id;
      float recoil_proton_id;

      for(int v = 0; v < good_trk_id.size(); v++){ //loop through the tracks again
	float value = good_trk_id.at(v);
	if(chi2p_3D->at(value) < 70){ //this is our PID cut. Indicates track is a proton
	  protons_vector.push_back(value);
	}else if(chi2p_3D->at(value) > 70){
	  muon_vector.push_back(value); //id of the other track
	}else{
	  other_vector.push_back(value);
	}
      }

      if(_debug) std::cout<<"Length of protons_vector for this event: "<<protons_vector.size()<<std::endl;
      if(_debug) std::cout<<"Length of muon_vector for this event: "<<muon_vector.size()<<std::endl;
      if(_debug) std::cout<<"Length of other_vector for this event: "<<other_vector.size()<<std::endl; 
      
      if(muon_vector.size() == 2) {
	events_2mu++;
      }
      if(other_vector.size() == 2) {
	events_2other++;
      }

      //PID CUTS
      if(other_vector.size() != 0) continue; //cut out events that have the other vector
      pid_cut0++;
      if(protons_vector.size() + muon_vector.size() != 3) continue; //cut out events that don't have 3 tracks
      pid_cut1++;
      if(protons_vector.size() != 2) continue; //cut out events that don't have 2 protons
      proton_cut++;
      if(muon_vector.size() != 1) continue; //cut out events that don't have 1 muon
      muon_cut++;

      //First Identify which track is the muon, leading proton, and recoil protons
      muon_id = muon_vector[0];  
      float mom0 = reco_mom_proton->at(protons_vector[0]);
      float mom1 = reco_mom_proton->at(protons_vector[1]);
      if (abs(mom0) > abs(mom1)){
	leading_proton_id = protons_vector[0];
	recoil_proton_id = protons_vector[1];
      }else{
	  leading_proton_id = protons_vector[1];
	  recoil_proton_id = protons_vector[0];
      }
      std::cout<<"muon id: "<<muon_id<<std::endl; 
      std::cout<<"leading proton id: "<<leading_proton_id<<std::endl;
      std::cout<<"reocil proton id: "<<recoil_proton_id<<std::endl;

      //remove events in which we can't reconstruct the momentum                                                                                                                                                                                                           
      if(reco_mom_muon->at(muon_id) == -9999.) continue;
      n_mom_mu++;
      if(reco_mom_proton->at(leading_proton_id) == -9999.) continue;
      n_mom_p1++;
      if(reco_mom_proton->at(recoil_proton_id) == -9999.) continue;
      n_mom_p2++;

      //Filling stuff for my efficinecy study
      ////////////////////////////////////
      if(mc_ccnc == 0 && mc_nupdg == 14 && mc_n_proton == 2 && mc_n_pion0 == 0 && mc_n_pionpm == 0 && fv == true) {
	if(_debug) std::cout<<"Num"<<std::endl;
	h_num_overlay[0]->Fill(mc_enu,wgt);
	h_num_overlay[1]->Fill(mc_q2,wgt);
	h_num_overlay[2]->Fill(mc_X,wgt);
	h_num_overlay[3]->Fill(mc_Y,wgt);
	for(int i=0; i < mc_g4_pdg->size(); i++){
	  if(mc_g4_pdg->at(i) == 13){
	    for(int j=0; j < mc_g4_mom_muon->size(); j++){
	      h_num_overlay[4]->Fill(mc_g4_mom_muon->at(j),wgt);
	    }
	  }
	  if(mc_g4_pdg->at(i) == 2212){
	    for(int j=0; j < mc_g4_mom_proton->size(); j++){
	      h_num_overlay[5]->Fill(mc_g4_mom_proton->at(j),wgt);
	      //h_num_overlay[6]->Fill(mc_g4_mom_proton->at(i),wgt);                                                                                                                                                                                             
	    }
	  }
	}
	num++;
	num_cc2p++;
      }
 
	if(_debug) std::cout<<"[DEBUG] Finish Processing Run: "<<run<<" Subrun: "<<subrun<<" Event: "<<event<<std::endl;
	if(_debug) std::cout<<"-----------------------------------"<<std::endl;
	
	good_trk.clear();
	events_remaining++;

   }// End of loop through each event

   std::cout<<"-----MODULE SUMMARY-----"<<std::endl;
   std::cout << "[EFFICIENCY] Initial Number of Events: "<<nentries<<std::endl;
   std::cout << "[EFFICIENCY] Number of Events with Vertex in FV: "<<fvcntr<<std::endl;
   std::cout << "[EFFICIENCY] Number of Events with Reco-Truth of the Vertex Location < 5cm: "<<truthcntr<<std::endl;
   if(_debug) std::cout << "[EFFICIENCY] How Many PFPs are in the Nu Slice?: "<<isfromnucntr<<std::endl;
   if(_debug) std::cout << "[EFFICIENCY] Number of Events with 3 PFPs: "<<has3pfp<<std::endl;
   if(_debug) std::cout << "[EFFICIENCY] Number of Events with 0 Showers: "<<has3pfp/3<<std::endl;
   std::cout << "[EFFICIENCY] Number of Events with 3 Tracks: "<<threetrkcntr/3<<std::endl;
   std::cout << "[EFFICIENCY] Number of Events with the 3 Track's start within 5cm of the Vertex: "<<vectorsize3<<std::endl;
   if(_debug) std::cout << "[EFFICIENCY] Number of Events with the Vector Size Equal to 3: "<<vectorsize3<<std::endl;
   std::cout<<"[ANALYZER] Number of Events with the Second Shortest Track Contained: "<<secondtrkgood<<std::endl;
   std::cout<<"[ANALYZER] Number of Events with the Shortest Track Contained: "<<shortesttrkgood<<std::endl;
   std::cout<<"[ANALYZER] Number of Events with The Other Vector Larger than 0: "<<pid_cut0<<std::endl;
   std::cout<<"[ANALYZER] Number of Events with More than 3 tracks: "<<pid_cut1<<std::endl;
   std::cout<<"[ANALYZER] Number of Events with 2 Protons: "<<proton_cut<<std::endl; 
   std::cout<<"[ANALYZER] Number of Events with 1 Muon: "<<muon_cut<<std::endl;
   //   std::cout<<"[ANALYZER] Number of Events After Momentum Cut: "<<mom_cut<<std::endl;
   std::cout<<"[ANALYZER] Muon Momentum Quality Cut: "<<n_mom_mu<<std::endl;
   std::cout<<"[ANALYZER] Leading Proton Momentum Quality Cut: "<<n_mom_p1<<std::endl;
   std::cout<<"[ANALYZER] Recoil Proton Momentum Quality Cut: "<<n_mom_p2<<std::endl;
   std::cout << "[ANALYZER] Sanity Check of the Total Number of Events Remaining: "<<events_remaining<<std::endl;
   std::cout << "[EFFICIENCY] Number of Generated 3 Prong Events that were Reconstructed: "<<num<<std::endl;
   std::cout << "[EFFICIENCY] Number of Generated 3 Prong Events: "<<denom<<std::endl;
   std::cout << "[EFFICIENCY] Rough Efficiency Estimate: "<<float(100.*(float(num)/float(denom)))<<"%"<<std::endl;
   std::cout << "[EFFICIENCY] Rough Efficiency Estimate For CC2p: "<<float(100.*(float(num_cc2p)/float(denom_cc2p)))<<"%"<<std::endl;
   std::cout << "[PURITY] Rough Purity Estimate For CC2p: "<<float(100.*(float(num_cc2p)/float(events_remaining)))<<"%"<<std::endl;
   std::cout <<"-----CLOSING TIME. YOU DON'T HAVE TO GO HOME, BUT YOU CAN'T STAY HERE-----"<<std::endl;
   
   //Don't forget to write all of your histograms before you leave!
   ///////////////////////////////////////////////////////////////
   tfile->cd();
   for(int i = 0; i < h_list.size(); i++){
     h_list[i]->Write();
   }
   tfile->Close(); //write the root file that contains our histograms                                    

} //End of program
