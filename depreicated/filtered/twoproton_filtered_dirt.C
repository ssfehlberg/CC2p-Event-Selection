#define twoproton_filtered_dirt_cxx
#include "twoproton_filtered_dirt.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void twoproton_filtered_dirt::Loop()
{
  //Making a new Root File that will contain all the histograms that we will want to plot:                                                                ///////////////////////////////////////////////////////////////////////////////////////                                                              
  TFile *tfile = new TFile("histograms_filtered_dirt.root","RECREATE");

  //File with RSE's in them                                                                                                                           
  ofstream myfile;//File that will contain RSE of good events                                                                                          
  myfile.open("files_filtered_dirt.list");
  myfile<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;

  //Define all the histograms I am going to fill                                                                                                       
  ////////////////////////////////////////////                                                                                                         
  Define_Histograms();

  //Defining all the constans we will use later                                                                                                        
  //////////////////////////////                                                                                                                       
  bool _debug = false; //debug statements                                                                                                                
  //Counters                                                                                                                           
  int fvcntr = 0; //Number of events with reconstructed vertex within the FV                                                          
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
  int events_remaining = 0; //sanity check for number of events remaining                                                             
  int events_chi2p = 0;
  int events_chi2mu = 0;
  int events_2mu = 0; //events with 2 muons
  int events_2other = 0; //events with 2 others
  int n_mom_mu = 0;
  int n_mom_p1 = 0;
  int n_mom_p2 = 0;

  //neutrino counters                                                                                                                                  
  int neutrinos_0 = 0;
  int neutrinos_1 = 0;
  int neutrinos_else = 0;

  //FV Stuff                                                                                                                                           
  float_t FV_edge = 10.0;
  float_t xmin = 0.0 + FV_edge;
  float_t xmax = 256.35 - FV_edge;
  float_t ymin = -116.5 + FV_edge;
  float_t ymax = 116.5 - FV_edge;
  float_t zmin = 0.0 + FV_edge;
  float_t zmax = 1036.8 - FV_edge;

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
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
      
      /////////////////////////////////////////////////////////////////                                                                                
      //NOW TO FILL ALL THE HISTOGRAMS THAT ARE BEFORE WE MAKE ALL OUR CUTS                                                                                   //////////////////////////////////////////////////////////////////                                                                               
      Fill_Histograms(0);

      //Neutrinos shit                                                                                                                                        ////////////                                                                                                                                     
      if(n_neutrinos_per_event == 0){
        neutrinos_0++;
      }else if(n_neutrinos_per_event == 1){
        neutrinos_1++;
      }else{
        neutrinos_else++;
      }

      ///////////////////////////////////////////////////////////////                                                                                         //NOW TO GET STARTED WITH OUR SELECTION                                                                                                                 ////////////////////////////////////////////////////////////////                                                                                 
      std::vector<float> good_trk; //vector used to define good daughters                                        
      std::vector<int> good_trk_id;                                 
      std::vector<float> good_trk_length; //vector of the good daughters lengths
      std::vector<float> good_trk_start_x; //vector of the good tracks start
      std::vector<float> good_trk_end_x; //vector of the good tracks ends
      std::vector<float> good_trk_start_y; //vector of the good tracks start                                                                                                                                                                                                            
      std::vector<float> good_trk_end_y; //vector of the good tracks ends    
      std::vector<float> good_trk_start_z; //vector of the good tracks start                                                                                                                                                                                                            
      std::vector<float> good_trk_end_z; //vector of the good tracks ends

      //First: Find events that have a reco vtx w/in FV. The MicroBooNE FV is: 0 < x < 256: -116 < y < 116: 0 < z < 1056                                      ///////////////////////////////////////////////////////                                                                                          
      if(_debug) std::cout<<"-----NOW TO CHECK THAT THE RECONSTRUCTED VERTEX IS IN THE FV-----"<<std::endl;
      if(_debug) std::cout<<"[DEBUG] Location of the Vertex: "<<reco_nu_vtxx<<" "<<reco_nu_vtxy<<" "<<reco_nu_vtxz<<std::endl;
      if(_debug) std::cout<<"-----------------------------------"<<std::endl;

      //Here is where the cut is actually applied                                                                                                      
      if((reco_nu_vtxx <= xmin || reco_nu_vtxx >= xmax) || (reco_nu_vtxy <= ymin || reco_nu_vtxy >= ymax) || (reco_nu_vtxz <= zmin || reco_nu_vtxz >= zmax)) continue;
      fvcntr++;

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

        //If a PFP is deemed to come from the neutrino slice we can continue on our merry way                                                                   //NOTE: This is equivalent to: if (parentPDG->at(i) != 14 || parentPDG->at(i) !=-14) continue;                                                          ///////////////////////////////////////////////////////                                                                                        
	if(is_from_nu_slice->at(i) == 0) continue;
        isfromnucntr++;

	//Make sure there are only 3 pfps, that are tracks, and no showers                                                                              
	if(n_pfp->at(i) != 3) continue;
        has3pfp++;
        if(n_trk->at(i) != 3) continue;
        threetrkcntr++;
        if(n_shower->at(i) != 0) continue;
        has0shower++;

	//Final cut! Want to ensure that the 3D location of the start/end of PFP/track is < 5cm from  reco vertex                                               /////////////////////////////////////////////////////                                                                                          
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

	h_overlay[0]->Fill(n_pfp->at(i), wgt);
        h_overlay[1]->Fill(vtx_n_pfp, wgt);
        h_overlay[2]->Fill(n_trk->at(i), wgt);
        h_overlay[3]->Fill(n_shower->at(i), wgt);

      } //end of if from nu slice
      
      //Make sure that there are 3 good daughter tracks                                                                                   
      if(_debug) std::cout<<"Number of Good Tracks: "<<good_trk.size()<<std::endl;      
      if(good_trk.size() != 3) continue;
      vectorsize3++;
      
      //Fill Hitograms after 3 good tracks
      Fill_Histograms(1);

      //Final Cut. We have to check containment of the two shortest tracks
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

      //Fill Histograms after containment cut since it is the last cut in intial selection:
      ///////////////////////////////////////////////////////////
      Fill_Histograms(2);

      //make sure to grab 3 plane variables
      float chi2p;
      float chi2mu;
      float dEdx;

      //Filling variables that are 3 plane
      for(int i = 0; i < good_trk_id.size(); i ++){
	int value = good_trk_id.at(i);
	for(int y = 0; y < 3; y++){  
	 if(y == 0){
	    chi2p = chi2_p_0->at(value);
	    chi2mu = chi2_mu_0->at(value);
	    dEdx = total_dedx_0->at(value);
	  }else if (y == 1){
	    chi2p = chi2_p_1->at(value);
	    chi2mu = chi2_mu_1->at(value);
	    dEdx = total_dedx_1->at(value);
	  }else if (y == 2){
	    chi2p = chi2_p_2->at(value);
	    chi2mu = chi2_mu_2->at(value);
	    dEdx = total_dedx_2->at(value);
	 }
	 //if(chi2_p_2->at(value) == -9999.) continue;
	 events_chi2p++;
	 //if(chi2_mu_2->at(value) == -9999.) continue;	  
	 events_chi2mu++;
	 h_chi2p[y]->Fill(chi2p,wgt); 
	 h_chi2mu[y]->Fill(chi2mu,wgt); 
	 h_dEdx_total[y]->Fill(dEdx,wgt);
	} //end of loop over the plane
	//if(chi2p_3D->at(value) == -9999.) continue;
        //if(chi2mu_3D->at(value) == -9999.) continue;
	//if(chi2pi_3D->at(value) == -9999.) continue;
	//h_chi2p_3D[0]->Fill(chi2p_3D->at(value),wgt);
        //h_chi2mu_3D[0]->Fill(chi2mu_3D->at(value),wgt);
	//h_chi2pi_3D[0]->Fill(chi2pi_3D->at(value),wgt);
      } //end of loop over particle ids		    

      h_chi2p_3D[0]->Fill(chi2p_3D->at(good_trk_id[longest_trk_index]),wgt);
      h_chi2mu_3D[0]->Fill(chi2mu_3D->at(good_trk_id[longest_trk_index]),wgt);
      h_chi2pi_3D[0]->Fill(chi2pi_3D->at(good_trk_id[longest_trk_index]),wgt);

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
	//std::cout<<"value"<<value<<std::endl;
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

      //////////////////////////////////////
      //Filling histograms after PID cuts:
      /////////////////////////////////////
      Fill_Histograms(3);
      h_chi2p_3D[1]->Fill(muon_vector[0],wgt);
      h_chi2mu_3D[1]->Fill(muon_vector[0],wgt);
      h_chi2pi_3D[1]->Fill(muon_vector[0],wgt);

      //////////////////////////
      //Particle Specific Plots
      ////////////////////////
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

      //remove events in which we can't reconstruct the momentum                                                                                                                                                                                                       
      if(reco_mom_muon->at(muon_id) == -9999.) continue;
      n_mom_mu++;
      if(reco_mom_proton->at(leading_proton_id) == -9999.) continue;
      n_mom_p1++;
      if(reco_mom_proton->at(recoil_proton_id) == -9999.) continue;
      n_mom_p2++;

      Fill_Particles(muon_id,leading_proton_id,recoil_proton_id);

      if(_debug) std::cout<<"[DEBUG] Finish Processing Run: "<<run<<" Subrun: "<<subrun<<" Event: "<<event<<std::endl;
      if(_debug) std::cout<<"-----------------------------------"<<std::endl;

      //One last thing: Make sure to ssave the RSE for the good events before you end your loop                                 
      myfile << run << " " << subrun << " " << event << " " ;
      myfile << endl;

      good_trk.clear();
      events_remaining++;

   } //ends loop over all events

   std::cout<<"-----MODULE SUMMARY-----"<<std::endl;
   std::cout << "[ANALYZER] Initial Number of Events: "<<nentries<<std::endl;
   std::cout << "[ANALYZER] Number of Events with Vertex in FV: "<<fvcntr<<std::endl;
   if(_debug) std::cout << "[ANALYZER] How Many PFPs are in the Nu Slice?: "<<isfromnucntr<<std::endl;
   if(_debug) std::cout << "[ANALYZER] Number of Events with 3 PFPs: "<<has3pfp<<std::endl;
   if(_debug) std::cout << "[ANALYZER] Number of Events with 0 Showers: "<<has3pfp/3<<std::endl;
   std::cout << "[ANALYZER] Number of Events with 3 Tracks: "<<threetrkcntr/3<<std::endl;
   std::cout << "[ANALYZER] Number of Events with the 3 Track's start within 5cm of the Vertex: "<<vectorsize3<<std::endl;
   if(_debug) std::cout << "[ANALYZER] Number of Events with the Vector Size Equal to 3: "<<vectorsize3<<std::endl;
   std::cout<<"[ANALYZER] Number of Events with the Second Shortest Track Contained: "<<secondtrkgood<<std::endl;
   std::cout<<"[ANALYZER] Number of Events with the Shortest Track Contained: "<<shortesttrkgood<<std::endl;
   std::cout<<"[ANALYZER] Number of Events with The Other Vector Larger than 0: "<<pid_cut0<<std::endl;
   std::cout<<"[ANALYZER] Number of Events with More than 3 tracks: "<<pid_cut1<<std::endl;
   std::cout<<"[ANALYZER] Number of Events with 2 Protons: "<<proton_cut<<std::endl; 
   std::cout<<"[ANALYZER] Number of Events with 1 Muon: "<<muon_cut<<std::endl;
   std::cout<<"[ANALYZER] Muon Momentum Quality Cut: "<<n_mom_mu<<std::endl;
   std::cout<<"[ANALYZER] Leading Proton Momentum Quality Cut: "<<n_mom_p1<<std::endl;
   std::cout<<"[ANALYZER] Recoil Proton Momentum Quality Cut: "<<n_mom_p2<<std::endl;
   std::cout << "[ANALYZER] Sanity Check of the Total Number of Events Remaining: "<<events_remaining<<std::endl;
   std::cout <<"-----CLOSING TIME. YOU DON'T HAVE TO GO HOME, BUT YOU CAN'T STAY HERE-----"<<std::endl;

   std::cout<<"Neutrinos 0: "<<neutrinos_0<<std::endl;
   std::cout<<"Neutrinos 1: "<<neutrinos_1<<std::endl;
   std::cout<<"Neutrinos Else: "<<neutrinos_else<<std::endl;
   std::cout<<"events_chi2p++: "<<events_chi2p<<std::endl;
   std::cout<<"events_chi2mu++: "<<events_chi2mu<<std::endl;
   std::cout<<"events_2mu: "<<events_2mu<<std::endl; 
   std::cout<<"events_2other: "<<events_2other<<std::endl;

   //Don't forget to write all of your histograms before you leave!                                                                                    
   ///////////////////////////////////////////////////////////////                                                                                     
   tfile->cd();
   Write_Histograms(); //function that writes all our histograms                                                                                       
   tfile->Close(); //write the root file that contains our histograms                                                                                  
   myfile.close(); //Write the file that contains the RSE of good events                                                                                

} //ends program
