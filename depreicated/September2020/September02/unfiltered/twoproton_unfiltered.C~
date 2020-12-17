#define twoproton_unfiltered_cxx
#include "twoproton_unfiltered.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>

void twoproton_unfiltered::Loop()
{

  //Making a new Root File that will contain all the histograms that we will want to plot:
  ///////////////////////////////////////////////////////////////////////////////////////
  TFile *tfile = new TFile("histograms_unfiltered.root","RECREATE");

  //Files with RSE's in them                                                                            
  ofstream myfile;//File that will contain RSE of good events                                          
  ofstream cc2p; //File that will contain good cc2p events                                                                                                 
  myfile.open("files_unfiltered.list");
  cc2p.open("files_unfiltered_cc2p.list");
  myfile<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;
  cc2p<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;

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
  int events_remaining = 0; //sanity check for number of events remaining

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

  //shortest track and second shortest track stuff:
  bool second_trk;
  bool short_trk;

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

      /////////////////////////////////////////////////////////////////                                                                                    
      //NOW TO FILL ALL THE HISTOGRAMS THAT ARE BEFORE WE MAKE ALL OUR CUTS 
      //////////////////////////////////////////////////////////////////  
      Fill_Histograms_Mine(0);
      Fill_Histograms_Raquel(0);

      //Neutrinos shit
      ////////////
      if(n_neutrinos_per_event == 0){
	neutrinos_0++;
      }else if(n_neutrinos_per_event == 1){
	neutrinos_1++;
      }else{
	neutrinos_else++;
      }

      ///////////////////////////////////////////////////////////////
      //NOW TO GET STARTED WITH OUR SELECTION
      ////////////////////////////////////////////////////////////////
      std::vector<float> good_trk; //vector used to define good daughters
      std::vector<float> good_trk_length; //vector of the good daughters lengths
      std::vector<int> good_trk_id;
      std::vector<int> good_trk_pdg;
      std::vector<float> good_trk_start_x; //vector of the good tracks start
      std::vector<float> good_trk_end_x; //vector of the good tracks ends
      std::vector<float> good_trk_start_y; //vector of the good tracks start                                                                                                                                                          
      std::vector<float> good_trk_end_y; //vector of the good tracks ends    
      std::vector<float> good_trk_start_z; //vector of the good tracks start                                                                                                                                                          
      std::vector<float> good_trk_end_z; //vector of the good tracks ends    

      //First: Find events that have a reco vtx w/in FV. The MicroBooNE FV is: 0 < x < 256: -116 < y < 116: 0 < z < 1056
      ///////////////////////////////////////////////////////
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
	
	//If a PFP is deemed to come from the neutrino slice we can continue on our merry way                                  
	//NOTE: This is equivalent to: if (parentPDG->at(i) != 14 || parentPDG->at(i) !=-14) continue; 
	///////////////////////////////////////////////////////
	if(is_from_nu_slice->at(i) == 0) continue;
	isfromnucntr++;

	//Make sure there are only 3 pfps, that are tracks, and no showers
	if(n_pfp->at(i) != 3) continue;
	has3pfp++;
	if(n_trk->at(i) != 3) continue;
	threetrkcntr++;
	if(n_shower->at(i) != 0) continue;
	has0shower++;

	//Final cut! Want to ensure that the 3D location of the start/end of PFP/track is < 5cm from  reco vertex
	/////////////////////////////////////////////////////
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
	  good_trk_length.push_back(float(reco_length->at(i)));
	  good_trk_id.push_back(int(id_pfp->at(i))); //grab the id of all the tracks that pass this cut
	  good_trk_pdg.push_back(int(mc_pdg->at(i)));
	  good_trk_start_x.push_back(float(reco_start_x->at(i)));
	  good_trk_end_x.push_back(float(reco_end_x->at(i)));			
	  good_trk_start_y.push_back(float(reco_start_y->at(i)));
	  good_trk_end_y.push_back(float(reco_end_y->at(i)));
	  good_trk_start_z.push_back(float(reco_start_z->at(i)));
	  good_trk_end_z.push_back(float(reco_end_z->at(i)));
	
	  std::cout<<"MC PDG inide the resolution loop: "<<mc_pdg->at(i)<<std::endl;
	}
	   
	if(_debug) std::cout<<"-----SANITY CHECK 2: For Events with Exactly 3 PFP's attached to the Vertex-----" <<std::endl;
	if(_debug) std::cout<<"[DEBUG] Number of PFParticles in Event: " <<n_pfp_per_event << std::endl;
	if(_debug) std::cout<<"[DEBUG] Number of Tracks in Event: " << n_trk_per_event <<std::endl;
	if(_debug) std::cout<<"[DEBUG] Number of Showers in the Event: "<<n_shower_per_event<<std::endl;
	if(_debug) std::cout<<"[DEBUG] Number of PFParticles at i: "<<n_pfp->at(i)<<std::endl;
	if(_debug) std::cout<<"[DEBUG] Number of Tracks at i: "<<n_trk->at(i)<<std::endl;
	if(_debug) std::cout<<"[DEBUG] Number of Showers at i: "<<n_shower->at(i)<<std::endl;
	if(_debug) std::cout<<"-----------------------------------"<<std::endl;

	if(_debug) std::cout<<"[DEBUG] We finish with all the cuts. Yeah!"<<std::endl;

	h_overlay[0]->Fill(n_pfp->at(i), wgt);
	h_overlay[1]->Fill(vtx_n_pfp, wgt);
	h_overlay[2]->Fill(n_trk->at(i), wgt);
	h_overlay[3]->Fill(n_shower->at(i), wgt);
	
	 }// end of if from nu slice
      
      if(_debug) std::cout<<"Number of Good Tracks: "<<good_trk.size()<<std::endl;

      //Make sure that there are 3 good daughter tracks
      /////////////////////////////////////////////////
      if(good_trk.size() != 3) continue;
      vectorsize3++;

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

      for(int j=0; j < 3; j++){

	if(j != shortest_trk_index && j != longest_trk_index){ //second shortest track/second longest track
	  if((good_trk_start_x.at(j) <= xmin || good_trk_start_x.at(j) >= xmax || good_trk_end_x.at(j) <= xmin || good_trk_end_x.at(j) >= xmax) 
	     || (good_trk_start_y.at(j) <= ymin || good_trk_start_y.at(j) >= ymax || good_trk_end_y.at(j) <= ymin || good_trk_end_y.at(j) >= ymax) 
	     || (good_trk_start_z.at(j) <= zmin || good_trk_start_z.at(j) >= zmax || good_trk_end_z.at(j) <= zmin || good_trk_end_z.at(j) >= zmax)) {
	    second_trk = false;
	  } else{
	    second_trk = true;
	  }//end of the else

	}else if( j == shortest_trk_index) {//shortest track
	  if((good_trk_start_x.at(j) <= xmin || good_trk_start_x.at(j) >= xmax || good_trk_end_x.at(j) <= xmin || good_trk_end_x.at(j) >= xmax)
             || (good_trk_start_y.at(j) <= ymin || good_trk_start_y.at(j) >= ymax || good_trk_end_y.at(j) <= ymin || good_trk_end_y.at(j) >= ymax)
             || (good_trk_start_z.at(j) <= zmin || good_trk_start_z.at(j) >= zmax || good_trk_end_z.at(j) <= zmin || good_trk_end_z.at(j) >= zmax)) {
	    short_trk = false;
	  } else{
	    short_trk = true;
	  } //end of the else
	} //end of the shortest track loop
      } //end of the for loop

      //Here is where we make the actual cut
      if(second_trk != true) continue;
      secondtrkgood++;
      if(short_trk != true) continue;
      shortesttrkgood++;

      bool exists = std::find(std::begin(good_trk_pdg),std::end(good_trk_pdg),-9999) != std::end(good_trk_id);
      if(exists != 1){
	std::cout<<"EXXIST: "<<exists<<std::endl;
        for(int i =0; i < good_trk_pdg.size(); i++){
	  std::cout<<"MC_PDG: "<<good_trk_pdg.at(i)<<std::endl;
        }
      }

      /*
       //I am printig out some stuff cause I don't unnderstannd what is happening
      for(int i = 0; i < mc_pdg->size(); i++){
	std::cout<<"MC Pdg at "<<i<<" : "<<mc_pdg->at(i)<<std::endl;
	std::cout<<"MC Start at "<<i<<" : "<<mc_start_x_sce->at(i)<<std::endl;
	std::cout<<"MC End at "<<i<<" : "<<mc_end_x_sce->at(i)<<std::endl;
      }
      for(int i = 0; i < mc_g4_pdg->size(); i++){
	std::cout<<"MC G4 Pdg at "<<i<<" : "<<mc_g4_pdg->at(i)<<std::endl;
	std::cout<<"MC G4 Mom at "<<i<<" : "<<mc_g4_mom_all->at(i)<<std::endl;
	std::cout<<"MC G4 Start at "<<i<<" : "<<mc_g4_start_x_sce->at(i)<<std::endl;
	std::cout<<"MC G4 End at "<<i<<" : "<<mc_g4_end_x_sce->at(i)<<std::endl;
      }
      for(int i = 0; i < parentPDG->size(); i++){
	std::cout<<"Parent  Pdg at "<<i<<" : "<<parentPDG->at(i)<<std::endl;
      }
      for(int i = 0; i < good_trk_id.size(); i++){
	std::cout<<"Good Track ID at "<<i<<" : "<<good_trk_id.at(i)<<std::endl;
      }
      for(int i = 0; i < mc_mom->size(); i++){
	std::cout<<"MC Momentum at "<<i<<" : "<<mc_mom->at(i)<<std::endl;
      }
      for(int i = 0; i < reco_mom->size(); i++){
	std::cout<<"Reco Momentum at "<<i<<" : "<<reco_mom->at(i)<<std::endl;
	std::cout<<"Reco Muon Momentum at "<<i<<" : "<<reco_mom_muon->at(i)<<std::endl;
	std::cout<<"Reco Proton Momentum at "<<i<<" : "<<reco_mom_proton->at(i)<<std::endl;
      }
      

      //how many threshold particles are there?
      std::cout<<"Number of Threshold Muons: "<<mc_n_threshold_muon<<std::endl;
      std::cout<<"Number of Threshold Protons: "<<mc_n_threshold_proton<<std::endl;
      std::cout<<"Number of Threshold Pi+-: "<<mc_n_threshold_pionpm<<std::endl;
      std::cout<<"Number of Threshold Pi0: "<<mc_n_threshold_pion0<<std::endl;
      std::cout<<"Number of Threshold Electron: "<<mc_n_threshold_electron<<std::endl;
      std::cout<<"Number of Threshold Neutrons: "<<mc_n_threshold_neutron<<std::endl;

      std::cout<<"Length_ID: "<<good_trk_id.size()<<std::endl;
      std::cout<<"Length of PDG G4: "<<mc_g4_pdg->size()<<std::endl;
      std::cout<<"Length of PDG MC: "<<mc_pdg->size()<<std::endl;
      std::cout<<"MC Mode: "<<mc_mode<<std::endl;
      std::cout<<"CCNC: "<<mc_ccnc<<std::endl;
      */

      //My sad attempt to try and get the PID stuff working
      //////////////////////////////////////////////////////
      float chi2p;
      float chi2mu;
      for(int v = 0; v < good_trk_id.size(); v++){

	float value = good_trk_id.at(v);
	std::cout<<"Value of Value: "<<value<<std::endl;
	std::cout<<"value of good_trk_id: "<<good_trk_id.at(v)<<std::endl;

	for(int y = 0; y < 3; y++){

	  if(y == 0){
	    chi2p = chi2_p_0->at(value);
	    chi2mu = chi2_mu_0->at(value);
	  }else if (y == 1){
	    chi2p = chi2_p_1->at(value);
	    chi2mu = chi2_mu_1->at(value);
	  }else if (y == 2){
	    chi2p = chi2_p_2->at(value);
	    chi2mu = chi2_mu_2->at(value);
	  }

	    h_chi2p[y][0]->Fill(chi2p,wgt);
	    h_chi2mu[y][0]->Fill(chi2mu,wgt);

	    //if(mc_pdg->at(value) == -9999.) continue;
	    //if(chi2p == -9999.) continue;
	    //if(chi2mu == -9999.) continue;
	
	    if(mc_pdg->at(value) == 2212 || mc_pdg->at(value) == -2212){
	      h_chi2p[y][1]->Fill(chi2p,wgt);
	      h_chi2mu[y][1]->Fill(chi2mu,wgt);

	    } else if(mc_pdg->at(value) == 211 || mc_pdg->at(value) == -211 || mc_pdg->at(value) == 111) {
	      h_chi2p[y][3]->Fill(chi2p,wgt);
	      h_chi2mu[y][3]->Fill(chi2mu,wgt);

	    } else if(mc_pdg->at(value) == 13 || mc_pdg->at(value) == -13){
	      h_chi2p[y][2]->Fill(chi2p,wgt);
	      h_chi2mu[y][2]->Fill(chi2mu,wgt);

	    } else if(mc_pdg->at(value) == 11 || mc_pdg->at(value) == -11){
	      h_chi2p[y][4]->Fill(chi2p,wgt);
	      h_chi2mu[y][4]->Fill(chi2mu,wgt);
        
	    } else {
	      h_chi2p[y][5]->Fill(chi2p,wgt);
	      h_chi2mu[y][5]->Fill(chi2mu,wgt);
	    }
	
	} //end of loop over the plane
      } //end of loop over the particle ids

      /////////////////////////////////////////////////////////////////
      //NOW TO FILL ALL THE HISTOGRAMS THAT ARE AFTER WE MADE ALL OUR CUTS
      //////////////////////////////////////////////////////////////////
      h_correlation_overlay[0]->Fill(reco_nu_vtxx,reco_nu_vtxy,wgt);
      h_correlation_overlay[1]->Fill(mc_nu_vtxx,mc_nu_vtxy,wgt);
      h_correlation_overlay[2]->Fill(mc_nu_vtxx_sce,mc_nu_vtxy_sce,wgt);

      Fill_Histograms_Mine(1);
      Fill_Histograms_Raquel(1);

      if(_debug) std::cout<<"[DEBUG] Finish Processing Run: "<<run<<" Subrun: "<<subrun<<" Event: "<<event<<std::endl;
      if(_debug) std::cout<<"-----------------------------------"<<std::endl;

      //One last thing: Make sure to ssave the RSE for the good events before you end your loop                                                     
      myfile << run << " " << subrun << " " << event << " " ;
      myfile << endl;

      good_trk.clear();
      events_remaining++;

   }// End of loop through each event


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
   std::cout << "[MC_RECO_RAQUEL] Number of CCQEL Events: "<<qel[1]<<" Fraction of the Total: "<<float(100.*(float(qel[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO_RAQUEL] Number of CCRES Events: "<<res[1]<<" Fraction of the Total: "<<float(100.*(float(res[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO_RAQUEL] Number of CCMEC Events: "<<mec[1]<<" Fraction of the Total: "<<float(100.*(float(mec[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO_RAQUEL] Number of CCCOH Events: "<<coh[1]<<" Fraction of the Total: "<<float(100.*(float(coh[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO_RAQUEL] Number of CCDIS Events: "<<dis[1]<<" Fraction of the Total: "<<float(100.*(float(dis[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO_RAQUEL] Number of CCNue Events: "<<ccnue_raquel[1]<<" Fraction of the Total: "<<float(100.*(float(ccnue_raquel[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO_RAQUEL] Number of OUTFV Events: "<<outfv_raquel[1]<<" Fraction of the Total: "<<float(100.*(float(outfv_raquel[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO_RAQUEL] Number of NC Events: "<<nc_raquel[1]<<" Fraction of the Total: "<<float(100.*(float(nc_raquel[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO_RAQUEL] Number of Else Events: "<<other_raquel[1]<<" Fraction of the Total: "<<float(100.*(float(other_raquel[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout <<"-----ONE FOR THE DAGGER, AND ONE FOR THE ONE YOU BELIEVE!!-----"<<std::endl;

   std::cout<<"-----MC RECO'D SUMMARY-----"<<std::endl;
   std::cout << "[MC_RECO] Initial Number of Events That were Reconstructed: "<<events_remaining<<std::endl;
   std::cout << "[MC_RECO] Number of CCOpOpi Events: "<<cc0p0pi[1]<<" Fraction of the Total: "<<float(100.*(float(cc0p0pi[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO] Number of CC1p0pi Events: "<<cc1p0pi[1]<<" Fraction of the Total: "<<float(100.*(float(cc1p0pi[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO] Number of CC2p0pi Events: "<<cc2p0pi[1]<<" Fraction of the Total: "<<float(100.*(float(cc2p0pi[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO] Number of CCNp0pi Events: "<<ccNp0pi[1]<<" Fraction of the Total: "<<float(100.*(float(ccNp0pi[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO] Number of CCNp1pi Events: "<<ccNp1pi[1]<<" Fraction of the Total: "<<float(100.*(float(ccNp1pi[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO] Number of CCNpNpi Events: "<<ccNpNpi[1]<<" Fraction of the Total: "<<float(100.*(float(ccNpNpi[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO] Number of CCNue Events: "<<ccnue[1]<<" Fraction of the Total: "<<float(100.*(float(ccnue[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO] Number of OUTFV Events: "<<outfv[1]<<" Fraction of the Total: "<<float(100.*(float(outfv[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO] Number of NC Events: "<<nc[1]<<" Fraction of the Total: "<<float(100.*(float(nc[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout << "[MC_RECO] Number of Else Events: "<<other[1]<<" Fraction of the Total: "<<float(100.*(float(other[1])/float(events_remaining)))<<"%"<<std::endl;
   std::cout <<"-----NOTHING REALLY MATTERS. ANYONE CAN SEE. NOTHING REALLY MATTERS. NOTHING REALLY MATTERS TO ME-----"<<std::endl;

   std::cout<<"Neutrinos 0: "<<neutrinos_0<<std::endl;
   std::cout<<"Neutrinos 1: "<<neutrinos_1<<std::endl;
   std::cout<<"Neutrinos Else: "<<neutrinos_else<<std::endl;

   //Don't forget to write all of your histograms before you leave!
   ///////////////////////////////////////////////////////////////
   tfile->cd();
   Write_Histograms(); //function that writes all our histograms
   tfile->Close(); //write the root file that contains our histograms                                                                                 
   myfile.close(); //Write the file that contains the RSE of good events 
   cc2p.close(); //Write the file that contains the RSE of good 1mu2p events

} //End of program

