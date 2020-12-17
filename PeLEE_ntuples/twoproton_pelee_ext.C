#define twoproton_pelee_ext_cxx
#include "twoproton_pelee_ext.h"
#include "histogram_funcs.h"

void twoproton_pelee_ext::Loop()
{
  //Define objects of classes
  ////////////////////////////
  histogram_funcs hist; //histogram_funcs.h
  
  //Making a new Root File that will contain all the histograms that we will want to plot:                                    
  ///////////////////////////////////////////////////////////////////////////////////////                                      
  TFile *tfile = new TFile("root_files/histograms_filtered_ext.root","RECREATE");

  //File with RSE's in them                                                                                                   
  ofstream myfile;//File that will contain RSE of good events                                                                 
  myfile.open("lists/files_filtered_ext.list");
  myfile<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;

  //Define all the histograms I am going to fill
 ////////////////////////////////////////////                                                                                
  hist.Define_Histograms("ext");

  //Defining all the constans we will use later                                                                                          
  //////////////////////////////                                                                                                      
  bool _debug = false; //debug statements                                                                         
  double pot_wgt = 0.2993; //POT weight                                    
  
  
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
      std::cout<<"BEGINNING TO PROCESS RUN: " <<run << "  SUBRUN: "<< sub << "  EVENT: " << evt <<std::endl;
      std::cout<<"-----------------------------------"<<std::endl;

      if(_debug) std::cout<<"Entry: "<<jentry<<std::endl;

      //Okay. The CCInclusive Selection requires the following things:
      // 1) The reconstructed neutrino vertex is inside the FV 
      // 2) The start/end point of every pfp is within the FV (i.e. contained)
      // 3) The topological score of every pfp is above 0.1 (filled per event): topological_score > 0.1
      // 4) The cosmic impact parameter is greater than 10 cm for each PFP (filled per event): CosmicIP > 10 
      // We are now goinnng to make plots of those cut variables

      hist.h_topological_score->Fill(topological_score,pot_wgt);
      hist.h_cosmic_impact_parameter->Fill(CosmicIP,pot_wgt);



      /* if(_debug) std::cout<<"-----WHAT KIND OF THINGS ARE IN THE EVENT-----"<<std::endl;
      if(_debug) std::cout<<"[DEBUG] Number of PFParticles in Event:" <<n_pfp_per_event << std::endl;
      if(_debug) std::cout<<"[DEBUG] Number of Tracks in Event:" << n_trk_per_event <<std::endl;
      if(_debug) std::cout<<"[DEBUG] Number of Vertex PFPs:" << vtx_n_pfp <<std::endl;
      if(_debug) std::cout<<"-----------------------------------"<<std::endl;
      */ 
   } //end of Loop over events

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
   hist.Write_Histograms(false); //function that writes all our histograms                                                              
   tfile->Close(); //write the root file that contains our histograms                                                         
   myfile.close(); //Write the file that contains the RSE of good events                                                     

} //end of progrm
