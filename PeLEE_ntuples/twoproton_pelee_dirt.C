#define twoproton_pelee_dirt_cxx
#include "twoproton_pelee_dirt.h"
#include "histogram_funcs.h"
#include "helper_funcs.h"

void twoproton_pelee_dirt::Loop()
{

  //Define objects of classes
  ////////////////////////////
  histogram_funcs hist; //histogram_funcs.h
  helper_funcs cuts; //helper_funcs.h  
  
  //Making a new Root File that will contain all the histograms that we will want to plot:                                    
  ///////////////////////////////////////////////////////////////////////////////////////                                      
  TFile *tfile = new TFile("root_files/histograms_filtered_dirt_wgt.root","RECREATE"); //wgt indicates using the CV MC values

  //File with RSE's in them                                                                                                   
  ofstream myfile;//File that will contain RSE of good events                                                                 
  myfile.open("lists/files_filtered_dirt_wgt.list");
  myfile<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;

  //Define all the histograms I am going to fill
 ////////////////////////////////////////////                                                                                
  hist.Define_Histograms("dirt_wgt");

  //Defining all the constans we will use later                                                                                          
  //////////////////////////////                                                                                                      
  bool _debug = false; //debug statements                                                                        
  double pot_wgt = 0.141; //POT weight
  double mc_wgt; //spline times tuned cv weight
                     
  //Counters                                                                                                       
  int fvcntr = 0; //Number of events with reconstructed vertex within the FV                                      
  int pfp_starts_contained = 0; //how many pfp's start within the FV                                      
  int toposcore = 0; //how many events have topological score above 0.1                                           
  int cosmicip = 0;//how many events with cosmic IP above 10cm                                                    
                   
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

      //Defining the MC Weight cause it is dumb
      /////////////////////////////
      if(std::isfinite(weightTune) && weightTune <= 100.) {
	mc_wgt = weightSplineTimesTune;
      } else {
	mc_wgt = 1 * weightSpline;
      }

    //Filling histograms before any selection is made
    ////////////////////////////////////////////////
    hist.Fill_Histograms(0, TVector3(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z), CosmicIP, topological_score,pot_wgt*mc_wgt);
   
    //Just casually checking how many neutrino slices we have
    if(nslice == 0){
      neutrinos_0++;
    }else if(nslice == 1){
      neutrinos_1++;
    }else{
      neutrinos_else++;
    }

    //Okay. The CCInclusive Selection requires the following things:
    // 1) The reconstructed neutrino vertex is inside the FV 
    // 2) The start point of every pfp is within the FV (i.e. contained)
    // 3) The topological score of every pfp is above 0.1 (filled per event): topological_score > 0.1
    // 4) The cosmic impact parameter is greater than 10 cm for each PFP (filled per event): CosmicIP > 10 
    // We are now goinnng to make plots of those cut variables


    //1) Check that the event is in the FV
    //////////////////////////////////////
    if(cuts.In_FV(10,10,10,10,10,50,reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z) == false) continue; //10 cm border except from backend of detector
    fvcntr++;

    //Fill Histograms
    hist.Fill_Histograms(1, TVector3(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z), CosmicIP, topological_score,pot_wgt*mc_wgt);

    //2) The start point of every pfp is within the FV
    ///////////////////////////////////////////////////
    for(int p = 0; p < n_pfps; p++){
      unsigned int generation = pfp_generation_v->at( p ); //only check direct neutrino daughters (gen ==2)
      if ( generation != 2u ) continue;
      float tscore = trk_score_v->at( p ); //get the track score for the current PFParticle
      // A PFParticle is considered a track if its score is above the track score
      // cut. Get the track or shower start coordinates as appropriate.
      float x, y, z;
      if ( tscore > cuts.TRACK_SCORE_CUT ) {
	x = trk_sce_start_x_v->at( p );
	y = trk_sce_start_y_v->at( p );
	z = trk_sce_start_z_v->at( p );
      }
      else {
	x = shr_start_x_v->at( p );
	y = shr_start_y_v->at( p );
	z = shr_start_z_v->at( p );
      }
      if(cuts.In_FV(10,10,10,10,10,50,x,y,z) == false) continue;
    } //end of Loop through PFP's
    pfp_starts_contained++; 

    //Fill Histograms
    hist.Fill_Histograms(2,TVector3(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z), CosmicIP, topological_score,pot_wgt*mc_wgt);

    //3) The topoloogical score of every neutrino slice is above 0.1
    ///////////////////////////////////////////////////////
    if(topological_score < cuts.TOPO_SCORE_CUT) continue;
    toposcore++;

    //Fill Histograms  
    hist.Fill_Histograms(3, TVector3(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z), CosmicIP, topological_score,pot_wgt*mc_wgt);

    //4) The cosmic impact parameter is greater than 10 cm for every neutrino slice. Honestly a dumb cut. Will remove later
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(CosmicIP < cuts.COSMIC_IP_CUT) continue;
    cosmicip++;

    //Fill Histograms  
    hist.Fill_Histograms(4, TVector3(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z), CosmicIP, topological_score,pot_wgt*mc_wgt);

    //Now to apply the ve and NC rejection cuts. These are slightly modified to match our 1mu2p needs 
    /////////////////////////////////////////////////////////////////////////////////////////////////

    //Next: The Muon Selection
    //////////////////////////

   } //end of Loop over events

   std::cout<<"-----MODULE SUMMARY-----"<<std::endl;
   std::cout << "[ANALYZER] Initial Number of Events: "<<nentries<<std::endl;
   std::cout << "[ANALYZER] Number of Events with Vertex in FV: "<<fvcntr<<std::endl;
   std::cout << "[ANALYZER] How Many Events with All PFP Starts within the FV: "<<pfp_starts_contained<<std::endl;
   std::cout << "[ANALYZER] How Many Events with Topological Score above 0.1: "<<toposcore<<std::endl;
   std::cout << "[ANALYZER] Number of Events with CosmicIP above 10cm: "<<cosmicip<<std::endl;
   std::cout << "[ANALYZER] Number of Events with 3 Tracks: "<<threetrkcntr/3<<std::endl;
   std::cout<<  "[ANALYZER] Number of Events with the Second Shortest Track Contained: "<<secondtrkgood<<std::endl;
   std::cout<<  "[ANALYZER] Number of Events with the Shortest Track Contained: "<<shortesttrkgood<<std::endl;
   std::cout<<  "[ANALYZER] Number of Events with The Other Vector Larger than 0: "<<pid_cut0<<std::endl;
   std::cout<<  "[ANALYZER] Number of Events with More than 3 tracks: "<<pid_cut1<<std::endl;
   std::cout<<  "[ANALYZER] Number of Events with 2 Protons: "<<proton_cut<<std::endl; 
   std::cout<<  "[ANALYZER] Number of Events with 1 Muon: "<<muon_cut<<std::endl;
   std::cout<<  "[ANALYZER] Muon Momentum Quality Cut: "<<n_mom_mu<<std::endl;
   std::cout<<  "[ANALYZER] Leading Proton Momentum Quality Cut: "<<n_mom_p1<<std::endl;
   std::cout<<  "[ANALYZER] Recoil Proton Momentum Quality Cut: "<<n_mom_p2<<std::endl;
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
