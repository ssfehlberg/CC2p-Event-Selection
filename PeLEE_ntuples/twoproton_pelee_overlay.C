#define twoproton_pelee_overlay_cxx
#include "twoproton_pelee_overlay.h"

void twoproton_pelee_overlay::Loop()
{
  //Define objects of classes
  ////////////////////////////
  histogram_funcs hist; //histogram_funcs.h
  helper_funcs cuts; //helper_funcs.h   

  //Making a new Root File that will contain all the histograms that we will want to plot:
  ///////////////////////////////////////////////////////////////////////////////////////
  TFile *tfile = new TFile("root_files/histograms_filtered_wgt.root","RECREATE"); //wgt indicates applying cenntral value MC value

  //Files with RSE's in them                                                                            
  ofstream myfile;//File that will contain RSE of good events                                          
  ofstream cc2p; //File that will contain good cc2p events                                                                          
  //ofstream ccNp0pi_file;     
  myfile.open("lists/files_filtered_wgt.list");
  cc2p.open("lists/files_filtered_wgt_cc2p.list");
  //ccNp0pi_file.open("lists/files_filtered_wgt_ccNp0pi.list"); 
  myfile<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;
  cc2p<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;
  //ccNp0pi_file<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;  

  //Define all the histograms I am going to fill                                
  ////////////////////////////////////////////
  Define_Histograms();

  //Defining all the constans we will use later
  //////////////////////////////
  bool _debug = false; //debug statements
  double pot_wgt = 0.0347; //POT weight
  double mc_wgt;
  double TRACK_SCORE_CUT = 0.5;
  double TOPO_SCORE_CUT = 0.3;
  double COSMIC_IP_CUT = 10.0;
  double MUON_MOM_CUT = 0.1;
  double PROTON_MOM_CUT = 0.25;
  double CHARGED_PI_MOM_CUT = 0.065;
  double PION0_MOM_CUT = 0.065;
  double MASS_PROTON = 0.93827208;
  double MASS_MUON = 0.10565837;
  double MASS_PION0 = 0.13497666;
  double MASS_PIONPM =0.13957000;
 
  //Counters
  int fvcntr = 0; //Number of events with reconstructed vertex within the FV                                      
  //int pfp_starts_contained = 0; //how many pfp's start within the FV                                      
  int toposcore = 0; //how many events have topological score above 0.1                                           
  //int cosmicip = 0;//how many events with cosmic IP above 10cm                                                    

  int isfromnucntr = 0; //how many pfp's are from the neutrino slice
  int has3pfp = 0; //how many events have exactly 3 pfps
  int has0shower = 0;//how many events has 0 showers (i.e. three tracks)
  int threetrkcntr = 0; //Number of events with three tracks    
  int vectorsize3 = 0; //Number of events with 3 tracks whose start is < 5cm from the reco vertex
  int secondtrkgood = 0; //Number of events where the second shortest/longest track is contained
  int shortesttrkgood=0; //Number of events where the shortest track is contained
  int events_remaining = 0; //sanity check for number of events remaining
  int pid_cut0 = 0; //sanity pid cut
  int pid_cut1 = 0; //sannity pid cut
  int proton_cut = 0; //number of events with 2 protons in final state
  int muon_cut = 0; //number of events with 1 muon in final state
  int events_chi2p = 0; //sanity checks on chi2
  int events_chi2mu = 0; //sanity cheks on chi2
  int events_mc = 0; //sanity checks o chi2
  int events_2mu = 0; //events with two muons
  int events_2other = 0; //events with 2 others
  int n_mom_mu = 0;
  int n_mom_p1 = 0;
  int n_mom_p2 = 0;

  //neutrino counters
  int neutrinos_0 = 0;
  int neutrinos_1 = 0;
  int neutrinos_else = 0;

  //stupid counters
  int a = 0;
  int b = 0;

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

    //Defining the MC Weight cause it is dumb                                                                                                                                                                                                        
    /////////////////////////////                                                                                                                                                                                                                    
    if(std::isfinite(weightTune) && weightTune <= 100.) {
      mc_wgt = weightSplineTimesTune;
    } else {
      mc_wgt = 1 * weightSpline;
    }

    //Just casually checking how many neutrino slices we have
    if(nslice == 0){
      neutrinos_0++;
    }else if(nslice == 1){
      neutrinos_1++;
    }else{
      neutrinos_else++;
    }

    int mc_n_threshold_muon = 0;
    int mc_n_threshold_proton = 0;
    int mc_n_threshold_pion0 = 0;
    int mc_n_threshold_pionpm = 0;

    //ugh we have to fill the damn mc values:
    for ( size_t p = 0u; p < mc_pdg->size(); ++p ) {
      int pdg = mc_pdg->at( p );
      float energy = mc_E->at( p );
      if ( pdg == 13) {
	double mom = cuts.real_sqrt( std::pow(energy, 2) - std::pow(MASS_MUON, 2) );
	if ( mom > MUON_MOM_CUT ) {
	  mc_n_threshold_muon++;
	}
      }
      else if ( pdg == 2212 ) {
	double mom = cuts.real_sqrt( std::pow(energy, 2) - std::pow(MASS_PROTON, 2) );
	if ( mom > PROTON_MOM_CUT ) 
	  mc_n_threshold_proton++;
      }
      else if ( pdg == 111 ) {
	double mom = cuts.real_sqrt( std::pow(energy, 2) - std::pow(MASS_PION0, 2) );
	if ( mom > PION0_MOM_CUT) 
	mc_n_threshold_pion0++;
      }
      else if ( std::abs(pdg) == 211 ) {
	double mom = cuts.real_sqrt( std::pow(energy, 2) - std::pow(MASS_PIONPM, 2) );
	if ( mom > CHARGED_PI_MOM_CUT ) {
	  mc_n_threshold_pionpm++;
	}
      }
    }

    //Filling histograms before any selection is made
    ////////////////////////////////////////////////
    cuts.Overlay_In_FV(10,10,10,10,10,10,reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z);  // to fill the fv bool

    Fill_Histograms_Mine(0, pot_wgt*mc_wgt, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0,mc_n_threshold_pionpm,cuts.fv); //need to fix the weight
    Fill_Histograms_Raquel(0, pot_wgt*mc_wgt,cuts.fv);

    //Okay. The CCInclusive Selection requires the following things:
    // 1) The reconstructed neutrino vertex is inside the FV 
    // 2) The start point of every pfp is within the FV (i.e. contained)
    // 3) The topological score of every pfp is above 0.1 (filled per event): topological_score > 0.1
    // 4) The cosmic impact parameter is greater than 10 cm for each PFP (filled per event): CosmicIP > 10 
    // We are now goinnng to make plots of those cut variables

    //1) Check that the event is in the FV
    //////////////////////////////////////
    if(cuts.In_FV(10,10,10,10,10,10,reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z) == false) continue; //10 cm border except from backend of detector
    fvcntr++;

    //Fill Histograms
    Fill_Histograms_Mine(1, pot_wgt*mc_wgt, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0,mc_n_threshold_pionpm,cuts.fv); //need to fix the weight
    Fill_Histograms_Raquel(1, pot_wgt*mc_wgt,cuts.fv);

    //2) The start point of every pfp is within the FV
    ///////////////////////////////////////////////////
    /*for(int p = 0; p < n_pfps; p++){
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
    Fill_Histograms_Mine(2, pot_wgt*mc_wgt, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0,mc_n_threshold_pionpm,cuts.fv); //need to fix the weight
    Fill_Histograms_Raquel(2, pot_wgt*mc_wgt, cuts.fv);
    */
    //3) The topoloogical score of every neutrino slice is above 0.1
    ///////////////////////////////////////////////////////
    if(topological_score < cuts.TOPO_SCORE_CUT) continue;
    toposcore++;

    //Fill Histograms  
    Fill_Histograms_Mine(2, pot_wgt*mc_wgt, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0,mc_n_threshold_pionpm,cuts.fv); //need to fix the weight
    Fill_Histograms_Raquel(2, pot_wgt*mc_wgt, cuts.fv);
    
    //4) The cosmic impact parameter is greater than 10 cm for every neutrino slice. Honestly a dumb cut. Will remove later
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*    if(CosmicIP < cuts.COSMIC_IP_CUT) continue;
    cosmicip++;

    //Fill Histograms  
    Fill_Histograms_Mine(3, pot_wgt*mc_wgt, mc_n_threshold_muon, mc_n_threshold_proton, mc_n_threshold_pion0,mc_n_threshold_pionpm,cuts.fv); //need to fix the weight
    Fill_Histograms_Raquel(3, pot_wgt*mc_wgt, cuts.fv);
    */

    //Filling Some Cut Variables to be Used in Optimizing Cuts
    ///////////////////////////////////////////////////////////
    std::vector<int> testVector;
    for(int i = 0; i < mc_pdg->size(); i++){
      testVector.push_back(mc_pdg->at(i));
    }
    for(int i = 0; i < n_pfps; i++){
      if(i >= mc_pdg->size()){
	testVector.resize(i+1);
      }
      int track_pdg = testVector.at(i);
      Fill_Track_Plots(i,track_pdg,pot_wgt*mc_wgt);
    }
    testVector.clear();

    //Next: The Muon Selection
    //////////////////////////
   
    events_remaining++;

  } //end of Loop over events

  //Before we finish, we need to make the efficiency and purity plots:
  ///////////////////////////////////////////////////////////////////
  std::vector<int> cut_values = {static_cast<int>(nentries),fvcntr,toposcore};
  for(int i = 0; i < number; i++){
    double eff = double(cc2p0pi[i]) / double(cc2p0pi[0]);
    double purity = double(cc2p0pi[i]) / double(cut_values[i]);
    std::cout<<"Value of Efficinecy After Cut "<<i<<": "<<eff<<std::endl;
    std::cout<<"Value of Purity After Cut "<<i<<": "<<purity<<std::endl;
    eff_graph->SetPoint(i,i+1,eff);
    pur_graph->SetPoint(i,i+1,purity);
  }

  std::cout<<"-----MODULE SUMMARY-----"<<std::endl;
  std::cout << "[ANALYZER] Initial Number of Events: "<<nentries<<std::endl;
  std::cout << "[ANALYZER] Number of Events with Vertex in FV: "<<fvcntr<<std::endl;
  //std::cout << "[ANALYZER] How Many Events with All PFP Starts within the FV: "<<pfp_starts_contained<<std::endl;
  std::cout << "[ANALYZER] How Many Events with Topological Score above 0.3: "<<toposcore<<std::endl;
  //std::cout << "[ANALYZER] Number of Events with CosmicIP above 10cm: "<<cosmicip<<std::endl;
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
  std::cout<<"events_chi2mu++: "<<events_chi2mu<<std::endl;
  std::cout<<"events_chi2p++: "<<events_chi2p<<std::endl;
  std::cout<<"events_mc: "<<events_mc<<std::endl;
  std::cout<<"events_2mu: "<<events_2mu<<std::endl;    
  std::cout<<"events_2other: "<<events_2other<<std::endl;
  std::cout<<"1mu2p"<<res_count[0]<<std::endl;
  std::cout<<"1mu1p1pi"<<res_count[1]<<std::endl;
  std::cout<<"1muNp"<<res_count[2]<<std::endl;
  std::cout<<"else"<<res_count[number-1]<<std::endl;

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

  
} //end of progrm
