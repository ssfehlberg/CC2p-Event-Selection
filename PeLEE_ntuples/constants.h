#ifndef CONSTANTS_H
#define CONSTANTS_H
namespace Constants{

  //debug statements
  bool _debug = false; 

  //Cut Values
  double TRACK_SCORE_CUT = 0.8;
  double TRACK_DIST_CUT = 4;
  double PID_CUT = 0.2;

  //Momentum Cut Values (only used by the overlay)
  double MUON_MOM_CUT = 0.1;
  double PROTON_MOM_CUT_LOW = 0.25;
  double PROTON_MOM_CUT_HIGH = 1.2;
  double CHARGED_PI_MOM_CUT = 0.065;
  double PION0_MOM_CUT = 0.065;

  //Useful masses
  double MASS_PROTON = 0.93827208; //GeV
  double MASS_MUON = 0.10565837; //GeV
  double MASS_PION0 = 0.13497666; //GeV
  double MASS_PIONPM =0.13957000; //GeV
  double MASS_NEUTRON = 0.93956541; // GeV

  //Counters                                                                                               
  ////////////////////////////////////////
  int fvcntr = 0; //Number of events with reconstructed vertex within the FV                                      
  int threepfps = 0; //Number of Events with Three PFPs
  int threetrkcntr = 0; //Number of events with three tracks    
  int threetrk_connected = 0; //Number of Events with three tracks attached to the vertex
  int pid = 0; //Number of events with 1 track with PID > 0.6 and 2 tracks with PID < 0.6
  int events_remaining = 0; //sanity check for number of events remaining

  //Neutrino counters                                                                                     
  int neutrinos_0 = 0; //# of 0 neutrino slice events
  int neutrinos_1 = 0; //# of 1 neutrino slice events
  int neutrinos_else = 0; ///# of other neutrino slice events

  //stupid counters used only by overlay
  int nue = 0; //checking the number of nue's
  int uhoh = 0; //helps to diagnose the proton id
  int contained = 0;
  int uncontained = 0;
  int denom_contained =0; //checking number of events that are uncontained and contained
  int denom_uncontained = 0;
  int num_contained = 0;
  int num_uncontained = 0;

}
#endif
