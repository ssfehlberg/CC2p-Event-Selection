#ifndef CONSTANTS_H
#define CONSTANTS_H
namespace Constants{

  bool _debug = false; //debug statements                                                                       
  
  //useful masses
  double MASS_PROTON = 0.93827208;
  double MASS_MUON = 0.10565837;

  //Counters                                                                                               
  int fvcntr = 0; //Number of events with reconstructed vertex within the FV                                      
  int threepfps = 0; //Number of Events with Three PFPs
  int threetrkcntr = 0; //Number of events with three tracks    
  int threetrk_connected = 0; //Number of Events with three tracks attached to the vertex
  int pid = 0; //Number of events with 1 track with PID > 0.6 and 2 tracks with PID < 0.6
  int events_remaining = 0; //sanity check for number of events remaining

  //neutrino counters                                                                                     
  int neutrinos_0 = 0;
  int neutrinos_1 = 0;
  int neutrinos_else = 0;

}
#endif
