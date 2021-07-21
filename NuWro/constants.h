#ifndef CONSTANTS_H
#define CONSTANTS_H
namespace Constants{

  //debug statements
  bool _debug = false; 

  //Are we adding protons together in the STVs?
  bool add_protons = true;

  //Are we using the xsecc binning
  bool use_xsec_binning = false;

  //Cut Values
  double TRACK_SCORE_CUT = 0.8;
  double TRACK_DIST_CUT = 4;
  double PID_CUT = 0.2;

  //Momentum Cut Values (only used by the overlay)
  double MUON_MOM_CUT_LOW = 0.1;
  double MUON_MOM_CUT_HIGH = 2.5;
  double PROTON_MOM_CUT_LOW = 0.25;
  double PROTON_MOM_CUT_HIGH = 1.2;
  double CHARGED_PI_MOM_CUT = 0.065;
  double PION0_MOM_CUT = 0.065;

  //Argon 40 properties
  double MASS_TARGET = 37.215526; //GeV
  // This binding energy value is used in GENIE v3.0.6
  // double BINDING_ENERGY = 0.0295; // 40Ar, GeV
  // This value is the shell-occupancy-weighted mean of the $E_{\alpha}$ values
  // listed for 40Ar in Table II of arXiv:1609.03530. MINERvA uses an identical
  // procedure for 12C to obtain the binding energy value of 27.13 MeV, which is
  // adopted in their STV analysis described in arXiv:1910.08658.
  double BINDING_ENERGY = 0.02478; // 40Ar, GeV

  //useful masses
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
  int reco_muon_mom_cut = 0; //Number of events where reco muon momentum is < 0.1 and > 2.5
  int reco_lead_mom_cut = 0; //Number of events where reco lead momentum is < 0.25 and > 1.2  
  int reco_recoil_mom_cut = 0; //Number of events where reco recoil momentum is < 0.25 and > 1.2  
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
