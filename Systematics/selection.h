#include "constants.h"
using namespace Constants;
class Selection{

 public:

  virtual bool Overlay_In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z); //is the true vertex in the FV?
  virtual bool In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z); //is the reco vertex in the FV?
  virtual bool True_CC2p(float x, float y, float z, int ccnc, int nu_pdg, int mc_n_threshold_muon, int mc_n_threshold_proton, int mc_n_threshold_pion0, int mc_n_threshold_pionpm); //is this a true 1mu2p event
  virtual bool Event_Selection(int n_pfps, int tracks_w_good_score,int tracks_w_good_distance, int muons, int protons); //3 pfps, who are tracks, with 1mu and 2 protons
  virtual bool Reco_Momentum(TVector3 vMuon, TVector3 vLead, TVector3 vRec); //is the reco momentum of the particles okay?
  virtual bool Reco_Event(); //create the reco event boolean

  bool true_fv; //is the true vertex in the FV?
  bool reco_fv; //is the reconstructed vertex in the FV?
  bool true_cc2p; //is the event a true cc2p?
  bool three_pfps; //are there 3 pfps
  bool three_tracks; //are the 3 pfps tracks?
  bool three_tracks_4cm_vertex; //are the 3 tracks within 4cm of the vertex?
  bool good_pid; //are there 1 muon and 2 protons
  bool good_muon_mom; //is the reco momentum within range for muon
  bool good_lead_mom; //for lead proton?
  bool good_recoil_mom; //for recoil proton?
  bool reconstructed_event; //is this a reconstrcuted event?

}; //end of class definition

//Used in the Overlay to help with the MC breakdown definitions. Does same thing as above but returns fv variable
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Selection::Overlay_In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z){
  float_t xmin = 0.0 + x_low_edge;
  float_t xmax = 256.35 - x_up_edge;
  float_t ymin = -116.5 + y_low_edge;
  float_t ymax = 116.5 - y_up_edge;
  float_t zmin = 0.0 + z_low_edge;
  float_t zmax = 1036.8 - z_up_edge;

  if((x <= xmin || x >= xmax) || (y <= ymin || y >= ymax) || (z <= zmin || z >= zmax)){
    true_fv = false;
  } else{
    true_fv = true;
  } 
  return true_fv;
} //end of overlay in FV

//Determines if the input x,y,z is wihtin the FV. The distance from any edge can be defined by VAR_low(high)_edge
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Selection::In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z){
  float_t xmin = 0.0 + x_low_edge;
  float_t xmax = 256.35 - x_up_edge;
  float_t ymin = -116.5 + y_low_edge;
  float_t ymax = 116.5 - y_up_edge;
  float_t zmin = 0.0 + z_low_edge;
  float_t zmax = 1036.8 - z_up_edge;

  if((x <= xmin || x >= xmax) || (y <= ymin || y >= ymax) || (z <= zmin || z >= zmax)){
    reco_fv = false;
  } else{
    reco_fv = true;
  } 
  return reco_fv;
} //end of In_FV

//Determines if the Event is a True CC2p                                                                                                                                                                                                                                   ////////////////////////////////////////                                                                                                                                                                                                                                  
bool Selection::True_CC2p(float x, float y, float z, int ccnc, int nu_pdg, int mc_n_threshold_muon, int mc_n_threshold_proton, int mc_n_threshold_pion0, int mc_n_threshold_pionpm){
  Overlay_In_FV(10,10,10,10,10,10,x,y,z); //check is the true vertex is in the FV      
  if(ccnc == 0 && nu_pdg == 14 && mc_n_threshold_proton == 2 && mc_n_threshold_muon == 1 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && true_fv == true){
    true_cc2p = true;
  }else{
    true_cc2p = false;
  }
  return true_cc2p;
} //end of true CC2p


bool Selection::Event_Selection(int n_pfps, int tracks_w_good_score,int tracks_w_good_distance, int muons, int protons){
  
  //are there 3 PFPs?
  if(n_pfps == 3){
    three_pfps = true;
  } else{
    three_pfps = false;
  }

  //Are the 3 pfps tracks?
  if(n_pfps == 3 && tracks_w_good_score == 3){
    three_tracks = true;
  } else{
    three_tracks = false;
  }

  //Are the 3 tracks within 4cm of the vertex?
  if(n_pfps == 3 && tracks_w_good_score == 3 && tracks_w_good_distance == 3){
    three_tracks_4cm_vertex = true;
  }else {
    three_tracks_4cm_vertex = false;
  }

  //Are there 1 muon and 2 protons?
  if(n_pfps == 3 && tracks_w_good_score == 3 && tracks_w_good_distance == 3 && muons == 1 && protons == 2){
    good_pid = true;
  }else{
    good_pid = false;
  }

  return three_pfps && three_tracks && three_tracks_4cm_vertex && good_pid;

}//end of event Selection

bool Selection::Reco_Momentum(TVector3 vMuon, TVector3 vLead, TVector3 vRec){

    if(vMuon.Mag() > MUON_MOM_CUT_LOW && vMuon.Mag() < MUON_MOM_CUT_HIGH){
      good_muon_mom = true;
    }else{
      good_muon_mom = false;
    }

    if(vLead.Mag() > PROTON_MOM_CUT_LOW && vLead.Mag() < PROTON_MOM_CUT_HIGH){
      good_lead_mom = true;
    }else{
      good_lead_mom = false;
    }

    if(vRec.Mag() > PROTON_MOM_CUT_LOW && vRec.Mag() < PROTON_MOM_CUT_HIGH){
      good_recoil_mom = true;
    }else{
      good_recoil_mom = false;
    }
    
    return good_muon_mom && good_lead_mom && good_recoil_mom;

} //end of Reco_Momentum

bool Selection::Reco_Event(){
  if(good_pid == true && good_muon_mom == true && good_lead_mom == true && good_recoil_mom == true){
    reconstructed_event = true;
  } else {
    reconstructed_event = false;
  }
  return reconstructed_event;

}//end of reco event

