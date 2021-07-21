#include "constants.h"
using namespace Constants;
class Cuts{

 public:

  virtual bool Overlay_In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z); //is the true vertex in the FV?
  virtual bool In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z); //is the reco vertex in the FV?
  virtual bool True_CC2p(float x, float y, float z, int ccnc, int nu_pdg, int mc_n_threshold_muon, int mc_n_threshold_proton, int mc_n_threshold_pion0, int mc_n_threshold_pionpm); //is this a true 1mu2p event
  virtual void True_Breakdown(float x, float y, float z, int ccnc, int nu_pdg, int mc_n_threshold_muon, int mc_n_threshold_proton, int mc_n_threshold_pion0, int mc_n_threshold_pionpm, int interaction, int i); //overlay breakdown for plots
  virtual bool Event_Selection(int n_pfps, int tracks_w_good_score,int tracks_w_good_distance, int muons, int protons); //3 pfps, who are tracks, with 1mu and 2 protons
  virtual bool Reco_Momentum(TVector3 vMuon, TVector3 vLead, TVector3 vRec); //is the reco momentum of the particles okay?
  virtual bool Reco_Event(); //create the reco event boolean
  virtual void Fill_Counters(bool overlay,float x, float y, float z, int ccnc, int nu_pdg, int mc_n_threshold_muon, int mc_n_threshold_proton, int mc_n_threshold_pion0, int mc_n_threshold_pionpm,int interaction); //Fill all the Counters!

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

  //Counters for the MC Backgrounds
  /////////////////////////////////////
  static const int number = 10;

  //number of generated event/channel                                                                             
  int cc0p0pi[number+1] = {0}; //number-1 is after pid. number then would be particle specifics...should be same as number-1
  int cc1p0pi[number+1] = {0};
  int cc2p0pi[number+1] = {0};
  int ccNp0pi[number+1] = {0};
  int ccNp1pi[number+1] = {0};
  int ccNpNpi[number+1] = {0};
  int ccnue[number+1] = {0};
  int outfv[number+1] = {0};
  int nc[number+1] = {0};
  int other[number+1] = {0};
   
  //number of generated event/channel                                                                             
  int qel[number+1] = {0};
  int res[number+1] = {0};
  int mec[number+1] = {0};
  int coh[number+1] = {0};
  int dis[number+1] = {0};
  int ccnue_raquel[number+1] = {0};
  int outfv_raquel[number+1] = {0};
  int nc_raquel[number+1] = {0};
  int other_raquel[number+1] = {0};

}; //end of class definition

//Used in the Overlay to help with the MC breakdown definitions. Does same thing as above but returns fv variable
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Cuts::Overlay_In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z){
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
bool Cuts::In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z){
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
bool Cuts::True_CC2p(float x, float y, float z, int ccnc, int nu_pdg, int mc_n_threshold_muon, int mc_n_threshold_proton, int mc_n_threshold_pion0, int mc_n_threshold_pionpm){
  Overlay_In_FV(10,10,10,10,10,10,x,y,z); //check is the true vertex is in the FV      
  if(ccnc == 0 && nu_pdg == 14 && mc_n_threshold_proton == 2 && mc_n_threshold_muon == 1 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && true_fv == true){
    true_cc2p = true;
  }else{
    true_cc2p = false;
  }
  return true_cc2p;
} //end of true CC2p

//Determines if the MC Breakdown
////////////////////////////////////////
void Cuts::True_Breakdown(float x, float y, float z, int ccnc, int nu_pdg, int mc_n_threshold_muon, int mc_n_threshold_proton, int mc_n_threshold_pion0, int mc_n_threshold_pionpm, int interaction, int i){

  Overlay_In_FV(10,10,10,10,10,10,x,y,z); //check is the true vertex is in the FV

  //Now to fill My MC Breakdowns
  //////////////////////////////
  if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 0 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && true_fv == true){
    cc0p0pi[i]++;
    
    //cc1p0pi                                                                                                                      
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 1 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && true_fv == true){
    cc1p0pi[i]++;
  
    //cc2p0pi                                                                                                           
  } else if (ccnc == 0 && nu_pdg == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 2 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && true_fv == true){
    cc2p0pi[i]++;
    
    //ccNp0pi                                                                                                                       
  } else if (ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton > 2 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && true_fv == true){
    ccNp0pi[i]++;
    
    //ccNp1pi                                                                                                                      
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton >= 0 && (mc_n_threshold_pion0 == 1 || mc_n_threshold_pionpm == 1) && true_fv == true){
    ccNp1pi[i]++;
    
    //ccNpNpi                                                                                                                      
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton >= 0 && (mc_n_threshold_pion0 > 1 || mc_n_threshold_pionpm > 1) && true_fv == true){
    ccNpNpi[i]++;
    
    //CC NUE                                                                                                                       
  } else if(ccnc == 0 && abs(nu_pdg) == 12 && true_fv == true){
    ccnue[i]++;
  
    //OUT OF TRUE_FV                                                                                                                    
  } else if(true_fv == false){
    outfv[i]++;
           
  } else if(ccnc == 1 && true_fv == true){
    nc[i]++;

    //else                                                                                                                         
  } else{
    other[i]++;
  }

  //Now to Fill Raquel's MC Breakdown
  ///////////////////////////////////
   
  //CCQE
  if(ccnc == 0 && interaction == 0 && abs(nu_pdg) == 14 && true_fv==true){
    qel[i]++;
  
    //CCCoh                                                                                                                        
  } else if(ccnc == 0 && interaction == 3 && abs(nu_pdg) == 14 && true_fv == true){
    coh[i]++;

    //CCMEC                                                                                                                        
  } else if(ccnc == 0 && interaction == 10 && abs(nu_pdg) == 14 && true_fv==true){
    mec[i]++;

    //CCRES                                                                                                                        
  } else if(ccnc == 0 && interaction == 1 && abs(nu_pdg) == 14 && true_fv==true){
    //if(mc_n_threshold_muon == 1 && mc_n_threshold_proton == 2){
    //  res_count[0]++;
    //}else if (mc_n_threshold_muon == 1 && mc_n_threshold_proton == 1 && (mc_n_threshold_pion0 == 1 || mc_n_threshold_pionpm == 1)){
    //  res_count[1]++;
    //}else if(mc_n_threshold_muon == 1 && mc_n_threshold_proton > 2){
    //  res_count[2]++;
    //}else{
    //  res_count[3]++;
    //}
    res[i]++;
  
    //CCDIS                                                                                                                        
  } else if(ccnc == 0 && interaction == 2 && abs(nu_pdg) == 14 && true_fv==true){
    dis[i]++;
  
    //CCNue                                                                                                                        
  } else if(ccnc == 0  && abs(nu_pdg) == 12 && true_fv ==true){
    ccnue_raquel[i]++;
    
    //NC                                                                                                                           
  } else if(ccnc == 1 && true_fv == true){
    nc_raquel[i]++;
    
    //OUT OF TRUE_FV                                                                                                                   
  } else if(true_fv == false){
    outfv_raquel[i]++;
  
    //Other                                                                                                                        
  }else{
    other_raquel[i]++;
  } 

} //end of true_breakdown

bool Cuts::Event_Selection(int n_pfps, int tracks_w_good_score,int tracks_w_good_distance, int muons, int protons){
  
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

}//end of 1mu2p

bool Cuts::Reco_Momentum(TVector3 vMuon, TVector3 vLead, TVector3 vRec){

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

bool Cuts::Reco_Event(){
  if(pid == true && good_muon_mom == true && good_lead_mom == true && good_recoil_mom == true){
    reconstructed_event = true;
  } else {
    reconstructed_event = false;
  }
  return reconstructed_event;

}//end of reco event

void Cuts::Fill_Counters(bool overlay,float x, float y, float z, int ccnc, int nu_pdg, int mc_n_threshold_muon, int mc_n_threshold_proton, int mc_n_threshold_pion0, int mc_n_threshold_pionpm,int interaction){

  if(overlay == true){
    True_Breakdown(x,y,z,ccnc,nu_pdg,mc_n_threshold_muon,mc_n_threshold_proton,mc_n_threshold_pion0,mc_n_threshold_pionpm,interaction,0);
  }

  //fv
  if(reco_fv == true){
    fvcntr++;
    if(overlay == true){
      True_Breakdown(x,y,z,ccnc,nu_pdg,mc_n_threshold_muon,mc_n_threshold_proton,mc_n_threshold_pion0,mc_n_threshold_pionpm,interaction,1);
    }
  }

  //three pfps
  if(reco_fv && three_pfps == true){
    threepfps++;
    if(overlay == true){
      True_Breakdown(x,y,z,ccnc,nu_pdg,mc_n_threshold_muon,mc_n_threshold_proton,mc_n_threshold_pion0,mc_n_threshold_pionpm,interaction,2);
    }
  }

  //three tracks
  if(reco_fv==true && three_pfps == true && three_tracks == true){
    threetrkcntr++;
    if(overlay == true){
      True_Breakdown(x,y,z,ccnc,nu_pdg,mc_n_threshold_muon,mc_n_threshold_proton,mc_n_threshold_pion0,mc_n_threshold_pionpm,interaction,3);
    }
  }

  //three connected tracks
  if(reco_fv == true && three_pfps == true && three_tracks == true && three_tracks_4cm_vertex == true){
    threetrk_connected++;
    if(overlay == true){
      True_Breakdown(x,y,z,ccnc,nu_pdg,mc_n_threshold_muon,mc_n_threshold_proton,mc_n_threshold_pion0,mc_n_threshold_pionpm,interaction,4);
    }
  }
    
  //pid
  if(reco_fv == true && three_pfps == true && three_tracks == true && three_tracks_4cm_vertex == true && good_pid == true){
    if(overlay == true){
      True_Breakdown(x,y,z,ccnc,nu_pdg,mc_n_threshold_muon,mc_n_threshold_proton,mc_n_threshold_pion0,mc_n_threshold_pionpm,interaction,5);
    }
    pid++;
  }

  //muon_mom_cut
  if(good_pid == true && good_muon_mom == true){
    reco_muon_mom_cut++;
    if(overlay == true){
      True_Breakdown(x,y,z,ccnc,nu_pdg,mc_n_threshold_muon,mc_n_threshold_proton,mc_n_threshold_pion0,mc_n_threshold_pionpm,interaction,6);
    } 
  }

  //lead_mom_cut                                                                                                                                                                                                                
  if(good_pid == true && good_muon_mom == true && good_lead_mom == true){
    reco_lead_mom_cut++;
    if(overlay == true){
      True_Breakdown(x,y,z,ccnc,nu_pdg,mc_n_threshold_muon,mc_n_threshold_proton,mc_n_threshold_pion0,mc_n_threshold_pionpm,interaction,7);
    }
  }

  //recoil_mom_cut                                                                                                               
  if(good_pid == true && good_muon_mom == true && good_lead_mom == true && good_recoil_mom == true){
    reco_recoil_mom_cut++;
    if(overlay == true){
      True_Breakdown(x,y,z,ccnc,nu_pdg,mc_n_threshold_muon,mc_n_threshold_proton,mc_n_threshold_pion0,mc_n_threshold_pionpm,interaction,8);
    }
  }

}//end of fill counters
