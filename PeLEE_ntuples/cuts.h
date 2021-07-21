#include "constants.h"
using namespace Constants;
class cuts{

 public:
  virtual bool Overlay_In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z);
  virtual bool In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z);
  virtual bool Three_PFPs(int n_pfps);
  virtual bool 1mu2p(int n_pfps, int tracks_w_good_score,int tracks_w_good_distance, int muons, int protons);
  virtual bool Reco_Momentum(TVector3 vMuon, TVector3 vLead, TVector3 vRec);

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



}

//Used in the Overlay to help with the MC breakdown definitions. Does same thing as above but returns fv variable
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void cuts::Overlay_In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z){
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
}

//Determines if the input x,y,z is wihtin the FV. The distance from any edge can be defined by VAR_low(high)_edge
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool cuts::In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z){
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
}

//Determines if the Event is a True CC2p
////////////////////////////////////////
bool cuts::True_Breakdown(float x, float y, float z, int ccnc, int nu_pdg, int num_muon, int num_proton, int num_pion0, int num_pionpm, int i){

  Overlay_In_FV(10,10,10,10,10,10,x,y,z); //check is the true vertex is in the FV

  //Make sure to fill your true_cc2p bool
  if(ccnc == 0 && nu_pdg == 14 && mc_n_threshold_proton == 2 && mc_n_threshold_muon == 1 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && true_fv == true){
    true_cc2p = true;
  }else{
    true_cc2p = false;
  }

  //Now to fill My MC Breakdowns
  //////////////////////////////
  if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 0 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    Fill_Mine(i,1,wgt);
    cc0p0pi[i]++;
    
    //cc1p0pi                                                                                                                      
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 1 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    cc1p0pi[i]++;
  
    //cc2p0pi                                                                                                           
  } else if (ccnc == 0 && nu_pdg == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton == 2 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    cc2p0pi[i]++;
    
    //ccNp0pi                                                                                                                       
  } else if (ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton > 2 && mc_n_threshold_pion0 == 0 && mc_n_threshold_pionpm == 0 && fv == true){
    ccNp0pi[i]++;
    
    //ccNp1pi                                                                                                                      
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton >= 0 && (mc_n_threshold_pion0 == 1 || mc_n_threshold_pionpm == 1) && fv == true){
    ccNp1pi[i]++;
    
    //ccNpNpi                                                                                                                      
  } else if(ccnc == 0 && abs(nu_pdg) == 14 && mc_n_threshold_muon == 1 && mc_n_threshold_proton >= 0 && (mc_n_threshold_pion0 > 1 || mc_n_threshold_pionpm > 1) && fv == true){
    ccNpNpi[i]++;
    
    //CC NUE                                                                                                                       
  } else if(ccnc == 0 && abs(nu_pdg) == 12 && fv == true){
    ccnue[i]++;
  
    //OUT OF FV                                                                                                                    
  } else if(fv == false){                                                                                                               outfv[i]++;
           
  } else if(ccnc == 1 && fv == true){
    nc[i]++;

    //else                                                                                                                         
  } else{
    other[i]++;
  
  }

  //Now to Fill Raquel's MC Breakdown
  ///////////////////////////////////
   
  //CCQE
  if(ccnc == 0 && interaction == 0 && abs(nu_pdg) == 14 && fv==true){
    qel[i]++;
  
    //CCCoh                                                                                                                        
  } else if(ccnc == 0 && interaction == 3 && abs(nu_pdg) == 14 && fv == true){
    coh[i]++;

    //CCMEC                                                                                                                        
  } else if(ccnc == 0 && interaction == 10 && abs(nu_pdg) == 14 && fv==true){
    mec[i]++;

    //CCRES                                                                                                                        
  } else if(ccnc == 0 && interaction == 1 && abs(nu_pdg) == 14 && fv==true){
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
  } else if(ccnc == 0 && interaction == 2 && abs(nu_pdg) == 14 && fv==true){
    dis[i]++;
  
    //CCNue                                                                                                                        
  } else if(ccnc == 0  && abs(nu_pdg) == 12 && fv ==true){
    ccnue_raquel[i]++;
    
    //NC                                                                                                                           
  } else if(ccnc == 1 && fv == true){
    nc_raquel[i]++;
    
    //OUT OF FV                                                                                                                   
  } else if(fv == false){
    outfv_raquel[i]++;
  
    //Other                                                                                                                        
  }else{
    other_raquel[i]++;
  } 

} //end of true_breakdown

bool cuts::1mu2p(int n_pfps, int tracks_w_good_score,int tracks_w_good_distance, int muons, int protons){
  
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
    good_pid == true;
  }else{
    good_pid == false;
  }

}//end of 1mu2p

bool cuts::Reco_Momentum(TVector3 vMuon, TVector3 vLead, TVector3 vRec){

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

} //end of Reco_Momentum

bool cuts::Fill_Counters(){

    //fv
  if(reco_in_FV == true){
    fvcntr++;
    if(true_cc2p == true){
      cc2p0pi[1]++;
    }
  }

  //three pfps
  if(reco_in_FV && three_pfps == true){
    threepfps++;
    if(true_cc2p == true){
      cc2p0pi[2]++;
    }
  }

  //three tracks
  if(reco_in_FV==true && three_pfps == true && three_tracks == true){
    threetrkcntr++;
    if(true_cc2p == true){
      cc2p0pi[3]++;
    }
  }

  //three connected tracks
  if(reco_in_FV == true && three_pfps == true && three_tracks == true && three_tracks_4cm_vertex == true){
    threetrk_connected++;
    if(true_cc2p == true){
      cc2p0pi[4]++;
    }
  }
    
  //pid
  if(reco_in_FV == true && three_pfps == true && three_tracks == true && three_tracks_4cm_vertex == true && reco_cc2p_event == true){
    if(true_cc2p == true){
      cc2p0pi[5]++;
    }
    pid++;
  }

  //muon_mom_cut
  if(reconstructed_event == true && good_muon_mom == true){
    reco_muon_mom_cut++;
    if(true_cc2p == true){
      cc2p0pi[6]++;
    } 
  }

  //lead_mom_cut                                                                                                                                                                                                                
  if(reconstructed_event == true && good_muon_mom == true && good_lead_mom == true){
    reco_lead_mom_cut++;
    if(true_cc2p == true){
      cc2p0pi[7]++;
    }
  }

  //recoil_mom_cut                                                                                                               
  if(reconstructed_event == true && good_muon_mom == true && good_lead_mom == true && good_recoil_mom == true){
    reco_recoil_mom_cut++;
    if(true_cc2p == true){
      cc2p0pi[8]++;
    }
  }


}//end of fill counters
