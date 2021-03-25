#include "constants.h"
using namespace Constants;
class helper_funcs{
 
 public:
  virtual void Overlay_In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z);
  virtual bool In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z);
  virtual double real_sqrt(double x);
  virtual void Calculate_Variables(bool add_protons, TVector3 vmuon, TVector3 vlead, TVector3 vrec);

  bool fv; //is event in FV?
  double muon_mom;
  double lead_mom;
  double recoil_mom;
  double muon_KE;
  double lead_KE;
  double recoil_KE;
  double muon_phi;
  double lead_phi;
  double recoil_phi;
  double muon_theta;
  double lead_theta;
  double recoil_theta;

  //Now for the interesting stuff                                                                                                                                                                                                                                   
  //////////////////////////////                                                                                                                                                                                                                                    
  //double EMuon = vMuon[3];
  //double ELead = vLead[3];
  //double ERec = vRec[3];
  //double E_tot = (EMuon + MASS_MUON) + ELead + ERec;
  /* private:

  Constants constant; 
  */

};//end of class definition


bool helper_funcs::In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z){
  //Just defining the x,y, and z locations of the FV                                                              
  float_t xmin = 0.0 + x_low_edge;
  float_t xmax = 256.35 - x_up_edge;
  float_t ymin = -116.5 + y_low_edge;
  float_t ymax = 116.5 - y_up_edge;
  float_t zmin = 0.0 + z_low_edge;
  float_t zmax = 1036.8 - z_up_edge;

  if((x <= xmin || x >= xmax) || (y <= ymin || y >= ymax) || (z <= zmin || z >= zmax)){
    return false;
  } else{
    return true;
  } 
} //end of In_FV

void helper_funcs::Overlay_In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z){
  //Just defining the x,y, and z locations of the FV                                                              
  float_t xmin = 0.0 + x_low_edge;
  float_t xmax = 256.35 - x_up_edge;
  float_t ymin = -116.5 + y_low_edge;
  float_t ymax = 116.5 - y_up_edge;
  float_t zmin = 0.0 + z_low_edge;
  float_t zmax = 1036.8 - z_up_edge;

  if((x <= xmin || x >= xmax) || (y <= ymin || y >= ymax) || (z <= zmin || z >= zmax)){
    fv = false;
  } else{
    fv = true;
  } 
} //end of Overlay_In_FV

double helper_funcs::real_sqrt( double x ) {
  if ( x < 0. ) return 0.;
  else return std::sqrt( x );
}

 void helper_funcs::Calculate_Variables(bool add_protons, TVector3 vmuon, TVector3 vlead, TVector3 vrec){

   TLorentzVector vMuon(vmuon[0],vmuon[1],vmuon[2],std::sqrt(vmuon.Mag2() + std::pow(MASS_MUON,2)));
   TLorentzVector vLead(vlead[0],vlead[1],vlead[2],std::sqrt(vlead.Mag2() + std::pow(MASS_PROTON,2)));
   TLorentzVector vRec(vrec[0],vrec[1],vrec[2],std::sqrt(vrec.Mag2() + std::pow(MASS_PROTON,2)));

    //momentum, KE, phi, theta of particles
    muon_mom = vmuon.Mag();
    lead_mom = vlead.Mag();
    recoil_mom = vrec.Mag();

    muon_KE = std::sqrt(vmuon.Mag2() + std::pow(MASS_MUON,2))-MASS_MUON;
    lead_KE = std::sqrt(vlead.Mag2() + std::pow(MASS_PROTON,2))-MASS_PROTON;
    recoil_KE = std::sqrt(vrec.Mag2() + std::pow(MASS_PROTON,2))-MASS_PROTON;

    muon_phi = vmuon.Phi();
    lead_phi = vlead.Phi();
    recoil_phi = vrec.Phi();

    muon_theta = cos(vmuon.Theta());
    lead_theta = cos(vlead.Theta());
    recoil_theta = cos(vrec.Theta());

    /*
    //Now for the interesting stuff
    //////////////////////////////
    double EMuon = vMuon[3];
    double ELead = vLead[3];
    double ERec = vRec[3];
    double E_tot = (EMuon + MASS_MUON) + ELead + ERec;
    
    //Beam Stuff
    TVector3 PT_miss(vMuon[0]+vLead[0]+vRec[0],vMuon[1]+vRec[1]+vLead[1],0);
    double Eneutrino = (EMuon+MASS_MUON) + ELead + ERec +((PT_miss.Mag2())/(2.0*35.37)) + 0.0304;
    TVector3 vBeam(0.,0.,Eneutrino); // z-direction is defined along the neutrino direction                            
    TVector3 vq = vBeam - vmuon; // Momentum transfer                                                                  
    TVector3 vmiss = vlead - vq; // Missing momentum        
    double E_tot_minus_beam = (E_tot - Eneutrino) * 1000;

    //Struck nucleon Momentum:                                                                                                                                                                                                                                      
    TVector3 vector_sum(vMuon[0] + vLead[0] + vRec[0], vMuon[1] + vLead[1] + vRec[1], vMuon[2] + vLead[2] + vRec[2]);
    TVector3 p_struck_nuc_vector(vector_sum[0], vector_sum[1], 0);
    p_struck_nuc = p_struck_nuc_vector.Mag();
    pz_tot = vlead[2] + vrec[2];

    //Stuff for STVs
    TVector3 vProton;
    if(add_protons){
      vProton.SetXYZ(vLead[0]+vRec[0],vLead[1]+vRec[1],vLead[2]+vRec[2]);
    }else{
      vProton.SetXYZ(vLead[0],vLead[1],vLead[2]);
    }
    double delta_pT = (vmuon + vProton).Perp(); //perp takes the magnitude as well. stv delta pt                                                                                                                                                                    
    double delta_phiT = std::acos( (-vmuon.X()*vProton.X() - vmuon.Y()*vProton.Y()) / (vmuon.XYvector().Mod() * vProton.XYvector().Mod())); //stv delta phi t                                                                                                       
    TVector2 delta_pT_vec = (vmuon + vProton).XYvector();
    double  delta_alphaT = std::acos( (-vmuon.X()*delta_pT_vec.X()- vmuon.Y()*delta_pT_vec.Y()) / (vmuon.XYvector().Mod() * delta_pT_vec.Mod()) ); //stv delta alpha T   
    
    //opening angles
    double open_angle_protons_lab = ((vLead[0]*vRec[0])+(vLead[1]*vRec[1])+(vLead[2]*vRec[2]))/(vlead.Mag()*vrec.Mag()); //vLead.Angle(vRec);  //note this is the cos(opening angle) opening angle between the protons                             
    double open_angle_mu_leading = ((vLead[0]*vMuon[0])+(vLead[1]*vMuon[1])+(vLead[2]*vMuon[2]))/(vlead.Mag()*vmuon.Mag()); //note this is the cos(opening angle) of the angle between the leading proton and the muon   
    double open_angle_mu_both_protona = ((vProton[0]*vMuon[0])+(vProton[1]*vMuon[1])+(vProton[2]*vMuon[2]))/(vProton.Mag()*vmuon.Mag()); //cos(opening angle) between the total proton momentum vector and the muon
    double En = std::sqrt(std::pow(NEUTRON_MASS,2) + vmiss.Mag2()); //energy of struck nucleon   
    
    //stuff for the boost to the COM
    //TLorentzVector betacm(vmiss[0] + vRec[0] + vBeam[0],vmiss[1] + vRec[1] + vBeam[1], vmiss[2] + vRec[2]+ vBeam[2], En + ERec + Eneutrino); //beta for CM              
    TLorentz muon(vMuon[0],vMuon[1],vMuon[2],vMuon[3]);
    TLorentz lead(vLead[0],vLead[1],vLead[2],vLead[3]);
    TLorentz rec(vRec[0],vRec[1],vRec[2],vRec[3]);
    TLorentzVector betacm(vRec[0]+vLead[0]+vMuon[0],vRec[1]+vLead[1]+vMuon[1],vRec[2]+vLead[2]+vMuon[2],ERec+ELead+EMuon); 
    TVector3 boost = betacm.BoostVector(); //the boost vector                                                          
    lead.Boost(-boost); //boost leading proton                                                                         
    rec.Boost(-boost); //boost recoil proton                                                                           
    muon.Boost(-boost);   
    double open_angle_protons_com = cos(lead.Angle(rec.Vect()));

    //Sanity Check
    //TLorentzVector betacm(vRec[0]+vLead[0],vRec[1]+vLead[1],vRec[2]+vLead[2],ERec+ELead);
    //TVector3 boost = betacm.BoostVector(); //the boost vector                                                                                                                                                                                                          //lead.Boost(-boost); //boost leading proton 
    //rec.Boost(-boost); //boost recoil proton                                                                                                                                                                                                                       
    //std::cout<<"Value of the added vectors x : "<<lead[0]+rec[0]<<std::endl;
    //std::cout<<"Value of the added vectors y : "<<lead[1]+rec[1]<<std::endl;
    //std::cout<<"Value of the added vectors z : "<<lead[2]+rec[2]<<std::endl;
  
    (if _debug) std::cout<<"[CALCULATE_VARIABLES] Value of PT_miss magnitude: "<<PT_miss.Mag()<<std::endl;
    (if _debug) std::cout<<"[CALCULATE_VARIABLES] Value of PT_miss magnitude2: "<<PT_miss.Mag2()<<std::endl;
    (if _debug) std::cout<<"[CALCULATE_VARIABLES] Value of PT_miss magnitude2 divided by 2*35.37: "<<(PT_miss.Mag2())/(2.0*35.37)<<std::endl;
    (if _debug) std::cout<<"[CALCULATE_VARIABLES] Value of Eneutrino: "<<Eneutrino<<std::endl;
    (if _debug) std::cout<<"[CALCULATE_VARIABLES] Value of E_tot_minus_beam: "<<E_tot_minus_beam<<std::endl;
    (if _debug) std::cout<<"[CALCULATE_VARIABLES] Value of the added vectors x : "<<lead[0]+rec[0]+muon[0]<<std::endl;
    (if _debug) std::cout<<"[CALCULATE_VARIABLES] Value of the added vectors y : "<<lead[1]+rec[1]+muon[1]<<std::endl;
    (if _debug) std::cout<<"[CALCULATE_VARIABLES] Value of the added vectors z : "<<lead[2]+rec[2]+muon[2]<<std::endl;
    */
}
