class helper_funcs{
 
 public:
  virtual void Overlay_In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z);
  virtual bool In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z);
  virtual double real_sqrt(double x);

  bool fv; //is event in FV?

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
