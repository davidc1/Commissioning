//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class larlite::SimpleClusterer+;
#pragma link C++ class larlite::LinearHitRemoval+;
#pragma link C++ class larlite::LinearClusterRemoval+;
#pragma link C++ class larlite::LinearClusterSubsetRemoval+;
#pragma link C++ class larlite::LinearClusterLocalRemoval+;
#pragma link C++ class larlite::MakeClusterPolygon+;
#pragma link C++ class larlite::MakeHits+;
#pragma link C++ class Stats+;
#pragma link C++ class larlite::ClusterFilter+;
#pragma link C++ class Linearity+;
#pragma link C++ class larlite::PhotonClusterer+;
//ADD_NEW_CLASS ... do not change this line
#endif







