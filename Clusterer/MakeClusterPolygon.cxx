#ifndef MAKECLUSTERPOLYGON_CXX
#define MAKECLUSTERPOLYGON_CXX

#include "MakeClusterPolygon.h"
#include "DataFormat/hit.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"

namespace larlite {

  Polygon2D MakeClusterPolygon::MakePolygon(const std::vector<larlite::hit>& hits, const double& frac){
    
    auto helper = ::larutil::GeometryHelper::GetME();
    auto geo  = ::larutil::Geometry::GetME();
    double w2cm = helper->WireToCm();
    double t2cm = helper->TimeToCm();
    
    UChar_t plane = geo->ChannelToPlane(hits[0].Channel());
    
    // make vector of PxHits to then construct the cluster
    std::vector<::larutil::PxHit> pxhits;
    
    for (auto const& h : hits){
      ::larutil::PxHit pxh;
      pxh.t = h.PeakTime() * t2cm;
      pxh.w = h.WireID().Wire * w2cm;
      pxh.plane = plane;
      pxh.charge = h.Integral();
      pxh.peak = h.PeakAmplitude();
      pxhits.push_back(pxh);
    }

    // get list of points that will make up the polygon edges
    std::vector<const ::larutil::PxHit*> polygonEdges;
    helper->SelectPolygonHitList(pxhits,polygonEdges,frac);
    //now making Polygon Object
    std::pair<float,float> tmpvertex;
    //make Polygon Object as in mac/PolyOverlap.cc
    std::vector<std::pair<float,float> > vertices;
    for (unsigned int i=0; i < polygonEdges.size(); i++){
      tmpvertex = std::make_pair( polygonEdges.at(i)->w,
				  polygonEdges.at(i)->t );
      vertices.push_back( tmpvertex );
    }
    auto poly = Polygon2D( vertices );

    return poly;
  }

}
#endif
