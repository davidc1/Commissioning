#ifndef LARLITE_LASERFINDER_CXX
#define LARLITE_LASERFINDER_CXX

#include "LaserFinder.h"
#include "DataFormat/hit.h"
#include "DataFormat/cluster.h"
#include "Base/GeoConstants.h"

namespace larlite {

  bool LaserFinder::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //

    return true;
  }
  
  bool LaserFinder::analyze(storage_manager* storage) {
  
    //
    // Do your event-by-event analysis here. This function is called for 
    // each event in the loop. You have "storage" pointer which contains 
    // event-wise data. To see what is available, check the "Manual.pdf":
    //
    // http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
    // 
    // Or you can refer to Base/DataFormatConstants.hh for available data type
    // enum values. Here is one example of getting PMT waveform collection.
    //
    // event_fifo* my_pmtfifo_v = (event_fifo*)(storage->get_data(DATA::PMFIFO));
    //
    // if( event_fifo )
    //
    //   std::cout << "Event ID: " << my_pmtfifo_v->event_id() << std::endl;
    //
  
    // Get the hits from David C.'s hit finding:
    auto hits = storage->get_data<event_hit>("rawhit");
    auto out_cluster_v = storage->get_data<event_cluster>("laserCluster");

    // Make exactly one cluster:
    out_cluster_v -> push_back(larlite::cluster());
    out_cluster_v -> back().set_integral(hits->front().Integral(),0,0);
    out_cluster_v -> back().set_id(0);
    out_cluster_v -> back().set_view(hits->front().View());


    float timeToCm = 0.06;
    float wireToCm = 0.3;

    float laser_start_w = 3448.0;
    float laser_start_t = 5775.0;

    float laser_end_w = 2546.0;
    float laser_end_t = 8501.0;

    TVector3 pointOnLine(laser_start_w*wireToCm, laser_start_t*timeToCm,0.0);
    TVector3 directionOfLine( (laser_start_w - laser_end_w)*wireToCm, (laser_start_t - laser_end_t)*timeToCm, 0.0);

    std::vector<unsigned int > laserHits;

    // Now look at all the hits, and keep only ones that are on the line of the laser track:
    unsigned int hit_index = 0;
    for (auto & hit : * hits){
      if (hit.View() != geo::kZ) {
        hit_index ++;
        continue;
      }
      TVector3 targetPoint(hit.WireID().Wire*wireToCm, hit.PeakTime()*timeToCm,0.0);

      // Next, determine the distance of this hit from laser track:
      float dist = DistanceToLine3D(pointOnLine,directionOfLine,targetPoint);
      if (dist < 3.5){
        laserHits.push_back(hit_index);
      }
      hit_index ++;


    }
    // Make new associations
    AssSet_t hit_ass;
    // AssUnit_t new_association;  
    // new_association.push_back(laserHits);
    hit_ass.push_back(laserHits);

    auto out_ass = storage->get_data<event_ass>(out_cluster_v->name());
    out_ass->set_association(out_cluster_v->id(), hits->id(), hit_ass);


    storage -> set_id(1,0,hits->event_id());

    return true;
  }

  bool LaserFinder::finalize() {

    // This function is called at the end of event loop.
    // Do all variable finalization you wish to do here.
    // If you need, you can store your ROOT class instance in the output
    // file. You have an access to the output file through "_fout" pointer.
    //
    // Say you made a histogram pointer h1 to store. You can do this:
    //
    // if(_fout) { _fout->cd(); h1->Write(); }
    //
    // else 
    //   print(MSG::ERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");
    //
  
    return true;
  }


  float LaserFinder::DistanceToLine3D( const TVector3 & pointOnLine, const TVector3 & directionOfLine, 
                                          const TVector3 & targetPoint) const
  {
    // This algorithm finds the distance between a point and a line by finding the closest point on the line
    // Using minimization techniques from calculus.

    // Line is defined by the vectors pointOnLine and directionOfLine.
    // So, any point on the line can be parametrized as (pointOnLine + t * directionOfLine)
    // This distance between any point on the line and the target point is thus:
    // L = |(pointOnLine + t*directionOfLine) . targetPoint |
    // 
    // Using this, minimize the distance with respect to t (actually, minimize the distance squared since it's easier):
    // d(L^2)/dt = 2 * ( (pointOnLine + t*directionOfLine) . targetPoint ) * directionOfLine
    // 
    // Set equal to 0 and solve for t:
    // pointOnLine . directionOfLine + t * directionOfLine^2 - targetPoint . directionOfLine = 0;
    // Therefore:
    float t = ( targetPoint.Dot(directionOfLine) - pointOnLine.Dot(directionOfLine) ) / (directionOfLine.Dot(directionOfLine));

    // Now, construct the closest approach point, subtract the target point, and return the mag
    TVector3 closestApproach = pointOnLine + t * directionOfLine;

    closestApproach -= targetPoint;
    return closestApproach.Mag();

  }


}
#endif
