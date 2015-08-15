#ifndef LARLITE_SIMPLECLUSTERER_CXX
#define LARLITE_SIMPLECLUSTERER_CXX

#include "SimpleClusterer.h"
#include "LArUtil/GeometryUtilities.h"
#include "LArUtil/Geometry.h"

namespace larlite {

  SimpleClusterer::SimpleClusterer(){

    _name        = "SimpleClusterer";
    _fout        = 0;
    _verbose     = false;
    _hitProducer = "gaushit";
    _radius      = 0.5;

  }

  bool SimpleClusterer::initialize() {

    _wire2cm  = larutil::GeometryUtilities::GetME()->WireToCm();
    _time2cm  = 0.06;//larutil::GeometryUtilities::GetME()->TimeToCm();

    return true;
  }
  
  bool SimpleClusterer::analyze(storage_manager* storage) {

    auto evt_hits = storage->get_data<event_hit>(_hitProducer);

    if (!evt_hits){
      std::cout << "No hits!" << std::endl;
      return false;
    }

    MakeHitMap(evt_hits,_plane);

    // a map to connect hit index wih a cluster index
    // each hit gets a cluster index
    std::map<size_t, size_t> _clusterMap;

    // keep track of largest cluster ID created
    size_t maxClusterID = 0;

    // iterator for hit cell map
    std::map<std::pair<int,int>, std::vector<size_t> >::iterator it;

    // loop through hits in each cell to find matches
    for (it = _hitMap.begin(); it != _hitMap.end(); it++){

      auto const& pair    = it->first;
      // wire-space cell index


      // prepare a hit list of all neighboring cells
      
      // _________
      // |__|__|__|
      // |__|__|__|
      // |__|__|__|

      auto const neighborhits = getNeighboringHits(pair);

      // h is the hit index
      for (auto const & h : neighborhits){

	// get this hit's cluster ID?
	// if it is not in any cluster
	// create a new one
	if ( _clusterMap.find(h) == _clusterMap.end() ){
	  _clusterMap[h] = maxClusterID+1;
	  maxClusterID += 1;
	}

	


      

    return true;
  }

  bool SimpleClusterer::finalize() {

    return true;
  }

  // get all hits from neighboring cells
  std::vector<size_t> SimpleClusterer::getNeighboringHits(const std::pair<int,int> pair){
   
    std::vector<size_t> hitlist;

    auto const& i       = pair->first;
    // time-space cell index
    auto const& j       = pair->second;
    // list of hits in this cell
    cellhits = _hitMap[pair];    

    for (auto &h : hitslist)
      hitlist.push_back(h);

    std::pair<int,int> pp;

    // now look at neighboring cells, if they exist
    if (i > 0){
      pp =  std::make_pair(i-1,j);
      
    }
    
  }

  // if two hits are further apart then the set distance -> not compatible
  bool SimpleClusterer::HitsCompatible(const hit h1, const hit h2){

    double dt = (h1.PeakTime()-h2.PeakTime())*_time2cm;
    double dw = (h1.WireID().Wire-h2.WireID().Wire)*_wire2cm;

    double d = dt*dt + dw*d2;

    if (d > (_radius*_radius))
      return false;

    return true;
  }

  void SimpleClusterer::MakeHitMap(const event_hit* hitlist, int plane){

    _hitMap.clear();
    // temporary pair
    std::pair<int,int> tmpPair;
    
    for (size_t h=0; h < hitlist->size(); h++){
      auto const& hit = hitlist->at(h);
      // skip if not of plane we want
      if (hit.View() != plane)
	continue;
      double t = hit.PeakTime()*_time2cm;
      double w = hit.WireID().Wire*_wire2cm;
      // map is (i,j) -> hit list
      // i : ith bin in wire of some width
      // j : jth bin in time of some width
      int i = int(w/_cellSize);
      int j = int(t/_cellSize);
      tmpPair = std::make_pair(i,j);
      // does this entry exist in the map?
      // if yes -> append to vector
      // if no create new vector and add to map
      if (_hitMap.find(tmpPair) == _hitMap.end()){
	std::vector<size_t> aaa = {h};
	_hitMap[tmpPair] = aaa;
      }
      else
	_hitMap[tmpPair].push_back(h);
    }// for all hits

    return;
  }

}
#endif
