#ifndef LARLITE_SIMPLECLUSTERER_CXX
#define LARLITE_SIMPLECLUSTERER_CXX

#include "SimpleClusterer.h"
#include "LArUtil/GeometryUtilities.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"

namespace larlite {

  SimpleClusterer::SimpleClusterer(){

    _name        = "SimpleClusterer";
    _fout        = 0;
    _verbose     = false;
    _hitProducer = "gaushit";
    _radius      = 2.0;
    _cellSize    = 2;

  }

  bool SimpleClusterer::initialize() {

    _wire2cm  = larutil::GeometryUtilities::GetME()->WireToCm();
    _time2cm  = 0.06;//larutil::GeometryUtilities::GetME()->TimeToCm();

    return true;
  }
  
  bool SimpleClusterer::analyze(storage_manager* storage) {

    auto evt_hits = storage->get_data<event_hit>(_hitProducer);
    auto ev_clusters = storage->get_data<event_cluster>("rawcluster");
    auto cluster_ass_v = storage->get_data<event_ass>(ev_clusters->name());

    //set event ID through storage manager
    storage->set_id(storage->get_data<event_hit>(_hitProducer)->run(),
		    storage->get_data<event_hit>(_hitProducer)->subrun(),
		    storage->get_data<event_hit>(_hitProducer)->event_id()); 

    if (!evt_hits){
      std::cout << "No hits!" << std::endl;
      return false;
    }

    // a map to connect hit index wih a cluster index
    // each hit gets a cluster index
    // _clusterMap[hit_index] -> cluster_index
    std::map<size_t, size_t> _clusterMap;
    // a map to connect the cluster index with the vector of hit indices for that cluster
    // _clusters[index] -> vector of hit indices for that cluster
    std::map<size_t,std::vector<size_t> > _clusters;

    // keep track of largest cluster ID created
    size_t maxClusterID = 0;

    for (int pl=0; pl < 3; pl++){
      
      MakeHitMap(evt_hits,pl);
      

      
      // iterator for hit cell map
      std::map<std::pair<int,int>, std::vector<size_t> >::iterator it;
      
      // loop through hits in each cell to find matches
      for (it = _hitMap.begin(); it != _hitMap.end(); it++){

	// keep track of compatible hits
	int numCompat = 0;
	
	auto const& pair = it->first;
	
	// wire-space cell index
	// prepare a hit list of all neighboring cells
	// _________
	// |__|__|__|
	// |__|__|__|
	// |__|__|__|
	std::vector<size_t> neighborhits;
	getNeighboringHits(pair,neighborhits);
	for (size_t h1=0; h1 < neighborhits.size(); h1++){
	  // has this hit been added to a cluster?
	  // if so not necessary to look at
	  auto const& hit1 = neighborhits[h1];
	  // keep track if the hit will ever be matched to another
	  bool matched = false;
	  // if not find hits it should be clustered with and add it to the appropriate cluster
	  for (size_t h2=h1+1; h2 < neighborhits.size(); h2++){
	    auto const& hit2 = neighborhits[h2];
	    // are the hits compatible?
	    bool compat = HitsCompatible(evt_hits->at(h1), evt_hits->at(h2));
	    // should the hits go in the same cluster?
	    if (compat){
	      //std::cout << "hits compatible!" << std::endl;
	      matched = true;
	      numCompat += 1;
	      // if both hits have already been assigned to a cluster then we can merge the cluster indices!
	      if ( (_clusterMap.find(hit1) != _clusterMap.end()) and
		   (_clusterMap.find(hit2) != _clusterMap.end()) ){
		// if they are in different clusters:
		if (_clusterMap[hit1] != _clusterMap[hit2]){
		  auto idx1 = _clusterMap[hit1];
		  auto idx2 = _clusterMap[hit2];
		  // hit indices for 1st cluster:
		  auto hits1 = _clusters[idx1];
		  auto hits2 = _clusters[idx2];
		  // append hits2 to hits1
		  for (auto h : hits2){
		    hits1.push_back(h);
		    // also change the index that the hit goes to (idx1 instead of idx2)
		    _clusterMap[h] = idx1;
		  }
		  _clusters[idx1] = hits1;
		  // erase cluster @ index2
		  _clusters.erase(idx2);
		}
	      }
	      // if compatible and the 2nd hit has been added to a cluster
	      // add hit1 to the same cluster
	      else if (_clusterMap.find(hit2) != _clusterMap.end()){
		auto clusIdx = _clusterMap[hit2];
		_clusterMap[hit1] = clusIdx;
		_clusters[clusIdx].push_back(hit1);
	      }
	      // otherwise, add both to a new cluster
	      else if (_clusterMap.find(hit1) != _clusterMap.end()){
		auto clusIdx = _clusterMap[hit1];
		_clusterMap[hit2] = clusIdx;
		_clusters[clusIdx].push_back(hit2);
	      }
	      // if neither has a cluster yet
	      else{
		// create a new cluster for this match
		_clusterMap[hit1] = maxClusterID;
		_clusterMap[hit2] = maxClusterID;
		std::vector<size_t> cl = {hit1,hit2};
		_clusters[maxClusterID] = cl;
		maxClusterID += 1;
	      }
	      //std::cout << std::endl;
	    }// if the two hits are compatible
	  }// 2nd loop through hits in the cell
	  // has this hit been matched? if not we still need to add it as its own cluster
	  /*
	    if (matched == false){
	    _clusterMap[hit1] = maxClusterID;
	    maxClusterID += 1;
	    }
	  */
	}// 1st loop through hits in the cell
      }// loop through all cells

    }

    // make a vector for the clusters
    std::vector<std::vector<unsigned int> > cluster_vector;
    for (auto it = _clusters.begin(); it != _clusters.end(); it++){
      auto indices = it->second;
      // if there are enough indices, make a cluster
      if (indices.size() > 2){
	std::vector<unsigned int> clus;
	for (auto idx : indices)
	  clus.push_back(idx);
	cluster_vector.push_back(clus);
      }// if there are 2 hits in cluster
    }
    
    //std::cout << "number of clusters: " << _clusters.size() << std::endl;
    //std::cout << "number of clusters: " << cluster_vector.size() << std::endl;
    
    // vector for assocaitions
    std::vector<std::vector<unsigned int> > _cluster_hit_ass;
    // for each cluster create a larlite::cluster
    for (size_t i=0; i < cluster_vector.size(); i++){
      if (cluster_vector[i].size() > 0){
	larlite::cluster clus;
	// vector for associations
	ev_clusters->push_back(clus);
	_cluster_hit_ass.push_back(cluster_vector[i]);
      }
    }
    
    //std::cout << "number of larlite clusters: " << ev_clusters->size() << std::endl;
    //std::cout << "number of clusters to be saved: " << _cluster_hit_ass.size() << std::endl;
    cluster_ass_v->set_association(ev_clusters->id(),product_id(data::kHit,evt_hits->name()),_cluster_hit_ass);    

    return true;
  }

  bool SimpleClusterer::finalize() {

    return true;
  }

  // get all hits from neighboring cells
  void SimpleClusterer::getNeighboringHits(const std::pair<int,int>& pair, std::vector<size_t>& hitIndices){
   
    auto const& i       = pair.first;
    // time-space cell index
    auto const& j       = pair.second;

    // _________
    // |__|__|__|
    // |__|XX|__|
    // |__|__|__|
    if (_hitMap.find(std::make_pair(i,j)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i,j)])
	hitIndices.push_back(h);
    }

    // now look at neighboring cells, if they exist
    // _________
    // |__|__|__|
    // |XX|__|__|
    // |__|__|__|
    if (_hitMap.find(std::make_pair(i-1,j)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i-1,j)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|__|
    // |__|XX|__|
    if (_hitMap.find(std::make_pair(i,j-1)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i,j-1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|__|
    // |XX|__|__|
    if ( _hitMap.find(std::make_pair(i-1,j-1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i-1,j-1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|XX|__|
    // |__|__|__|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i,j+1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|XX|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i+1,j)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|XX|
    // |__|__|__|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i+1,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j+1)])
	hitIndices.push_back(h);
    }
    // _________
    // |XX|__|__|
    // |__|__|__|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i-1,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i-1,j+1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|__|
    // |__|__|XX|
    if ( _hitMap.find(std::make_pair(i+1,j-1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j-1)])
	hitIndices.push_back(h);
    }
  }

  // if two hits are further apart then the set distance -> not compatible
  bool SimpleClusterer::HitsCompatible(const hit& h1, const hit& h2){

    if (h1.WireID().Plane != h2.WireID().Plane)
      return false;

    double dt = (h1.PeakTime()-h2.PeakTime())*_time2cm;
    double dw = ((double)h1.Channel()-(double)h2.Channel())*_wire2cm;
    double d = dt*dt + dw*dw;
    if (d > (_radius*_radius))
      return false;

    return true;
  }

  void SimpleClusterer::MakeHitMap(const event_hit* hitlist, int plane){

    _hitMap.clear();
    // temporary pair
    std::pair<int,int> tmpPair;

    // set maximum cell size
    _maxI = 0;
    _maxJ = 0;
    
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
      if (i > _maxI) _maxI = i;
      if (j > _maxJ) _maxJ = j;
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
