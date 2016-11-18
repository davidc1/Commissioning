#ifndef LARLITE_SIMPLECLUSTERER_CXX
#define LARLITE_SIMPLECLUSTERER_CXX

#include "SimpleClusterer.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"
#include "DataFormat/vertex.h"

namespace larlite {

  SimpleClusterer::SimpleClusterer(){

    _name        = "SimpleClusterer";
    _fout        = 0;
    _verbose     = false;
    _hitProducer = "gaushit";
    _vtxProducer = "";
    _out_clusterProducer = "rawclus";
    _useVtx      = false;
    _radius      = 2.0;
    _cellSize    = 2;
    _vtx_radius  = 0;
    _max_rms     = 100;
    _vtx_w_cm = {0,0,0};
    _vtx_t_cm = {0,0,0};
    _tick_min = 0;
    _tick_max = 9600;

  }

  bool SimpleClusterer::initialize() {

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();

    std::cout << "********************************" << std::endl;
    std::cout << "Wire -> cm conversion : " << _wire2cm << std::endl;
    std::cout << "Time -> cm conversion : " << _time2cm << std::endl;
    std::cout << "********************************" << std::endl;

    return true;
  }
  
  bool SimpleClusterer::analyze(storage_manager* storage) {

    auto evt_hits      = storage->get_data<event_hit>(_hitProducer);
    auto ev_clusters   = storage->get_data<event_cluster>(_out_clusterProducer);
    auto cluster_ass_v = storage->get_data<event_ass>(ev_clusters->name());
    auto evt_vtx       = storage->get_data<event_vertex>(_vtxProducer);
    
    //set event ID through storage manager
    storage->set_id(storage->run_id(),storage->subrun_id(),storage->event_id());

    if (!evt_hits){
      std::cout << "No hits!" << std::endl;
      return false;
    }

    if ( (_vtxProducer != "") and (!evt_vtx) ){
      std::cout << "No vertex even though one requested...quit" << std::endl;
      return false;
    }
    
    
    if ( (_vtxProducer != "") and (evt_vtx->size() == 1) ){
      auto const& vtx = evt_vtx->at(0);
      auto geoH = larutil::GeometryHelper::GetME();
      std::vector<double> xyz = {vtx.X(), vtx.Y(), vtx.Z()};
      for (size_t pl = 0; pl < 3; pl++){
	auto const& pt = geoH->Point_3Dto2D(xyz,pl);
	_vtx_w_cm[pl] = pt.w;
	_vtx_t_cm[pl] = pt.t + 800 * _time2cm;
      }
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

	// pair = (i,j) indices of this cell in the _hitMap
	auto const& pair = it->first;
	
	// wire-space cell index
	// prepare a hit list of all neighboring cells
	// _________
	// |__|__|__|
	// |__|__|__|
	// |__|__|__|
	std::vector<size_t> cellhits = it->second;

	std::vector<size_t> neighborhits;
	getNeighboringHits(pair,neighborhits);

	for (size_t h1=0; h1 < cellhits.size(); h1++){
	  // has this hit been added to a cluster?
	  // if so not necessary to look at
	  auto const& hit1 = cellhits[h1];
	  // keep track if the hit will ever be matched to another
	  bool matched = false;
	  // if not find hits it should be clustered with and add it to the appropriate cluster
	  for (size_t h2=0; h2 < neighborhits.size(); h2++){
	    auto const& hit2 = neighborhits[h2];
	    if (hit1 == hit2) continue;
	    // are the hits compatible?
	    bool compat = HitsCompatible(evt_hits->at(hit1), evt_hits->at(hit2));
	    // should the hits go in the same cluster?
	    if (compat){
	      matched = true;
	      // if both hits have already been assigned to a cluster then we can merge the cluster indices!
	      if ( (_clusterMap.find(hit1) != _clusterMap.end()) and
		   (_clusterMap.find(hit2) != _clusterMap.end()) ){
		// if in the same cluster -> do nothing
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
		}// if they are in different clusters
	      }
	      // if compatible and the 2nd hit has been added to a cluster
	      // add hit1 to the same cluster
	      else if ( (_clusterMap.find(hit2) != _clusterMap.end()) and
			(_clusterMap.find(hit1) == _clusterMap.end()) ){
		auto clusIdx = _clusterMap[hit2];
		_clusterMap[hit1] = clusIdx;
		_clusters[clusIdx].push_back(hit1);
	      }
	      // otherwise, add both to a new cluster
	      else if ( (_clusterMap.find(hit1) != _clusterMap.end()) and
			(_clusterMap.find(hit2) == _clusterMap.end()) ){
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

    }// loop through all planes

    // make a vector for the clusters
    std::vector<std::vector<unsigned int> > cluster_vector;
    for (auto it = _clusters.begin(); it != _clusters.end(); it++){
      auto indices = it->second;
      // if there are enough indices, make a cluster
      if (indices.size() >= 1){
	std::vector<unsigned int> clus;
	for (auto idx : indices)
	  clus.push_back(idx);
	cluster_vector.push_back(clus);
      }// if there are 2 hits in cluster
    }
    
    // vector for assocaitions
    std::vector<std::vector<unsigned int> > _cluster_hit_ass;
    // for each cluster create a larlite::cluster
    for (size_t i=0; i < cluster_vector.size(); i++){
      if (cluster_vector[i].size() > 0){
	larlite::cluster clus;
	clus.set_n_hits(cluster_vector[i].size());
	// grab wire/time bounds for this cluster
	float w_min = 9999;
	float t_min = 9999;
	float w_max = -1;
	float t_max = -1;
	float integral = 0;
	for (auto const& idx : cluster_vector[i]){
	  auto const& hit = evt_hits->at(idx);
	  float t = hit.PeakTime();
	  if (t > t_max) t_max = t;
	  if (t < t_min) t_min = t;
	  float w = hit.WireID().Wire;
	  if (w > w_max) w_max = w;
	  if (w < w_min) w_min = w;
	  integral += hit.Integral();
	}
	clus.set_start_wire(w_min,1.);
	clus.set_end_wire(w_max,1.);
	clus.set_start_tick(t_min,1.);
	clus.set_end_tick(t_max,1.);
	clus.set_integral(integral,0.,0.);
	clus.set_summedADC(integral, 0., 0.);
	clus.set_view( evt_hits->at( cluster_vector[i][0] ).View() );
	//clus.set_planeID( evt_hits->at( cluster_vector[i][0] ).WireID().Plane );
	//clus.set_planeID(2);
	// vector for associations
	ev_clusters->push_back(clus);
	_cluster_hit_ass.push_back(cluster_vector[i]);
      }
    }
    
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

    double dt = ( h1.PeakTime() - h2.PeakTime() ) * _time2cm;
    //  if the hit time-ranges overlap, this distnce should be 0
    if (TimeOverlap(h1,h2,dt) == true)
      dt = 0;
    double dw = ((double)h1.Channel()-(double)h2.Channel())*_wire2cm;
    double d = dt*dt + dw*dw;
    if (d > (_radius*_radius))
      return false;

    return true;
  }

  bool SimpleClusterer::TimeOverlap(const larlite::hit& h1,
				    const larlite::hit& h2,
				    double& dmin) const
    {
      bool overlap = false;

      auto T1 = h1.PeakTime() * _time2cm; // time of first hit
      auto T2 = h2.PeakTime() * _time2cm;
      auto W1 = h1.RMS() * _time2cm;
      auto W2 = h2.RMS() * _time2cm;

      if ( ( (T1 + W1) > (T2 - W2) ) and
	   ( (T1 - W1) < (T2 - W2) ) )
	return true;

      if ( ( (T1 + W1) > (T2 + W2) ) and
	   ( (T1 - W1) < (T2 + W2) ) )
	return true;

      if ( ( (T2 + W2) > (T1 - W1) ) and
	   ( (T2 - W2) < (T1 - W1) ) )
	return true;
      
      if ( ( (T2 + W2) > (T1 + W1) ) and
	   ( (T2 - W2) < (T1 + W1) ) )
	return true;

      double d = dmin;
      
      d = fabs( (T1+W1) - (T2-W2) );
      if (d < dmin)
	dmin = d;
      
      d = fabs( (T1+W1) - (T2+W2) );
      if (d < dmin)
	dmin = d;
      
      d = fabs( (T1-W1) - (T2+W2) );
      if (d < dmin)
	dmin = d;

      d = fabs( (T1-W1) - (T2-W2) );
      if (d < dmin)
	dmin = d;

      return false;
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

      // if RMS above threshold -> ignore
      if (hit.RMS() > _max_rms)
	continue;

      // remove hits with time-tick < _tick_min or > _tick_max
      if ( (hit.PeakTime() < _tick_min) or (hit.PeakTime() > _tick_max) )
	continue;
      
      double t = hit.PeakTime()*_time2cm;
      double w = hit.WireID().Wire*_wire2cm;

      // if hit too close to vertex -> ignore
      if ( (_vtxProducer != "") && (_useVtx == true) ){
	double d = sqrt( ( (t - _vtx_t_cm[plane]) * (t - _vtx_t_cm[plane]) ) +
			 ( (w - _vtx_w_cm[plane]) * (w - _vtx_w_cm[plane]) ) );
	if (d < _vtx_radius)
	  continue;
      }
      
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
