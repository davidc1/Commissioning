#ifndef LARLITE_LINEARHITREMOVAL_CXX
#define LARLITE_LINEARHITREMOVAL_CXX

#include "LinearHitRemoval.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"

namespace larlite {

  LinearHitRemoval::LinearHitRemoval()
    : _tree(nullptr)
  {

    _name        = "LinearHitRemoval";
    _fout        = 0;
    _verbose     = false;
    _clusProducer = "";
    _radius      = 2.0;
    _cellSize    = 2;
    _max_lin     = 0.7;

  }

  bool LinearHitRemoval::initialize() {
    
    if (_tree) delete _tree;
    _tree = new TTree("_tree","tree");
    _tree->Branch("_l",&_l,"l/D");
    _tree->Branch("_n_hits",&_n_hits,"n_hits/I");
    _tree->Branch("_pl",&_pl,"pl/I");
    
    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();

    std::cout << "********************************" << std::endl;
    std::cout << "Wire -> cm conversion : " << _wire2cm << std::endl;
    std::cout << "Time -> cm conversion : " << _time2cm << std::endl;
    std::cout << "********************************" << std::endl;

    return true;
  }
  
  bool LinearHitRemoval::analyze(storage_manager* storage) {

    auto ev_clus = storage->get_data<event_cluster>(_clusProducer);
    larlite::event_hit* ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());

    auto out_hits = storage->get_data<event_hit>("shrhits2");
    // produced clusters
    auto ev_clusters   = storage->get_data<event_cluster>("shrcluster2");
    auto cluster_ass_v = storage->get_data<event_ass>("shrcluster2");
    std::vector<std::vector<unsigned int> > cluster_hit_ass;
    
    //set event ID through storage manager
    storage->set_id(storage->run_id(),storage->subrun_id(),storage->event_id());

    if (!ev_hit){
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

      std::vector<unsigned int> ass_hit_idx_v;
      
      MakeHitMap(ass_cluster_hit_v,ev_hit,pl);
      
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

	  // get all hits within a small radius of this one
	  // measure local linarity
	  // if above some value, remove the hit
	  std::vector<double> hit_w_v;
	  std::vector<double> hit_t_v;

	  hit_w_v.push_back( ev_hit->at(hit1).Channel()  * _wire2cm );
	  hit_t_v.push_back( ev_hit->at(hit1).PeakTime() * _time2cm );
	  
	  for (size_t h2=0; h2 < neighborhits.size(); h2++){
	    auto const& hit2 = neighborhits[h2];
	    if (hit1 == hit2) continue;
	    if (HitsCompatible( ev_hit->at(hit1), ev_hit->at(hit2) ) == true){
	      hit_w_v.push_back( ev_hit->at(hit2).Channel()  * _wire2cm );
	      hit_t_v.push_back( ev_hit->at(hit2).PeakTime() * _time2cm );
	    }
	    
	  }// for all hits in neighboring cells

	  bool shrlike = true;

	  if (hit_w_v.size() > 2){
	    
	    double L = linearity(hit_w_v,hit_t_v);

	    _l = L;
	    _tree->Fill();

	    if (L > _max_lin)
	      shrlike = false;
	    
	  }// if enough hits
	  
	  if (shrlike == true){
	    ass_hit_idx_v.push_back( hit1 );
	    out_hits->emplace_back( ev_hit->at(hit1) );
	  }
	  
	}// 1st loop through hits in the cell
      }// loop through all cells

      larlite::cluster clus;
      clus.set_n_hits( ass_hit_idx_v.size() );
      cluster_hit_ass.push_back( ass_hit_idx_v );
      ev_clusters->emplace_back( clus );
      
    }// loop through all planes

    cluster_ass_v->set_association(ev_clusters->id(), product_id(data::kHit,ev_hit->name()),cluster_hit_ass);

    return true;
  }

  bool LinearHitRemoval::finalize() {

    _fout->cd();

    _tree->Write();

    return true;
  }

  // get all hits from neighboring cells
  void LinearHitRemoval::getNeighboringHits(const std::pair<int,int>& pair, std::vector<size_t>& hitIndices){
   
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

  // covariance
    double LinearHitRemoval::cov (const std::vector<double>& data1,
				  const std::vector<double>& data2) const
  {
    if(data1.size() == 0) { std::cout << "zero-vector!" << std::endl; return 0; }
    if(data2.size() == 0) { std::cout << "zero-vector!" << std::endl; return 0; }

    double result = 0.0;
    auto   mean1  = mean(data1);
    auto   mean2  = mean(data2);
    
    for(size_t i = 0; i < data1.size(); ++i)
      result += (data1[i] - mean1)*(data2[i] - mean2);
    
    return result/((double)data1.size());
      
  }
  
  double LinearHitRemoval::stdev(const std::vector<double>& data) const
  {
    if(data.size() == 0) { std::cout << "zero-vector!" << std::endl; return 0; }

    double result = 0.0;
    auto    avg   = mean(data);
    for(const auto& d: data)
      result += (d - avg)*(d - avg);
    
    return sqrt(result/((double)data.size()));
  }
  
  double LinearHitRemoval::mean(const std::vector<double>& data) const
  {
    if(data.size() == 0) { std::cout << "zero-vector!" << std::endl; return 0; }
	
    double result = 0.0;

    for(const auto& d : data) 
      result += d;
        
    return (result / ((double)data.size()));
  }


  // if two hits are further apart then the set distance -> not compatible
  bool LinearHitRemoval::HitsCompatible(const hit& h1, const hit& h2){

    if (h1.WireID().Plane != h2.WireID().Plane)
      return false;

    double dt = (h1.PeakTime()-h2.PeakTime())*_time2cm;
    double dw = ((double)h1.Channel()-(double)h2.Channel())*_wire2cm;
    double d = dt*dt + dw*dw;
    if (d > (_radius*_radius))
      return false;

    return true;
  }

  void LinearHitRemoval::MakeHitMap(const std::vector<std::vector<unsigned int> >& ass_v,
				    const event_hit* hitlist, int plane){

    _hitMap.clear();
    // temporary pair
    std::pair<int,int> tmpPair;

    for (auto const& hit_idx_v : ass_v){
      for (auto const& h : hit_idx_v){
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
    }// for all clusters
      
    return;
  }

  double LinearHitRemoval::linearity(const std::vector<double>& data1,
				     const std::vector<double>& data2) const
  {
    
    if (data1.size() < 2)
      return 1.;

    auto C  = cov(data1,data2);
    auto sW = cov(data1,data1);
    auto sT = cov(data2,data2);

    //return C / (sW*sT);

    double r_num = C;
    double r_den = sqrt( sW * sT );
    double r = 0.;

    if (r_den == 0)
      r = 0.;
    else
      r = r_num / r_den;
    if (r > 1.) r = 1.;
    if (r < -1) r = -1;

    return fabs(r);

    //double n = data1.size() - 2 ;

    //double lin = sqrt( (1-r*r) * sT / sW / n );

    //return lin;
  }


}
#endif
