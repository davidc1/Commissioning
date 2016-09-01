#ifndef LARLITE_LINEARCLUSTERLOCALREMOVAL_CXX
#define LARLITE_LINEARCLUSTERLOCALREMOVAL_CXX

#include "LinearClusterLocalRemoval.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"

namespace larlite {

  LinearClusterLocalRemoval::LinearClusterLocalRemoval()
    : _tree(nullptr)
  {

    _name        = "LinearClusterLocalRemoval";
    _fout        = 0;
    _verbose     = false;
    _radius      = 2.0;
    _clusterProducer = "rawcluster";
    _max_lin     = 1.0;
    _min_n_hits  = 100;
  }

  bool LinearClusterLocalRemoval::initialize() {

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
  
  bool LinearClusterLocalRemoval::analyze(storage_manager* storage) {

    auto ev_clus = storage->get_data<event_cluster>(_clusterProducer);

    larlite::event_hit* ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());
    
    auto out_hit = storage->get_data<event_hit>("shrhits2");
    
    //set event ID through storage manager
    storage->set_id(storage->run_id(),storage->subrun_id(),storage->event_id());

    if (!ev_hit){
      std::cout << "No hits!" << std::endl;
      return false;
    }

    // loop trhough each cluster and calculate linaerity
    // if above some thresdhold, remove cluster

    for (size_t i=0; i < ass_cluster_hit_v.size(); i++){

      auto hit_idx_v = ass_cluster_hit_v[i];

      _n_hits = hit_idx_v.size();

      if (hit_idx_v.size() < _min_n_hits){
	for (auto const& hit_idx : hit_idx_v)
	  out_hit->emplace_back( ev_hit->at( hit_idx ) );
	continue;
      }

      std::cout << "n hits : " << hit_idx_v.size() << std::endl;

      // for each hit in the cluster, determine the local linearity
      for (auto const& hit_idx : hit_idx_v){
	
	// first find the local neighbors
	std::vector<unsigned int> local_hit_idx_v;
	getNeighboringHits(hit_idx, hit_idx_v, local_hit_idx_v, ev_hit);
      
	// get coordinates of hits to calculate linearity
	std::vector<double> hit_w_v;
	std::vector<double> hit_t_v;

	_pl = ev_hit->at(hit_idx).WireID().Plane;

	for (auto const& hit_idx2 : local_hit_idx_v){
	  hit_w_v.push_back( ev_hit->at(hit_idx2).WireID().Wire  * _wire2cm );
	  hit_t_v.push_back( ev_hit->at(hit_idx2).PeakTime() * _time2cm );
	}

	// calculate covariance
	double L = linearity(hit_w_v,hit_t_v);

	bool remove = false;

	_l = log(L);
	_tree->Fill();

	
	if ( log(fabs(L)) < _max_lin)
	  remove = true;
	
	if (remove == false)
	  out_hit->emplace_back( ev_hit->at( hit_idx ) );
	
      }// loop through all hits in cluster

    }// loop through all clusters
    
    return true;
  }

  bool LinearClusterLocalRemoval::finalize() {

    _fout->cd();

    std::cout << "tree entries : " << _tree->GetEntries() << std::endl;
    _tree->Write();

    return true;
  }

  // covariance
  double LinearClusterLocalRemoval::cov (const std::vector<double>& data1,
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
  
  double LinearClusterLocalRemoval::stdev(const std::vector<double>& data) const
  {
    if(data.size() == 0) { std::cout << "zero-vector!" << std::endl; return 0; }

    double result = 0.0;
    auto    avg   = mean(data);
    for(const auto& d: data)
      result += (d - avg)*(d - avg);
    
    return sqrt(result/((double)data.size()));
  }
  
  double LinearClusterLocalRemoval::mean(const std::vector<double>& data) const
  {
    if(data.size() == 0) { std::cout << "zero-vector!" << std::endl; return 0; }
	
    double result = 0.0;

    for(const auto& d : data) 
      result += d;
        
    return (result / ((double)data.size()));
  }

  double LinearClusterLocalRemoval::linearity(const std::vector<double>& data1,
					      const std::vector<double>& data2) const
  {

    if (data1.size() < 2)
      return 1.;

    auto C  = cov(data1,data2);
    auto sW = cov(data1,data1);
    auto sT = cov(data2,data2);

    double r_num = C;
    double r_den = sqrt( sW * sT );
    double r = 0.;

    if (r_den == 0)
      r = 0.;
    else
      r = r_num / r_den;
    if (r > 1.) r = 1.;
    if (r < -1) r = -1;

    double n = data1.size() - 2 ;

    double lin = sqrt( (1-r*r) * sT / sW / n );

    return lin;

  }


  void LinearClusterLocalRemoval::getNeighboringHits(const unsigned int& hit_idx,
						     const std::vector<unsigned int>& hit_idx_v,
						     std::vector<unsigned int>& out_hit_v,
						     const event_hit* ev_hit) const
  {

    out_hit_v.clear();

    out_hit_v.push_back( hit_idx );

    auto const& hit = ev_hit->at( hit_idx );
    double t = hit.PeakTime() * _time2cm;
    double w = hit.WireID().Wire * _wire2cm;

    for (auto const& hit_idx2 : hit_idx_v){

      if (hit_idx2 == hit_idx)
	continue;

      auto const& hit2 = ev_hit->at( hit_idx2 );

      double t2 = hit2.PeakTime() * _time2cm;
      double w2 = hit2.WireID().Wire * _wire2cm;
      
      double d = sqrt( ( (w2-w) * (w2-w) + (t2-t)*(t2-t) ) );

      if (d < _radius)
	out_hit_v.push_back( hit_idx2 );

    }// for all hits in cluster

    return;
  }
  
}
#endif
