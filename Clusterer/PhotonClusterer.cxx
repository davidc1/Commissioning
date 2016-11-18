#ifndef LARLITE_PHOTONCLUSTERER_CXX
#define LARLITE_PHOTONCLUSTERER_CXX

#include "PhotonClusterer.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"
#include "Linearity.h"

namespace larlite {

  PhotonClusterer::PhotonClusterer()
  {

    _name = "PhotonClusterer";
    _fout = 0;
    _clusterProducer = "";
    _maxHits = 0;
    _minHits = 0;
    _minQ    = 0;
    _max_lin_v = {0.0};
    _min_n_hits_v = {0};
    _verbose = false;
    
  }

  bool PhotonClusterer::initialize() {

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();

    std::cout << "********************************" << std::endl;
    std::cout << "Wire -> cm conversion : " << _wire2cm << std::endl;
    std::cout << "Time -> cm conversion : " << _time2cm << std::endl;
    std::cout << "********************************" << std::endl;
    
    // make sure _max_lin_v and _min_n_his_v have the same size
    if (_max_lin_v.size() != _min_n_hits_v.size()){
      std::cout << "_max_lin_v and _min_n_hits_v do not have the same size! exit..." << std::endl;
      return false;
    }

    // make sure values are increasing for nhits requirement and decreasing for linearity requirement
    for (size_t i=0; i < _max_lin_v.size() - 1; i++){
      if (_max_lin_v[i+1] < _max_lin_v[i]){
	std::cout << "_max_lin_v values decreasing! quit..." << std::endl;
	return false;
      }
    }
    
    for (size_t i=0; i < _min_n_hits_v.size() - 1; i++){
      if (_min_n_hits_v[i+1] < _min_n_hits_v[i]){
	std::cout << "_min_n_hits_v values decreasing! quit..." << std::endl;
	return false;
      }
    }

    return true;
  }
  
  bool PhotonClusterer::analyze(storage_manager* storage) {

    // load clusters
    auto ev_clus = storage->get_data<event_cluster>(_clusterProducer);

    // load associated hits
    larlite::event_hit* ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass( ev_clus->id(), ev_hit, ev_clus->name() );

    // output clusters
    auto out_clusters   = storage->get_data<event_cluster>("photon");
    auto cluster_ass_v  = storage->get_data<event_ass>(out_clusters->name());
    std::vector<std::vector<unsigned int> > cluster_hit_ass;

    //set event ID through storage manager
    storage->set_id(storage->run_id(),storage->subrun_id(),storage->event_id());

    // loop through all clusters
    for (size_t i=0; i < ev_clus->size(); i++){

      auto const& clus = ev_clus->at(i);

      auto const& ass_hits = ass_cluster_hit_v[i];

      // too many or too few hits?
      if ( (ass_hits.size() < _minHits) or (ass_hits.size() > _maxHits) )
	continue;

      // determine the linearity threshold for this cluster
      double max_lin = _max_lin_v[0];
      for (size_t n=0; n < _min_n_hits_v.size(); n++){
	auto const& min_n_hits = _min_n_hits_v[n];
	if ( ass_hits.size() > min_n_hits )
	  max_lin = _max_lin_v[n];
      }

      // grab hit list
      // get coordinates of hits to calculate linearity
      std::vector<double> hit_w_v;
      std::vector<double> hit_t_v;

      float Qtot = 0;

      bool widehits = false;
      
      for (auto const& hit_idx : ass_hits){

	auto const& hit = ev_hit->at(hit_idx);

	// ignore hits that are too wide
	if (hit.RMS() > 15) { widehits = true; break; }
	
	double w = hit.WireID().Wire  * _wire2cm;
	double t = hit.PeakTime() * _time2cm;
	double q = hit.Integral();
	Qtot += q;
	hit_w_v.push_back( w );
	hit_t_v.push_back( t );

      }// for all hits in cluster

      if (widehits) continue;

      if (Qtot < _minQ) continue;

      ::Linearity lin(hit_w_v,hit_t_v);

      // if cluster is too linear do not save
      if (lin._local_lin_truncated_avg < max_lin) continue;
      
      // made it this far. save this cluster as a photon cluster
      out_clusters->emplace_back(clus);
      cluster_hit_ass.push_back( ass_hits );
      
    }// for all input clusters
    
    cluster_ass_v->set_association(out_clusters->id(),product_id(data::kHit,ev_hit->name()), cluster_hit_ass);    
    
    return true;
  }

  bool PhotonClusterer::finalize() {

    return true;
  }

}
#endif
