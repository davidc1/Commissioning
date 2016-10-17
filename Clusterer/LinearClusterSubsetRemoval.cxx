#ifndef LARLITE_LINEARCLUSTERSUBSETREMOVAL_CXX
#define LARLITE_LINEARCLUSTERSUBSETREMOVAL_CXX

#include "LinearClusterSubsetRemoval.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"
#include "DataFormat/vertex.h"

namespace larlite {

  LinearClusterSubsetRemoval::LinearClusterSubsetRemoval()
    : _tree(nullptr)
  {

    _name        = "LinearClusterSubsetRemoval";
    _fout        = 0;
    _verbose     = false;
    _clusterProducer = "";
    _vtxProducer     = "";
    _out_hitProducer = "";
    _max_lin_v = {0.0};
    _min_n_hits_v = {0};

    _vtx_w_cm = {0,0,0};
    _vtx_t_cm = {0,0,0};

  }

  bool LinearClusterSubsetRemoval::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("linearclusterremoval","Linear Cluster Removal TTree");
    _tree->Branch("_nhits",&_nhits,"nhits/I");
    _tree->Branch("_lin",  &_lin  ,"lin/D"  );
    _tree->Branch("_local_lin_truncated",  &_local_lin_truncated  ,"local_lin_truncated/D"  );
    _tree->Branch("_local_lin_avg",  &_local_lin_avg  ,"local_lin_avg/D"  );

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
  
  bool LinearClusterSubsetRemoval::analyze(storage_manager* storage) {

    auto ev_clus = storage->get_data<event_cluster>( _clusterProducer );
    auto evt_vtx = storage->get_data<event_vertex> ( _vtxProducer     );

    larlite::event_hit* ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());
    
    auto out_hit = storage->get_data<event_hit>( _out_hitProducer );
    
    //set event ID through storage manager
    storage->set_id(storage->run_id(),storage->subrun_id(),storage->event_id());

    if (!ev_hit){
      std::cout << "No hits!" << std::endl;
      return false;
    }

    // get vertex position on each plane
    if ( (evt_vtx->size() == 1) ){
      auto const& vtx = evt_vtx->at(0);
      auto geoH = larutil::GeometryHelper::GetME();
      std::vector<double> xyz = {vtx.X(), vtx.Y(), vtx.Z()};
      for (size_t pl = 0; pl < 3; pl++){
	auto const& pt = geoH->Point_3Dto2D(xyz,pl);
	_vtx_w_cm[pl] = pt.w;
	_vtx_t_cm[pl] = pt.t + 800 * _time2cm;
      }
    }

    // keep track of all hit indices to be removed
    std::vector<size_t> hit_idx_to_be_removed_v;

    // loop trhough each cluster and calculate linaerity
    // if above some thresdhold, remove cluster

    for (size_t i=0; i < ass_cluster_hit_v.size(); i++){

      // store output cluster hit indices
      std::vector<unsigned int> out_cluster_hit_idx_v;
      larlite::cluster out_clus;

      auto hit_idx_v = ass_cluster_hit_v[i];

      bool remove = false;

      // determine the linearity threshold for this cluster
      double max_lin = _max_lin_v[0];
      for (size_t n=0; n < _min_n_hits_v.size(); n++){
	auto const& min_n_hits = _min_n_hits_v[n];
	if ( hit_idx_v.size() > min_n_hits )
	  max_lin = _max_lin_v[n];
      }

      int pl = ev_hit->at(hit_idx_v[0]).WireID().Plane;
      
      // get coordinates of hits to calculate linearity
      std::vector<double> hit_w_v;
      std::vector<double> hit_t_v;

      // minimum distance to vtx
      double dmin = 1000.;
      
      for (auto const& hit_idx : hit_idx_v){
	double w = ev_hit->at(hit_idx).WireID().Wire  * _wire2cm;
	double t = ev_hit->at(hit_idx).PeakTime() * _time2cm;
	hit_w_v.push_back( w );
	hit_t_v.push_back( t );
	double d = ( w - _vtx_w_cm[pl] ) * ( w - _vtx_w_cm[pl] ) + ( t - _vtx_t_cm[pl] ) * ( t - _vtx_t_cm[pl] );
	if (d < dmin) dmin = d;
      }
      dmin = sqrt(dmin);
      
      ::Linearity lin(hit_w_v,hit_t_v);

      _nhits = hit_w_v.size();
      _lin                 = lin._lin;
      _local_lin_avg       = lin._local_lin_avg;
      _local_lin_truncated = lin._local_lin_truncated_avg;

      // get dot product between slope line and mean to vtx line
      double _vtx_mean_w = _vtx_w_cm[pl] - lin._meanx;
      double _vtx_mean_t = _vtx_t_cm[pl] - lin._meany;
      double mag = sqrt( _vtx_mean_w * _vtx_mean_w + _vtx_mean_t * _vtx_mean_t );
      _vtx_mean_w /= mag;
      _vtx_mean_t /= mag;
      
      double slope     = lin._slope;
      double intercept = lin._intercept;

      // impact parameter to vertex:
      double x0 = _vtx_w_cm[pl];
      double y0 = _vtx_t_cm[pl];
      double IP = fabs( - slope * x0 + y0 - intercept ) / sqrt( slope * slope + 1 );

      _tree->Fill();

      if (_local_lin_truncated < max_lin and !( (dmin > 9.0) and (IP < 3.0) ) )
	remove = true;

      if ( (IP > 20) and (_nhits > 20) and (_local_lin_truncated < 0.3) )
	remove = true;

      if (remove){
	for (auto const& hit_idx : hit_idx_v)
	  hit_idx_to_be_removed_v.push_back( hit_idx );
      }
      
      
    }// loop through all planes

    // sort hits to be removed
    std::sort( hit_idx_to_be_removed_v.begin(), hit_idx_to_be_removed_v.end() );
    // current index in sorted vector where we stopped searching
    size_t current_idx = 0;
    size_t max_idx = hit_idx_to_be_removed_v.size() -1 ;
    // loop through all hits
    for (size_t hit_idx = 0; hit_idx < ev_hit->size(); hit_idx++){
      // if hit index is less than what is currently @ the vector index -> save
      if (hit_idx < hit_idx_to_be_removed_v[ current_idx ] ){
	out_hit->emplace_back( ev_hit->at( hit_idx ) );
	continue;
      }
      // if hit index is the same, remove and move to the next
      else if (hit_idx == hit_idx_to_be_removed_v[ current_idx ] ){
	if (current_idx < max_idx )
	  current_idx += 1;
	continue;
      }
      // if hit index is larger than current index, move counter forward and check again
      else if (hit_idx > hit_idx_to_be_removed_v[ current_idx ] ) {
	if (current_idx < max_idx )
	  current_idx += 1;
	if (hit_idx < hit_idx_to_be_removed_v[ current_idx ] ){
	  out_hit->emplace_back( ev_hit->at( hit_idx ) );
	  continue;
	}
	// if hit index is the same, remove and move to the next
	if (hit_idx == hit_idx_to_be_removed_v[ current_idx ] ){
	  if (current_idx < max_idx )
	    current_idx += 1;
	  continue;
	}
      }
    }// for all hits
    
    return true;
  }

  bool LinearClusterSubsetRemoval::finalize() {

    if (_fout){
      _fout->cd();
      _tree->Write();
    }
      

    return true;
  }

}
#endif
