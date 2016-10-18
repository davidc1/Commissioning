#ifndef LARLITE_LINEARCLUSTERREMOVAL_CXX
#define LARLITE_LINEARCLUSTERREMOVAL_CXX

#include "LinearClusterRemoval.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"

namespace larlite {

  LinearClusterRemoval::LinearClusterRemoval()
    : _tree(nullptr)
  {

    _name        = "LinearClusterRemoval";
    _fout        = 0;
    _verbose     = false;
    _clusterProducer = "gaushit";
    _out_clusterProducer = "shrcluster";
    _out_hitProducer     = "shrhits";
    _max_lin_v = {0.0};
    _min_n_hits_v = {0};

  }

  bool LinearClusterRemoval::initialize() {

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
  
  bool LinearClusterRemoval::analyze(storage_manager* storage) {

    auto ev_clus = storage->get_data<event_cluster>(_clusterProducer);

    larlite::event_hit* ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());
    
    auto out_hit = storage->get_data<event_hit>( _out_hitProducer );
    auto out_clusters   = storage->get_data<event_cluster>( _out_clusterProducer );
    auto out_cluster_ass_v = storage->get_data<event_ass>(out_clusters->name());
    std::vector<std::vector<unsigned int> > out_cluster_hit_ass_v;
    
    //set event ID through storage manager
    storage->set_id(storage->run_id(),storage->subrun_id(),storage->event_id());

    if (!ev_hit){
      std::cout << "No hits!" << std::endl;
      return false;
    }

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
      
      // get coordinates of hits to calculate linearity
      std::vector<double> hit_w_v;
      std::vector<double> hit_t_v;
      
      for (auto const& hit_idx : hit_idx_v){
	hit_w_v.push_back( ev_hit->at(hit_idx).WireID().Wire  * _wire2cm );
	hit_t_v.push_back( ev_hit->at(hit_idx).PeakTime() * _time2cm );
      }

      ::Linearity lin(hit_w_v,hit_t_v);

      _nhits = hit_w_v.size();
      _lin                 = lin._lin;
      _local_lin_avg       = lin._local_lin_avg;
      _local_lin_truncated = lin._local_lin_truncated_avg;

      _tree->Fill();

      if (_local_lin_truncated < max_lin)
	remove = true;

      /* OLD WAY
      // calculate covariance
      double L = linearity(hit_w_v,hit_t_v);

      if ( log(fabs(L)) < max_lin)
	remove = true;
      */ 
      
      if (remove == false){
	// for all hits, add them to output
	for (auto const& hit_idx : hit_idx_v){
	  out_hit->emplace_back( ev_hit->at( hit_idx ) );
	  out_cluster_hit_idx_v.push_back( out_hit->size() - 1);
	}
      }// if hits are not to be removed

      if ( out_cluster_hit_idx_v.size() >0 ){
	out_cluster_hit_ass_v.push_back( out_cluster_hit_idx_v );
      	out_clus.set_n_hits( out_cluster_hit_idx_v.size() );
	out_clus.set_view(larlite::geo::View_t::kW);
	out_clusters->emplace_back( out_clus );
      }
      
    }// loop through all planes

    out_cluster_ass_v->set_association(out_clusters->id(),product_id(data::kHit,out_hit->name()),out_cluster_hit_ass_v);    
    
    return true;
  }

  bool LinearClusterRemoval::finalize() {

    if (_fout){
      _fout->cd();
      _tree->Write();
    }
      

    return true;
  }

  // covariance
  double LinearClusterRemoval::cov (const std::vector<double>& data1,
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
  
  double LinearClusterRemoval::stdev(const std::vector<double>& data) const
  {
    if(data.size() == 0) { std::cout << "zero-vector!" << std::endl; return 0; }

    double result = 0.0;
    auto    avg   = mean(data);
    for(const auto& d: data)
      result += (d - avg)*(d - avg);
    
    return sqrt(result/((double)data.size()));
  }
  
  double LinearClusterRemoval::mean(const std::vector<double>& data) const
  {
    if(data.size() == 0) { std::cout << "zero-vector!" << std::endl; return 0; }
	
    double result = 0.0;

    for(const auto& d : data) 
      result += d;
        
    return (result / ((double)data.size()));
  }

  double LinearClusterRemoval::linearity(const std::vector<double>& data1,
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
  
}
#endif
