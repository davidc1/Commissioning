#ifndef LARLITE_PADDLETRACKOPFLASH_CXX
#define LARLITE_PADDLETRACKOPFLASH_CXX

#include "PaddleTrackOpflash.h"
#include "OpT0Finder/OpTrackAlg/TrackHypothesis.h"
#include "OpT0Finder/OpTrackAlg/OpTrackModule.h"

namespace larlite {
  
  bool PaddleTrackOpflash::initialize() {
    
    if(!_pe_dis_hist)
      _pe_dis_hist =  new TH1F("PE distribution","PE distribution",32,0,32);
    
    if (_tree) {delete _tree;}
    _tree = new TTree("PaddleTree", "PaddleTree");
    _tree->Branch("t_opflash","std::vector<double>",&_t_opflash);
    _tree->Branch("t_ophit","std::vector<double>",&_t_ophit);
    _tree->Branch("pe_ophit","std::vector<double>",&_pe_ophit);
    _tree->Branch("pe_mchit","std::vector<double>",&_pe_mchit);
    
    _tree->Branch("run",&_run,"_run/I");
    _tree->Branch("subrun",&_subrun,"_subrun/I");
    _tree->Branch("event",&_event,"_event/I");
    _tree->Branch("trk_id",&_trk_id,"_trk_id/I");
    
    _tree->Branch("trj_filt","std::vector<TVector3>",&_trj_filt);
                                      
    _length_xfiducial = larutil::Geometry::GetME()->DetHalfWidth();
    _length_yfiducial = larutil::Geometry::GetME()->DetHalfHeight();
    _length_zfiducial = larutil::Geometry::GetME()->DetLength();

    _vfiducial = ::geoalgo::AABox(0, -_length_yfiducial, 0, 2 * _length_xfiducial, _length_yfiducial,_length_zfiducial);
    _vmucs_top = ::geoalgo::AABox(-71.795, 393.941, 531.45, -23.795, 398.451, 579.45);
    _vmucs_bottom = ::geoalgo::AABox(-19.6948, 316.041, 533.25, 28.3052, 320.551, 581.25);

    _track_positions.open ("TrackPositions.txt");
    
    _n_evt = 0;
    
    _pe_ophit.resize(32);

    fPhotonLib = new ubphotonlib::PhotonLibrary();
    fPhotonLib->LoadDefault();
    
    return true;
  }
  
  bool PaddleTrackOpflash::analyze(storage_manager* storage) {
    
    _t_opflash.clear();
    _t_ophit.clear();
    std::fill(_pe_ophit.begin(), _pe_ophit.end(), 0);
    
    auto ev_reco = storage->get_data<event_track>("trackkalmanhit");
    if (!ev_reco) {
      std::cout<<"........Couldn't find reco track data product in this event...... "<<std::endl;
    }

    auto ev_ophit = storage->get_data<event_ophit>("opFlash");
    if (!ev_ophit) {
      std::cout<<"........Couldn't find ophit data product in this event...... "<<std::endl;
    }
    
    auto ev_opflash = storage->get_data<event_opflash>("opFlash");
    if (!ev_opflash) {
      std::cout<<"........Couldn't find opflash data product in this event...... "<<std::endl;
    }
    
    for(size_t i = 0; i <ev_reco->size(); i++ ){
      
      auto const& trk = ev_reco->at(i);

      if (trk.NumberTrajectoryPoints() > 1) {

        //a geoalgo::Trajectory            
        ::geoalgo::Trajectory trj;

        for (size_t pt = 0; pt < trk.NumberTrajectoryPoints(); pt++) {
	  
          auto const& pos = trk.LocationAtPoint(pt);
	  
	  trj.push_back(::geoalgo::Vector(pos[0], pos[1], pos[2]));
	  
	  
	}// for all points in track

	::geoalgo::HalfLine trj_prj(trj[0], trj[0]-trj[trj.size()/2]);
        ::geoalgo::HalfLine trj_prj_neg(trj[trj.size()-1], trj[trj.size()/2]-trj[trj.size()-1]);

        auto const&  intersection_trj_prj_top    = _geoAlgo.Intersection(_vmucs_top,trj_prj);
        auto const&  intersection_trj_prj_bottom = _geoAlgo.Intersection(_vmucs_bottom,trj_prj);
	
        auto const&  intersection_trj_prj_top_neg    = _geoAlgo.Intersection(_vmucs_top,trj_prj_neg);
        auto const&  intersection_trj_prj_bottom_neg = _geoAlgo.Intersection(_vmucs_bottom,trj_prj_neg);
	
	if(intersection_trj_prj_top.size()>0 && intersection_trj_prj_bottom.size()>0){
          _n_evt++;
	  //std::cout<<intersection_trj_prj_top.size()<<",";
	  _MuCS_ints_x_top = intersection_trj_prj_top.at(0).at(0);
          _MuCS_ints_z_top = intersection_trj_prj_top.at(0).at(2);
          _MuCS_ints_x_bottom =  intersection_trj_prj_bottom.at(0).at(0);
          _MuCS_ints_z_bottom =  intersection_trj_prj_bottom.at(0).at(2);
	  
	  {
	    for(size_t opf = 0; opf < ev_opflash->size(); opf++){
	      auto const& opflash = ev_opflash->at(opf);
	      _t_opflash.push_back(opflash.Time());
	      
	    }
	  }//for all opflashes
	  
	  {
	    //std::cout<<ev_ophit->size()<<std::endl;
            for(size_t oph = 0; oph < ev_ophit->size(); oph++){
              auto const& ophit = ev_ophit->at(oph);
              _t_ophit.push_back(ophit.PeakTime());
	      //_pe_ophit.push_back(ophit.PE());
	      if(ophit.PeakTime()<-0.8 &&ophit.PeakTime()>-1){
		_pe_ophit.at(ophit.OpChannel()) = _pe_ophit.at(ophit.OpChannel())+ophit.PE();
	      }
	      //if (ophit.OpChannel()>30)std::cout<<ophit.OpChannel()<<std::endl;
	      
	    }
	    //std::cout<<_pe_ophit.at(31)<<std::endl;
          }//for all ophit
	  
	  _run    = storage->get_data<event_track>("trackkalmanhit")->run();
          _subrun = storage->get_data<event_track>("trackkalmanhit")->subrun();
          _event  = storage->get_data<event_track>("trackkalmanhit")->event_id();
          _trk_id = trk.ID();
	  //Save positions into txt file
	  _track_positions<<"Run: "<<_run<<" , Subrun: "<<_subrun<<" , Event: "<<_event<<" , Track_ID: "<<_trk_id<<std::endl; 
	  _track_positions<<trj.front()<<std::endl;
	  _track_positions<<trj.back()<<std::endl;
	  
	  /*
	  optrackfit::TrackHypothesis hypo( trj.front(), trj.back() );
	  std::vector< std::vector<double> > midpoints;
	  std::vector< std::vector<float> > chphotons;
	  optrackfit::makeVoxelList( hypo, *fPhotonLib, 29000.0, 2.3, 32, midpoints, chphotons );
	  */
	  
	  std::vector< std::vector<double> > midpoints;
	  std::vector< std::vector<float> > chphotons;

	  for( size_t step=0; step < trj.size()-1; step++){
	    optrackfit::TrackHypothesis hypo( trj.at(step), trj.at(step+1) );
	    optrackfit::makeVoxelList( hypo, *fPhotonLib, 29000.0, 2.3, 32, midpoints, chphotons );
	  }
	  
	  _pe_mchit.clear();
	  _pe_mchit.resize( 32, 0.0 );
	  for (int ich=0; ich<32; ich++) {
	    double chpe = 0.0;
	    for ( int istep=0; istep<(int)chphotons.size(); istep++ ) {
	      chpe += chphotons.at(istep).at(ich)*0.0098*0.3;
	    }
	    if ( chpe>5.0 )
	      _pe_mchit.at(ich) += chpe;
	  }
	  
	  std::cout<<trk.NumberTrajectoryPoints()<<std::endl;
	  
	  for (size_t pt = 0; pt < trk.NumberTrajectoryPoints(); pt++) {

	    auto const& pos = trk.LocationAtPoint(pt);
	    _trj_filt.push_back(pos);
	  
	  }
	  
	  
	  _tree->Fill();
        }
        //correct for wrong recon direction
	/*
	if(intersection_trj_prj_top_neg.size()>0 && intersection_trj_prj_bottom_neg.size()>0){
          _n_evt++;
	  //std::cout<<intersection_trj_prj_top_neg.size();
          _MuCS_ints_x_top = intersection_trj_prj_top_neg.at(0).at(0);
          _MuCS_ints_z_top = intersection_trj_prj_top_neg.at(0).at(2);
          _MuCS_ints_x_bottom =  intersection_trj_prj_bottom_neg.at(0).at(0);
          _MuCS_ints_z_bottom =  intersection_trj_prj_bottom_neg.at(0).at(2);
	  _tree->Fill();
        }
	*/
      }
    }///for all tracks
   
    return true;
  }
  
  bool PaddleTrackOpflash::finalize() {
    
    /*std::cout<<_pe_ophit.size()<<std::endl;
    for(size_t i=0; i< _pe_ophit.size(); i++){
      double bincontent = _pe_ophit.at(i);
      std::cout<<_pe_ophit.at(i)<<std::endl;
      _pe_dis_hist->SetBinContent(i+1,bincontent);
      }*/
    
    _track_positions.close();
    
    if(_fout){
      _fout->cd();
      if(_tree) _tree->Write();
      std::cout<<"\nThrough-going muon found in "<<_n_evt<<" events!!!???"<<std::endl;
    }
    
    if(_save_histos){
      _pe_dis_hist->Write();
    }
    
    return true;
  }

}
#endif
