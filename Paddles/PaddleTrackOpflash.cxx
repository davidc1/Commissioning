#ifndef LARLITE_PADDLETRACKOPFLASH_CXX
#define LARLITE_PADDLETRACKOPFLASH_CXX

#include "PaddleTrackOpflash.h"

namespace larlite {
  
  bool PaddleTrackOpflash::initialize() {
    
    if(!_pe_dis_hist)
      _pe_dis_hist =  new TH1F("PE distribution","PE distribution",32,0,32);
    
    if (_tree) {delete _tree;}
    _tree = new TTree("PaddleTree", "PaddleTree");
    
    _tree->Branch("run",&_run,"_run/I");
    _tree->Branch("subrun",&_subrun,"_subrun/I");
    _tree->Branch("event",&_event,"_event/I");
    _tree->Branch("trk_id",&_trk_id,"_trk_id/I");
    _tree->Branch("zenith",&_theta,"_theta/D");
    _tree->Branch("QRatio_pl",&_qratio_pl,"_qratio_pl/D"); 
    _tree->Branch("QRatio_re",&_qratio_re,"_qratio_re/D");
    _tree->Branch("pe_mchit_sum",&_pe_mchit_sum,"_pe_mchit_sum/D");
    _tree->Branch("pe_ophit_sum",&_pe_ophit_sum,"_pe_ophit_sum/D");
    
    _tree->Branch("trj_filt","std::vector<TVector3>",&_trj_filt);
                                      
    _tree->Branch("t_opflash","std::vector<double>",&_t_opflash);
    _tree->Branch("t_ophit","std::vector<double>",&_t_ophit);
    _tree->Branch("pe_ophit","std::vector<double>",&_pe_ophit);
    _tree->Branch("pe_mchit","std::vector<double>",&_pe_mchit);
    
    _length_xfiducial = larutil::Geometry::GetME()->DetHalfWidth();
    _length_yfiducial = larutil::Geometry::GetME()->DetHalfHeight();
    _length_zfiducial = larutil::Geometry::GetME()->DetLength();

    _vfiducial = ::geoalgo::AABox(0, -_length_yfiducial, 0, 2 * _length_xfiducial, _length_yfiducial,_length_zfiducial);
    _vmucs_top = ::geoalgo::AABox(-71.795, 393.941, 531.45, -23.795, 398.451, 579.45);
    //_vmucs_top = ::geoalgo::AABox(-271.795, 393.941, 331.45, 223.795, 398.451, 779.45);
    _vmucs_bottom = ::geoalgo::AABox(-19.6948, 316.041, 533.25, 28.3052, 320.551, 581.25);
    
    _track_positions.open ("TrackPositions.txt");
    
    _n_evt = 0;
    
    //fPhotonLib = new ubphotonlib::PhotonLibrary();
    //fPhotonLib->LoadDefault();
    
    return true;
  }
  
  bool PaddleTrackOpflash::analyze(storage_manager* storage) {

    
    auto const geo = ::larutil::Geometry::GetME();
    
    _pe_ophit.resize(32,0.0);
    _pe_mchit.resize(32,0.0);
    _t_opflash.clear();
    _t_ophit.clear();
    std::fill(_pe_ophit.begin(), _pe_ophit.end(), 0);
    
    auto ev_reco = storage->get_data<event_track>("trackkalmanhit");
    if (!ev_reco) {
      std::cout<<"........Couldn't find reco track data product in this event...... "<<std::endl;
    }

    auto ev_ophit = storage->get_data<event_ophit>("OpHitFinder");
    if (!ev_ophit) {
      std::cout<<"........Couldn't find ophit data product in this event...... "<<std::endl;
    }
    /*
    auto ev_opflash = storage->get_data<event_opflash>("opFlash");
    if (!ev_opflash) {
      std::cout<<"........Couldn't find opflash data product in this event...... "<<std::endl;
    }
    */
    
    for(size_t i = 0; i <ev_reco->size(); i++ ){
      
      auto const& trk = ev_reco->at(i);
      
      if (trk.NumberTrajectoryPoints() > 1) {
	
        ::geoalgo::Trajectory trj;

        for (size_t pt = 0; pt < trk.NumberTrajectoryPoints(); pt++) {
	  
          auto const& pos = trk.LocationAtPoint(pt);
	  
	  trj.push_back(::geoalgo::Vector(pos[0], pos[1], pos[2]));
	  
	  
	}//Construct trajectory using points on track
	
	::geoalgo::HalfLine trj_prj(trj[0], trj[0]-trj[trj.size()/2]);
        ::geoalgo::HalfLine trj_prj_neg(trj[trj.size()-1], trj[trj.size()-1]-trj[trj.size()/2]);
	
	auto const&  intersection_trj_prj_top    = _geoAlgo.Intersection(_vmucs_top,trj_prj);
        auto const&  intersection_trj_prj_bottom = _geoAlgo.Intersection(_vmucs_bottom,trj_prj);
	
	//auto const&  intersection_trj_prj_top_neg    = _geoAlgo.Intersection(_vmucs_top,trj_prj_neg);
        //auto const&  intersection_trj_prj_bottom_neg = _geoAlgo.Intersection(_vmucs_bottom,trj_prj_neg);
	
	auto const&  intersection_trj_prj_fv        = _geoAlgo.Intersection(_vfiducial,trj_prj);
	auto const&  intersection_trj_prj_fv_neg    = _geoAlgo.Intersection(_vfiducial,trj_prj_neg);
	//auto const&  intersection_trj_fv            = _geoAlgo.Intersection(_vfiducial,trj);
	
	if(intersection_trj_prj_top.size()>0 && intersection_trj_prj_bottom.size()>0){
	  
	  //_n_evt++;
	  
	  std::vector<double> pt1,pt2,pt3,delta_p12,delta_p13,delta_p23;//1,2 for intersections on mucs_top/bottom, 3 for 1 project to mucs_bottom 
	  pt1.resize(3,0.0);
	  pt2.resize(3,0.0);
	  pt3.resize(3,0.0);
	  delta_p12.resize(3,0.0);
	  delta_p13.resize(3,0.0);
	  delta_p23.resize(3,0.0);
	  pt1 = intersection_trj_prj_top[0];
	  pt2 = intersection_trj_prj_bottom[0];
	  pt3 = {intersection_trj_prj_top[0].at(0),intersection_trj_prj_bottom[0].at(1),intersection_trj_prj_top[0].at(2)};
	  
	  for (size_t i = 0; i<3 ; i++){
	    delta_p12.push_back(pt1.at(i)-pt2.at(i));
	    delta_p13.push_back(pt1.at(i)-pt3.at(i));
	    delta_p23.push_back(pt2.at(i)-pt3.at(i));
	  }
	  double length_12 = sqrt(std::inner_product(begin(delta_p12), end(delta_p12), begin(delta_p12), 0.0));
	  double length_13 = sqrt(std::inner_product(begin(delta_p13), end(delta_p13), begin(delta_p13), 0.0));
	  double length_23 = sqrt(std::inner_product(begin(delta_p23), end(delta_p23), begin(delta_p23), 0.0));
	  
	  double cos_theta = (pow(length_12,2)+pow(length_13,2)-pow(length_23,2))/(2.*length_12*length_13); 
	  
	  if(length_12*length_13>0)_theta = acos (cos_theta) * 180.0 / M_PI;
	  
	  _MuCS_ints_x_top = intersection_trj_prj_top.at(0).at(0);
          _MuCS_ints_z_top = intersection_trj_prj_top.at(0).at(2);
          _MuCS_ints_x_bottom =  intersection_trj_prj_bottom.at(0).at(0);
          _MuCS_ints_z_bottom =  intersection_trj_prj_bottom.at(0).at(2);
	  
	  {
	    /*
	      for(size_t opf = 0; opf < ev_opflash->size(); opf++){
	      auto const& opflash = ev_opflash->at(opf);
	      _t_opflash.push_back(opflash.Time());
	      
	    }
	    */
	  }//loop over all opflashes
	  
	  {
	    for(size_t oph = 0; oph < ev_ophit->size(); oph++){
              auto const& ophit = ev_ophit->at(oph);
              _t_ophit.push_back(ophit.PeakTime());
	      //_pe_ophit.push_back(ophit.PE());
	      if(ophit.PeakTime()<-0.8 &&ophit.PeakTime()>-1){
		auto const pmt_id = geo->OpDetFromOpChannel(ophit.OpChannel());
		_pe_ophit[pmt_id] += ophit.PE();
	      }
	    }
	  }//loop over all ophit
	  
	  _run    = storage->get_data<event_track>("trackkalmanhit")->run();
          _subrun = storage->get_data<event_track>("trackkalmanhit")->subrun();
          _event  = storage->get_data<event_track>("trackkalmanhit")->event_id();
          _trk_id = trk.ID();
	  
	  //Save positions into txt file
	  _track_positions<<"Run: "<<_run<<" , Subrun: "<<_subrun<<" , Event: "<<_event<<" , Track_ID: "<<_trk_id<<std::endl; 
	  //_track_positions<<trj.front()<<std::endl;
	  //_track_positions<<trj.back()<<std::endl;

	  //
	  // Implement PhotonVisibility from OpT0Finder now
	  //
	  
	  ::flashana::QCluster_t tpc_obj;
	  ::flashana::LightPath LP;
	  ::flashana::Flash_t flash_obj;
	  ::flashana::PhotonLibHypothesis PL;
	  
	  tpc_obj = LP.FlashHypothesis(trj);
	  /*
	  for(size_t i=0;i<tpc_obj.size();i++){
	    std::cout<<"wtf?"<<tpc_obj.at(i).q<<std::endl;
	    }*/
	  flash_obj.pe_v.resize(32,0.0);
	  PL.FillEstimate(tpc_obj,flash_obj);
	  _pe_mchit =  flash_obj.pe_v;
	  //for(size_t i = 0; i<32;i++)std::cout<<flash_obj.pe_v.at(i);
	  	  
	  // Normalize
	  double pe_mchit_sum = 0;
	  double pe_ophit_sum = 0;
	  
	  pe_mchit_sum = std::accumulate(std::begin(_pe_mchit),std::end(_pe_mchit),0.0);
	  pe_ophit_sum = std::accumulate(std::begin(_pe_ophit),std::end(_pe_ophit),0.0);
	  
	  _pe_mchit_sum = pe_mchit_sum;
	  _pe_ophit_sum = pe_ophit_sum;
	  
	  /*
	  if(pe_ophit_sum<1) {
	    // Make rough guess on overall PE scale
	    pe_ophit_sum = 29000.0 * 2.3 * trj.Length();
	  }
	  for(auto& v : _pe_mchit) v *= ( pe_ophit_sum / pe_mchit_sum );
	  //std::cout<<trk.NumberTrajectoryPoints()<<std::endl;
	  for (size_t pt = 0; pt < trk.NumberTrajectoryPoints(); pt++) {
	    auto const& pos = trk.LocationAtPoint(pt);
	    _trj_filt.push_back(pos);
	  }
	  */
	  //Get Q ratio
	  std::sort(_pe_ophit.begin(),_pe_ophit.end());
	  std::sort(_pe_mchit.begin(),_pe_mchit.end());
	  double QRatio_pl = _pe_mchit.at(31)/pe_mchit_sum;//photon library
	  double QRatio_re = _pe_ophit.at(31)/pe_ophit_sum;//ophit
	  _qratio_pl = QRatio_pl;
	  _qratio_re = QRatio_re;
	  
	  /////////////////further select events
	  
	  {
            double length = trj.Length();
            _length_trj_prj_fv = 0;
            _length_trj_prj_fv_neg = 0;

            if(length>150){//1ST, track length larger than 50cm
	      std::vector<double> pt1,pt2,delta_p12;
              pt1.resize(3,0.0);
              pt2.resize(3,0.0);
              delta_p12.resize(3,0.0);
              if(intersection_trj_prj_fv.size()>0){//2ND, sum of gaps btw start/end of track and fv smaller than
                pt1 = intersection_trj_prj_fv[0];
                for(size_t i = 0; i<3 ; i++){
                  pt2[i] = trj.at(0).at(i);
                  delta_p12.push_back(pt1.at(i)-pt2.at(i));
                }
                _length_trj_prj_fv = sqrt(std::inner_product(begin(delta_p12), end(delta_p12), begin(delta_p12), 0.0));
              }
              pt1.resize(3,0.0);
              pt2.resize(3,0.0);
              delta_p12.resize(3,0.0);
              if(intersection_trj_prj_fv_neg.size()>0){
                pt1 = intersection_trj_prj_fv_neg[0];
                for(size_t i = 0; i<3 ; i++){
                  pt2[i] = trj.at(trj.size()-1).at(i);
                  delta_p12.push_back(pt1.at(i)-pt2.at(i));
                }
                _length_trj_prj_fv_neg = sqrt(std::inner_product(begin(delta_p12), end(delta_p12), begin(delta_p12), 0.0));
              }
              if(_length_trj_prj_fv+_length_trj_prj_fv_neg<100){
                _n_evt++;
		//std::cout<<_run<<" "<<_subrun<<" "<<_event<<std::endl;
		_tree->Fill();
	      }
	    }
          }
	  //_tree->Fill();
        }
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
