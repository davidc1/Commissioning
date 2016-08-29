#ifndef LARLITE_TRKMIPCALIBRATION_CXX
#define LARLITE_TRKMIPCALIBRATION_CXX

#include "TrkMIPCalibration.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"
#include "GeoAlgo/GeoTrajectory.h"
#include "DataFormat/track.h"
#include "DataFormat/hit.h"

namespace larlite {

  TrkMIPCalibration::TrkMIPCalibration()
    : _tree(nullptr)
  {
    _name="TrkMIPCalibration";
    _trk_producer = "pandoraNu";
    _fout=0;
  }

  bool TrkMIPCalibration::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("tree","tree");
    _tree->Branch("_trk_len",&_trk_len,"_trk_len/D");
    _tree->Branch("_plane",&_plane,"_plane/I");
    _tree->Branch("_pitch",&_pitch,"_pitch/D");
    _tree->Branch("_amp_v","std::vector<double>",&_amp_v);
    _tree->Branch("_dqdx_v","std::vector<double>",&_dqdx_v);
    _tree->Branch("_dedx_v","std::vector<double>",&_dedx_v);

    double _e_to_eV   = 23.6;   // ionization energy of Ar in eV
    double _eV_to_MeV = 1e-6; // eV -> MeV conversion
    double MeV_to_fC  = 1. / ( _e_to_eV * _eV_to_MeV );
    _recomb_factor = larutil::LArProperties::GetME()->ModBoxInverse( 2.0 ) / ( 2.0 * MeV_to_fC );
    std::cout << "Recombination factor : " << _recomb_factor << std::endl;

    return true;
  }
  
  bool TrkMIPCalibration::analyze(storage_manager* storage) {

    auto _geoH = larutil::GeometryHelper::GetME();

    // start with input tracks from a certain producer.
    auto* ev_trk = storage->get_data<event_track>(_trk_producer);

    // are there any tracks found?
    if (!ev_trk or (ev_trk->size() == 0)){
      print(msg::kERROR, __FUNCTION__,Form("No trk w/ producer %s",_trk_producer.c_str() ) );
      return false;
    }

    // grab the hits associated to these tracks
    event_hit* ev_hit = nullptr;
    // associations
    AssSet_t ass_trk_hit_v;
    ass_trk_hit_v = storage->find_one_ass(ev_trk->id(), ev_hit, ev_trk->name());

    // are there any hits found?
    if (!ev_hit or (ev_hit->size() == 0)){
      print(msg::kERROR, __FUNCTION__,Form("No hit associated to trk w/ producer %s",_trk_producer.c_str() ) );
      return false;
    }

    //for (size_t i=0; i < ev_hit->size(); i++)
    //std::cout << "plane " << ev_hit->at(i).WireID().Plane << std::endl;
    
    // grab tracks and their associated hits.
    for (size_t t=0; t < ev_trk->size(); t++){

      // grab this track
      auto const& trk = ev_trk->at(t);

      // produce a GeoAlgo trajectory
      geoalgo::Trajectory_t trj;
      for (size_t pt = 0; pt < trk.NumberTrajectoryPoints(); pt++)
	trj.push_back( geoalgo::Vector_t( trk.LocationAtPoint(pt) ) );

      // get a 3D direction for the track
      auto const& Dir3D = Get3DDir(trj);
      TVector3 TDir3D(Dir3D[0], Dir3D[1], Dir3D[2]);

      //std::cout << "This track has 3D direction " << Dir3D << std::endl;

      // get hits associated with this track
      auto const& ass_hit_v = ass_trk_hit_v[t];

      // prepare vector of hits per plane
      std::vector< std::vector<larlite::hit> > hit_v_v(3, std::vector<larlite::hit>());

      // loop through hits, for each plane apply plane-by-plane calibration constants
      for (size_t h=0; h < ass_hit_v.size(); h++){
	auto const& hit = ev_hit->at( ass_hit_v[h] );
	int pl  = hit.WireID().Plane;
	//std::cout << "Plane : " << pl << "\tArea : " << hit.Integral() << std::endl;
	hit_v_v[pl].push_back(hit);
      }// for all associated hits
	
      // loop throuhg the hits on each plane, and fill calorimetry values
      for (size_t pl=0; pl < hit_v_v.size(); pl++){

	_amp_v.clear();
	_dqdx_v.clear();
	_dedx_v.clear();

	_plane = pl;

	// calculate the pitch on this plane
	_pitch = _geoH->GetPitch( TDir3D, pl );

	for (auto const& hit : hit_v_v[pl]){

	  double amp = hit.Integral();
	  _amp_v.push_back(amp);

	  double dQ = _caloAlg.ElectronsFromADCArea(amp, pl);
	  double dQdx = dQ / _pitch;
	  _dqdx_v.push_back( dQdx );

	  double dEdx = (dQdx * 23.6) / 1e6;
	  dEdx /= _recomb_factor;
	  dEdx *= 1.15; // lifetime
	  _dedx_v.push_back(dEdx);
	  
	}// for all hits on this plane
	
	// fill tree
	_tree->Fill();
	
      }

    }// for all reconstructed tracks
  
    return true;
  }

  bool TrkMIPCalibration::finalize() {

    _tree->Write();

    return true;
  }

  geoalgo::Vector_t TrkMIPCalibration::Get3DDir(const geoalgo::Trajectory_t& trj){

    geoalgo::Vector_t dir(0,0,0);

    // total length added
    double lenTot = 0.;

    // average the directions of all segments in the track.
    for (size_t pt = 0 ; pt < trj.size() - 1; pt++){
      
      geoalgo::Vector_t thisdir = trj[pt+1] - trj[pt];
      // get the length
      double len = thisdir.Length();
      lenTot += len;
      // now normalize
      thisdir.Normalize();

      // add to length-weighted average 3D direction
      dir  = dir + thisdir * len;

    }// for all segments in the trajectory

    _trk_len = lenTot;

    dir = dir / lenTot;

    return dir;
  }
  
}
#endif
