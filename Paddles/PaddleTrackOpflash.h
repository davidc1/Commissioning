/**
 * \file PaddleTrackOpflash.h
 *
 * \ingroup Paddles
 * 
 * \brief Class def header for a class PaddleTrackOpflash
 *
 * @author Rui
 */

/** \addtogroup Paddles

    @{*/

#ifndef LARLITE_PADDLETRACKOPFLASH_H
#define LARLITE_PADDLETRACKOPFLASH_H

#include "Analysis/ana_base.h"
#include "DataFormat/track.h"
#include "DataFormat/ophit.h"
#include "DataFormat/opflash.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/calorimetry.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/Geometry.h"
#include "TTree.h"
#include "PaddleTrackAna.h"
#include "iostream"
#include "fstream"
#include <functional>
#include <algorithm>
#include <math.h>
#include "OpT0Finder/Base/FlashMatchManager.h"
#include "OpT0Finder/Algorithms/LightPath.h"
#include "OpT0Finder/App/MCQCluster.h"
#include "OpT0Finder/PhotonLibrary/PhotonVisibilityService.h"
#include "OpT0Finder/Algorithms/PhotonLibHypothesis.h"
#include <numeric>

namespace larlite {
  /**
     \class PaddleTrackOpflash
     User custom analysis class made by SHELL_USER_NAME
   */
  class PaddleTrackOpflash : public ana_base{
  
  public:

    /// Default constructor
    PaddleTrackOpflash()
      :_tree(nullptr)
      { _name="PaddleTrackOpflash";
	_pe_dis_hist = 0;
	_fout=0;
	_fout=0;
      }

    /// Default destructor
    virtual ~PaddleTrackOpflash(){}

    /** IMPLEMENT in PaddleTrackOpflash.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in PaddleTrackOpflash.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in PaddleTrackOpflash.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    //setter function
    bool UseData       (bool use) {_useData       = use ; _useSimulation = !use; return  _useData;}
    bool UseQCluster   (bool use) {_useQCluster   = use ; _useMCQCluster = !use ;return  _useQCluster;}
    
  protected:
    
    TTree* _tree;
    TH1F* _pe_dis_hist;
   
    std::ofstream _track_positions;
    
    bool _save_histos;
    bool _useData;
    bool _useSimulation;
    bool _useMCQCluster;
    bool _useQCluster;

    size_t _n_ev_reco;
    
    int _run;
    int _subrun;
    int _event;
    int _trk_id;
    
    int _test=0;
    int _n_evt;
    int _n_evt_paddle;
    int _n_evt_mc;
    int _n_tracks;
    int _opchannel_id;
    int _n_intersections_FV;
    int _n_intersections_mucs_top;
    int _n_intersections_mucs_bottom;
    
    std::vector<double> _t_opflash;
    std::vector<double> _t_ophit;
    std::vector<double> _pe_ophit;
    std::vector<double> _pe_mchit;
        
    //muon intersection w/ MuCS
    double _MuCS_ints_x_top;
    double _MuCS_ints_z_top;
    double _MuCS_ints_x_bottom;
    double _MuCS_ints_z_bottom;
    
    double _length_xfiducial;
    double _length_yfiducial;
    double _length_zfiducial;
    double _length_trj_prj_fv;
    double _length_trj_prj_fv_neg;
    
    double _theta;
    double _pe_mchit_sum;
    double _pe_ophit_sum;
    double _qratio_pl;
    double _qratio_re;
    
    double _mc_e;
    double _mc_e_dep;
    
    ::geoalgo::GeoAlgo _geoAlgo;
    ::geoalgo::AABox _vfiducial;
    ::geoalgo::AABox _vmucs_top;
    ::geoalgo::AABox _vmucs_bottom;
    std::vector<::geoalgo::Trajectory> _trj;
    //tracks contained in TPC
    std::vector<::geoalgo::Trajectory> _trj_con;
    //tracks pass MuCS
    std::vector<::geoalgo::Trajectory> _trj_mucs;
    std::vector<::geoalgo::HalfLine> _trj_prj;
    std::vector<::geoalgo::LineSegment> _prj_lineseg;

    std::vector<TVector3> _trj_filt;
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
