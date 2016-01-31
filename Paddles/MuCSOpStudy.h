/**
 * \file MuCSOpStudy.h
 *
 * \ingroup Paddles
 * 
 * \brief Class def header for a class MuCSOpStudy
 *
 * @author kazuhiro
 */

/** \addtogroup Paddles

    @{*/

#ifndef LARLITE_MUCSOPSTUDY_H
#define LARLITE_MUCSOPSTUDY_H

#include "Analysis/ana_base.h"
#include "OpT0Finder/Base/FlashMatchManager.h"
#include "OpT0Finder/Algorithms/LightPath.h"
#include "OpT0Finder/Algorithms/PhotonLibHypothesis.h"
#include "OpT0Finder/Algorithms/QLLMatch.h"
#include <TH1D.h>
#include <TH2D.h>
namespace larlite {
  /**
     \class MuCSOpStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class MuCSOpStudy : public ana_base{
  
  public:

    /// Default constructor
    MuCSOpStudy();

    /// Default destructor
    virtual ~MuCSOpStudy(){}

    /** IMPLEMENT in MuCSOpStudy.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MuCSOpStudy.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MuCSOpStudy.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void configure(const std::string config_file);

    ::flashana::FlashMatchManager& FlashMatchManager() { return _mgr; }

  protected:

    std::vector<geoalgo::Trajectory>  _cand_trj_v;
    std::vector<flashana::QCluster_t> _qcluster_v;
    std::vector<flashana::Flash_t> _flash_v;
    ::flashana::Flash_t _ophit_flash;

    TH2D* _hRatioMap;
    TH1D* _hHitFlashScore;
    TH1D* _hMatchTime;
    TH1D* _hMatchScore;
    TH2D* _hMatchScorePE;
    TH2D* _hMatchScoreTime;
    
    ::flashana::FlashMatchManager _mgr;
    ::flashana::LightPath _lpath;
    ::flashana::PhotonLibHypothesis _fhypo;
    ::flashana::QLLMatch _qll;

    double _ophit_tmin;
    double _ophit_tmax;

    std::string _ophit_producer;
    std::string _opflash_producer;
    std::string _track_producer;
    std::string _ctag_producer;

    bool _run_match;
    bool _use_ophit_flash;
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
