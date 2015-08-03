/**
 * \file RawHitMaker.h
 *
 * \ingroup RawHitMaker
 * 
 * \brief Class def header for a class RawHitMaker
 *
 * @author david
 */

/** \addtogroup RawHitMaker

    @{*/

#ifndef LARLITE_RAWHITMAKER_H
#define LARLITE_RAWHITMAKER_H

#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"

namespace larlite {
  /**
     \class RawHitMaker
     User custom analysis class made by SHELL_USER_NAME
   */
  class RawHitMaker : public ana_base{
  
  public:

    /// Default constructor
    RawHitMaker(){ _name="RawHitMaker"; _fout=0;}

    /// Default destructor
    virtual ~RawHitMaker(){}

    /** IMPLEMENT in RawHitMaker.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in RawHitMaker.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in RawHitMaker.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    // set RMS Max to cut on high noise channels
    void setRMSMax(double rms) { _rms_max = rms; }
    // set number of ticks to use to calculate baseline/rms
    void setNTicks(size_t n) { _nticks = n; }
    // set producer name for RawDigits
    void setRawDigitProducer(std::string s) { _producer = s; }
    // set cut on signal amplitude to have a hit (in sigmas based on rms)
    void setSigmaCut(double s) { _sigmacut = s; }

    /// function to calculate baseline
    double getBaseline(const std::vector<short> &wf);
    /// function to calculate RMS
    double getRMS(const std::vector<short> &wf, const double &baseline);
    /// function to find hits on waveform
    std::vector<larlite::hit> getHits(const std::vector<short> &wf,
				      const double &baseline,
				      const double &rms);

  protected:

    // producer name for RawDigits
    std::string _producer;
    // cut on noise level (only channels with noise less than this value will be searched for hits)
    double _rms_max;
    // number of ticks to use to calculate baseline & rms
    size_t _nticks;
    // sigma cut on the noise -> once waveform goes above/below this threshold defines the start/end of a hit
    double _sigmacut;
    
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
