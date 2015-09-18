/**
 * \file SaveWFs.h
 *
 * \ingroup MakeWFTTree
 * 
 * \brief Class def header for a class SaveWFs
 *
 * @author david
 */

/** \addtogroup MakeWFTTree

    @{*/

#ifndef LARLITE_SAVEWFS_H
#define LARLITE_SAVEWFS_H

#include "Analysis/ana_base.h"
#include "DataFormat/rawdigit.h"
#include "TTree.h"
#include "LArUtil/Geometry.h"

namespace larlite {
  /**
     \class SaveWFs
     User custom analysis class made by SHELL_USER_NAME
   */
  class SaveWFs : public ana_base{
  
  public:

    /// Default constructor
    SaveWFs();

    /// Default destructor
    virtual ~SaveWFs(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    /// Function to be called to add a channel to the list
    /// of channels for which a WF should be saved
    void AddChannel(int ch) { _chList.push_back(ch); }

  protected:

    /// List of channels for which to save waveform
    std::vector<int> _chList;

    /// TTree on which to store info
    TTree* _tree;
    int    _ch;
    int    _pl;
    int    _evt;
    int    _wire;
    std::vector<short> _ADCs;
    
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
