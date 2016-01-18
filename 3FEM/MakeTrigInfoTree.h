/**
 * \file MakeTrigInfoTree.h
 *
 * \ingroup 3FEM
 * 
 * \brief Class def header for a class MakeTrigInfoTree
 *
 * @author david
 */

/** \addtogroup 3FEM

    @{*/

#ifndef LARLITE_MAKETRIGINFOTREE_H
#define LARLITE_MAKETRIGINFOTREE_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class MakeTrigInfoTree
     User custom analysis class made by SHELL_USER_NAME
   */
  class MakeTrigInfoTree : public ana_base{
  
  public:

    /// Default constructor
    MakeTrigInfoTree();

    /// Default destructor
    virtual ~MakeTrigInfoTree(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:

    TTree* _tree;
    double _trig_time;
    std::vector<unsigned short> _ch0, _ch1;
    int _event, _event_frame_num, _fem_trig_sample_number_RAW;
    
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
