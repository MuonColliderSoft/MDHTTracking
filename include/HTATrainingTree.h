#pragma once

#include <vector>

#include "marlin/Processor.h"

#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

#include "TTree.h"

/**
 * Produces a flat ntuple for the HT array training
 *
 * @parameter MCParticleCollection Name of the Monte Carlo particle collection
 * @parameter TrackerHitCollections List of the tracker hits collections
 * @parameter TrackerSimHitCollections List of the tracker simhits collections
 * @parameter TrackerHitRelCollections List of the tracker hit relations collections
 * @parameter ParticleTypes Particle types to be saved in the tree
 *
 * @author M. Casarsa
 * @date  21 March 2025
 * @version $Id: $
 */
class HTATrainingTree : public marlin::Processor
{
public:
  virtual Processor*  newProcessor() { return new HTATrainingTree ; }

  HTATrainingTree(const HTATrainingTree &) = delete ;
  HTATrainingTree& operator =(const HTATrainingTree &) = delete ;
  HTATrainingTree() ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 

  virtual void check( LCEvent * evt ) ; 

  /** Called after data processing for clean up.
   */
  virtual void end() ;  

private:

  const dd4hep::rec::SurfaceMap* _surfMap {nullptr};

  TTree* _HTAtree {nullptr};
  
  int _nRun {0};
  int _nEvt {0};

  int   _part_pdg {0};
  float _part_px {0.};
  float _part_py {0.};
  float _part_pz {0.};
  float _part_e {0.};
  float _part_vx {0.};
  float _part_vy {0.};
  float _part_vz {0.};
  float _part_t0 {0.};
  float _part_q {0.};

  static const int MAX_HITS = 1000;
  
  int _n_hit {0};
  int _hit_index[MAX_HITS] = {0};
  int _hit_mcp[MAX_HITS] = {0};
  int _hit_id0[MAX_HITS] = {0};
  float _hit_x[MAX_HITS] = {0.0f};
  float _hit_y[MAX_HITS] = {0.0f};
  float _hit_z[MAX_HITS] = {0.0f};
  float _hit_t[MAX_HITS] = {0.0f};
  float _hit_xloc[MAX_HITS] = {0.0f};
  float _hit_yloc[MAX_HITS] = {0.0f};
  //std::vector<int>   _hit_index {};
  //std::vector<int>   _hit_mcp {};
  //std::vector<int>   _hit_id0 {};
  //std::vector<float> _hit_x {};
  //std::vector<float> _hit_y {};
  //std::vector<float> _hit_z {};
  //std::vector<float> _hit_t {};
  //std::vector<float> _hit_xloc {};
  //std::vector<float> _hit_yloc {};

  //! Input MC particle collection
  std::string _inputMCParticleCollection {};
  
  //! Input tracker hit collections
  std::vector<std::string> _inputTrackerHitCollections {};

  //! Input tracker simhit collections
  std::vector<std::string> _inputTrackerSimHitCollections {};

  //! Input tracker hit relation collections
  std::vector<std::string> _inputTrackerHitRelCollections {};

  //! Particle types
  std::vector<int> _particleTypes {};

  bool _saveOnlyPartHits = true;

};
