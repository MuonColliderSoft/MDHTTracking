#pragma once

#include "marlin/Processor.h"

#include "UTIL/BitField64.h"

#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

#include "include/HTArray.h"

namespace MarlinTrk {
  class IMarlinTrkSystem;
}

/**
 * Track finding based on Multi Dimensional Hough Transform (MDHT)
 *
 * @parameter TrackerHitCollections List of the tracker hits collections
 * @parameter OutputTrackCollection Name of the track output collection
 *
 * @author L.Ristori and M. Casarsa
 * @date  9 April 2025
 * @version $Id: $
 */

class HTATrackFinding : public marlin::Processor {
  public:
    virtual Processor* newProcessor() { return new HTATrackFinding; }

    HTATrackFinding(const HTATrackFinding&)            = delete;
    HTATrackFinding& operator=(const HTATrackFinding&) = delete;
    HTATrackFinding();

    /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
    virtual void init();

    /** Called for every run.
   */
    virtual void processRunHeader(LCRunHeader* run);

    /** Called for every event - the working horse.
   */
    virtual void processEvent(LCEvent* evt);

    virtual void check(LCEvent* evt);

    /** Called after data processing for clean up.
   */
    virtual void end();

    //  float fitHitCombination( std::vector<LCCollection*> *inputHits,
    //			   MarlinTrk::IMarlinTrkSystem *trkSystem );

  private:
    float _bField{0.f};

    const dd4hep::rec::SurfaceMap* _surfMap{nullptr};
    std::vector<int>               _trackerBarrelIDs{};
    std::vector<int>               _trackerEndcapIDs{};

    HTArray*    _HTA{nullptr};
    std::string _htaTrainingFile;

    std::mutex _htaMutex;

    std::shared_ptr<UTIL::BitField64> _encoder{};

    MarlinTrk::IMarlinTrkSystem* _trksystem = nullptr;

    bool _MSOn;
    bool _ElossOn;
    bool _SmoothOn;
    bool _extrapolateForward;

    float _particleMass;

    int _minLayersForFit;
    int _minTrackHits;
    float _maxTrackChi2;

    //! Azimuthal angle range for the tracker hits
    std::vector<float> _trackerHitPhiRange{};

    //! Input tracker hit collections
    std::vector<std::string> _inputTrackerHitCollections{};

    //! Output track collection
    std::string _outputTrackCollection{};

};
