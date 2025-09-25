#include "include/HTATrackFinding.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"

#include "EVENT/TrackerHitPlane.h"

#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/TrackStateImpl.h"
#include "IMPL/TrackImpl.h"

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/MarlinTrkUtils.h"

#include "UTIL/LCTrackerConf.h"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>


Parameters par;

HTATrackFinding aHTATrackFinding;

HTATrackFinding::HTATrackFinding() : Processor("HTATrackFinding") {

  _description = "Track finding algorithm based on a multi-dimensional Hough Transform";

  // --- Register the steering parameters

  registerProcessorParameter("TrackerHitCollections",
			     "Name of the tracker hit input collections",
			     _inputTrackerHitCollections, {});

  registerOutputCollection(LCIO::TRACK, "OutputTrackCollection",
			   "Name of the output track collection",
			   _outputTrackCollection, std::string("MDHT_tracks"));

  registerProcessorParameter("HTATrainingFile",
			     "Name of HTA training file",
			     _htaTrainingFile, std::string("HTAdata.txt"));

  registerProcessorParameter("TrackerHitPhiRange",
			     "Azimuthal angle range for the tracker hits",
			     _trackerHitPhiRange, {-M_PI, M_PI});

  registerProcessorParameter("MinLayersForFit",
			     "Minimum number of hits on different layers to fit the track",
			     _minLayersForFit, 5);
    
  registerProcessorParameter("MaxTrackChi2",
			     "Maximum chi2 value to keep a fitted track",
			     _maxTrackChi2, 45.f);
  
  registerProcessorParameter("MinTrackHits",
			     "Minimum number of hits to keep a fitted track",
			     _minTrackHits, 5);
  
  registerProcessorParameter("MultipleScatteringOn",
			     "Use multiple scattering in the track fit",
			     _MSOn, true);

  registerProcessorParameter("EnergyLossOn",
			     "Use energy loss in the track fit",
			     _ElossOn, true);

  registerProcessorParameter("SmoothOn",
			     "Smooth all mesurement sites in the track fit",
			     _SmoothOn, false);

  registerProcessorParameter("extrapolateForward",
			     "If true, extrapolate in the forward direction "
			     "(in-out), otherwise backward (out-in)",
			     _extrapolateForward, true);

  registerProcessorParameter("ParticleMass",
			     "Particle mass used in the track fit (default is the pion mass)",
			     _particleMass, 0.13957018f);
    
}

void HTATrackFinding::init() {

  // --- Print the initial parameters
  printParameters();

  // --- Create and initialize the HT array
  _HTA = new HTArray();
  _HTA->initHists();

  std::cout << " ======================================================================================= "
	    << std::endl;
  _HTA->print(std::cout);
  std::cout << " ======================================================================================= "
	    << std::endl;

  streamlog_out(MESSAGE) << "Opening the HTA training file: " << _htaTrainingFile << " ..." << std::endl;

  std::ifstream hta_file(_htaTrainingFile);
  if (!hta_file) {
    std::stringstream err_msg;
    err_msg << "HTA training file " << _htaTrainingFile << " not found!" << std::endl;
    throw std::runtime_error(err_msg.str());
  }

  int retcode = _HTA->read(hta_file);
  if (retcode < 0) {
    std::stringstream err_msg;
    err_msg << "The HTA training file does not match the current HTA configuration" << std::endl;
    throw EVENT::Exception(err_msg.str());
  }

  // --- Get the detector geometry and the surface map
  const dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

  dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension<dd4hep::rec::SurfaceManager>();
  _surfMap = surfMan.map("world");

  // --- Get the magnetic field value
  const double pos[3]       = {0., 0., 0.};
  double       bFieldVec[3] = {0., 0., 0.};
  theDetector.field().magneticField(pos, bFieldVec);
  _bField = bFieldVec[2] / dd4hep::tesla;

  // --- Get the tracker subdetector IDs
  _trackerBarrelIDs.push_back(theDetector.constant<unsigned int>("DetID_VXD_Barrel"));
  _trackerBarrelIDs.push_back(theDetector.constant<unsigned int>("DetID_IT_Barrel"));
  _trackerBarrelIDs.push_back(theDetector.constant<unsigned int>("DetID_OT_Barrel"));
  _trackerEndcapIDs.push_back(theDetector.constant<unsigned int>("DetID_VXD_Endcap"));
  _trackerEndcapIDs.push_back(theDetector.constant<unsigned int>("DetID_IT_Endcap"));
  _trackerEndcapIDs.push_back(theDetector.constant<unsigned int>("DetID_OT_Endcap"));

  // --- Define a CellID encoder for the tracker hits
  _encoder = std::make_shared<UTIL::BitField64>(lcio::LCTrackerCellID::encoding_string());

  // --- Create the Marlin tracking system
  _trksystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");

  _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useQMS, _MSOn);
  _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, _ElossOn);
  _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, _SmoothOn);
  _trksystem->init();
}

void HTATrackFinding::processRunHeader(LCRunHeader* /*run*/) {}

void HTATrackFinding::processEvent(LCEvent* evt) {

  // --- Reset the HT array
  _HTA->reset();

  // --- Set the tracking system configuration for this event
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS>       mson(_trksystem, _MSOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx>      elosson(_trksystem, _ElossOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(_trksystem, _SmoothOn);
    
  // --- Get the tracker hit collections
  const unsigned int         nTrackerHitCol = _inputTrackerHitCollections.size();
  std::vector<LCCollection*> inputHitColls(nTrackerHitCol);

  for (unsigned int icol = 0; icol < nTrackerHitCol; ++icol) {

    try {
      inputHitColls[icol] = evt->getCollection(_inputTrackerHitCollections[icol]);
    } catch (lcio::DataNotAvailableException& e) {
      streamlog_out(WARNING) << _inputTrackerHitCollections[icol] << " collection not available" << std::endl;
      return;
    }

  } // icol loop

 
    // --- Make the output track collection
  LCCollectionVec* trackCollection = new LCCollectionVec(LCIO::TRACK);
  // to store tracks instead of pointers in the new collection
  trackCollection->setSubset(false);
  // to point back to the hits
  LCFlagImpl trkFlag(0);
  trkFlag.setBit(LCIO::TRBIT_HITS);
  trackCollection->setFlag(trkFlag.getFlag());

  // --- Loop over the tracker hits and fill the HT array
  for (unsigned int icol=0; icol<inputHitColls.size(); ++icol){
    
    LCCollection* hit_col = inputHitColls[icol];
    if( !hit_col ) continue ;

    tbb::parallel_for(tbb::blocked_range<int>(0, hit_col->getNumberOfElements()),
		      [&](const tbb::blocked_range<int>& range) {
			for (int ihit = range.begin(); ihit < range.end(); ++ihit) {

			  const auto* hit = dynamic_cast<TrackerHitPlane*>(hit_col->getElementAt(ihit));
			  if (!hit) continue;

			  // --- Keep only hits in the given phi range
			  float hit_phi = std::atan2(hit->getPosition()[1], hit->getPosition()[0]);
			  if (hit_phi < _trackerHitPhiRange[0] || hit_phi > _trackerHitPhiRange[1]) continue;

			  const int cellID = hit->getCellID0();

			  UTIL::BitField64 localEncoder(lcio::LCTrackerCellID::encoding_string());
			  localEncoder.setValue(cellID);
			  const int system = localEncoder[lcio::LCTrackerCellID::subdet()];
			  const int side   = (localEncoder[lcio::LCTrackerCellID::side()] & 0x3);
			  const int layer  = localEncoder[lcio::LCTrackerCellID::layer()];

			  auto si = _surfMap->find(cellID);
			  if (si == _surfMap->end()) continue;
			  dd4hep::rec::ISurface* surf = si->second;

			  dd4hep::rec::Vector3D globalPoint(hit->getPosition()[0], hit->getPosition()[1],
							    hit->getPosition()[2]);
			  dd4hep::rec::Vector2D localPoint = surf->globalToLocal(dd4hep::mm * globalPoint);

			  dd4hep::rec::Vector3D u = surf->u();
			  dd4hep::rec::Vector3D v = surf->v();

			  char hitType = (std::find(_trackerBarrelIDs.begin(), _trackerBarrelIDs.end(), system) !=
					  _trackerBarrelIDs.end())
			    ? 'B'
			    : (std::find(_trackerEndcapIDs.begin(), _trackerEndcapIDs.end(), system) !=
			       _trackerEndcapIDs.end())
			    ? 'D'
			    : 'G';

			  int layerInd = 100 * system + 10 * side + layer;
			  int trackInd = 0;

			  Hit hta_hit(hitType, cellID, icol, ihit, layerInd, trackInd,
				      localPoint[0] / dd4hep::mm, localPoint[1] / dd4hep::mm,
				      hit->getTime() * SpeedOfLight, globalPoint[0], globalPoint[1],
				      globalPoint[2], u[0], u[1], u[2], v[0], v[1], v[2],
				      hit->getdU(), hit->getdV());

			  {
			    std::lock_guard<std::mutex> lock(_htaMutex);
			    _HTA->fill(hta_hit);
			  }
			}
		      }); // ihit parallelized loop
  } // icol loop


  // --- Get the best HTA cell

  unsigned iphi {};
  unsigned ieta {};
  unsigned ipt {};

  unsigned n_hits =_HTA->getBestCell(iphi, ieta, ipt);

  if ( n_hits > 0 )
    _HTA->ArrElem[iphi][ieta][ipt].printHits(std::cout);

  if ( n_hits < _minLayersForFit ) return;
    

  // --- Fit the cell hit candidates

  double chi2, phi, eta, invPt, z0, t0, beta;
  std::vector <long int> goodFitHitList;
  unsigned nHitsFit;
						
  int retCodeFit = _HTA->ArrElem[iphi][ieta][ipt].
    fitCandidate(&(_HTA->allHits),chi2, phi, eta, invPt, z0, t0, beta, nHitsFit, goodFitHitList);

  if(retCodeFit != 0){
    std::cout << "FIT FAILED with retCodeFit = " << retCodeFit << std::endl;
    return;
  }

  // --- Refit the HTA track

  auto marlin_trk = std::unique_ptr<MarlinTrk::IMarlinTrack>(_trksystem->createTrack());
  marlin_trk->setMass(_particleMass);
    
  TrackerHitVec trk_hits;
  for (size_t ihit = 0; ihit < goodFitHitList.size(); ++ihit) {

    Hit hta_hit = _HTA->allHits[goodFitHitList[ihit]];

    TrackerHit* hit = dynamic_cast<TrackerHit*>(inputHitColls[hta_hit.collInd]->getElementAt(hta_hit.hitInd));

    int status = marlin_trk->addHit(hit);
    if (status != 0) continue;

    trk_hits.push_back(hit);

  } // ihit loop
    
    
    // Initialize the MarlinTrkSystem fitter
    
  float initialTrackError_d0    = 10000.;
  float initialTrackError_phi0  = 4./(_HTA->phiStep*_HTA->phiStep);
  float initialTrackError_omega = 4./(_HTA->invPtStep*_HTA->invPtStep);
  float initialTrackError_z0    = 100;
  float initialTrackError_tanl  = 4./(cosh(_HTA->getCell_Eta(ieta))*cosh(_HTA->getCell_Eta(ieta)) *
				      _HTA->etaStep*_HTA->etaStep);
    
  const std::vector<float> cov_matrix = { initialTrackError_d0,
    0., initialTrackError_phi0,
    0., 0., initialTrackError_omega,
    0., 0., 0., initialTrackError_z0,
    0., 0., 0., 0., initialTrackError_tanl };
    
  const bool fit_direction = ( _extrapolateForward                ?
			       MarlinTrk::IMarlinTrack::forward   :
			       MarlinTrk::IMarlinTrack::backward );
    
  auto final_trk = std::unique_ptr<IMPL::TrackImpl>(new IMPL::TrackImpl());
    
  int fit_status = MarlinTrk::createFinalisedLCIOTrack(marlin_trk.get(), trk_hits, final_trk.get(),
						       fit_direction, cov_matrix, _bField, 100.);

  if ( fit_status==0 && final_trk->getChi2()<_maxTrackChi2 ){

    // --- Add the track to the output track collection
    trackCollection->addElement(final_trk.release());

  }

  streamlog_out(MESSAGE4) << " Final number of tracks " << trackCollection->getNumberOfElements()
			  << std::endl;
  evt->addCollection(trackCollection, _outputTrackCollection);
}

void HTATrackFinding::check(LCEvent* /*evt*/) {}

void HTATrackFinding::end() { delete _HTA; }
