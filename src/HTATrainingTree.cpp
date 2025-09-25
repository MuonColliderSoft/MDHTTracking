#include "include/HTATrainingTree.h"

#include <algorithm>

#include "marlin/ProcessorEventSeeder.h"
#include "marlin/Global.h"

#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/TrackerHitPlane.h"
#include "EVENT/SimTrackerHit.h"

#include "UTIL/LCRelationNavigator.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"

HTATrainingTree aHTATrainingTree;

HTATrainingTree::HTATrainingTree() : Processor("HTATrainingTree") {

    _description = "HTATrainingTree produces a flat ntuple for the HT array training.";

    // --- Register the steering parameters

    registerProcessorParameter("MCParticleCollection", "Name of the MCParticle collection", _inputMCParticleCollection,
                               std::string("MCParticle"));

    registerProcessorParameter("TrackerHitCollections", "Name of the tracker hit input collections",
                               _inputTrackerHitCollections, {});

    registerProcessorParameter("TrackerSimHitCollections", "Name of the tracker simhit input collections",
                               _inputTrackerSimHitCollections, {});

    registerProcessorParameter("TrackerHitRelationCollections", "Name of the tracker hit relations input collections",
                               _inputTrackerHitRelCollections, {});

    registerProcessorParameter("ParticleTypes", "Particle types to be saved in the tree", _particleTypes, {});

    registerProcessorParameter("TrackerHitPhiRange", "Azimuthal angle range for the tracker hits", _trackerHitPhiRange,
                               {-M_PI, M_PI});

    registerProcessorParameter("SaveOnlyPartHits", "Save only hits from the primary particles", _saveOnlyPartHits,
                               _saveOnlyPartHits);
}

void HTATrainingTree::init() {

    marlin::Global::EVENTSEEDER->registerProcessor(this);

    // --- Get the detector geometry and the surface map
    dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

    dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension<dd4hep::rec::SurfaceManager>();
    _surfMap                             = surfMan.map("world");

    // --- Create the hit vectors
    _hit_index = new std::vector<unsigned int>();
    _hit_mcp   = new std::vector<int>();
    _hit_id0   = new std::vector<unsigned int>();
    _hit_x     = new std::vector<float>();
    _hit_y     = new std::vector<float>();
    _hit_z     = new std::vector<float>();
    _hit_t     = new std::vector<float>();
    _hit_u     = new std::vector<float>();
    _hit_v     = new std::vector<float>();
    _hit_du    = new std::vector<float>();
    _hit_dv    = new std::vector<float>();
    _hit_ux    = new std::vector<float>();
    _hit_uy    = new std::vector<float>();
    _hit_uz    = new std::vector<float>();
    _hit_vx    = new std::vector<float>();
    _hit_vy    = new std::vector<float>();
    _hit_vz    = new std::vector<float>();

    // --- Define the HTA tree
    _HTAtree = new TTree("HTAtree", "HTA training tree");

    _HTAtree->Branch("nRun", &_nRun, "nRun/I");
    _HTAtree->Branch("nEvt", &_nEvt, "nEvt/I");

    _HTAtree->Branch("part_pdg", &_part_pdg, "part_pdg/I");
    _HTAtree->Branch("part_px", &_part_px, "part_px/F");
    _HTAtree->Branch("part_py", &_part_py, "part_py/F");
    _HTAtree->Branch("part_pz", &_part_pz, "part_pz/F");
    _HTAtree->Branch("part_e", &_part_e, "part_e/F");
    _HTAtree->Branch("part_vx", &_part_vx, "part_vx/F");
    _HTAtree->Branch("part_vy", &_part_vy, "part_vy/F");
    _HTAtree->Branch("part_vz", &_part_vz, "part_vz/F");
    _HTAtree->Branch("part_t0", &_part_t0, "part_t0/F");
    _HTAtree->Branch("part_q", &_part_q, "part_q/F");

    _HTAtree->Branch("n_hit", &_n_hit);
    _HTAtree->Branch("hit_index", _hit_index);
    _HTAtree->Branch("hit_mcp", _hit_mcp);
    _HTAtree->Branch("hit_id0", _hit_id0);
    _HTAtree->Branch("hit_x", _hit_x);
    _HTAtree->Branch("hit_y", _hit_y);
    _HTAtree->Branch("hit_z", _hit_z);
    _HTAtree->Branch("hit_t", _hit_t);
    _HTAtree->Branch("hit_u", _hit_u);
    _HTAtree->Branch("hit_v", _hit_v);
    _HTAtree->Branch("hit_du", _hit_du);
    _HTAtree->Branch("hit_dv", _hit_dv);
    _HTAtree->Branch("hit_ux", _hit_ux);
    _HTAtree->Branch("hit_uy", _hit_uy);
    _HTAtree->Branch("hit_uz", _hit_uz);
    _HTAtree->Branch("hit_vx", _hit_vx);
    _HTAtree->Branch("hit_vy", _hit_vy);
    _HTAtree->Branch("hit_vz", _hit_vz);

    // --- Print the initial parameters
    printParameters();
}

void HTATrainingTree::processRunHeader(LCRunHeader* /*run*/) {}

void HTATrainingTree::processEvent(LCEvent* evt) {

    // --- By default, match tracker reco and sim hits
    bool match_SimRecoHits = true;

    // --- Check if the simhit and relation collections are provided and if their numbers are consistent
    if (_inputTrackerSimHitCollections.size() == 0 || _inputTrackerHitRelCollections.size() == 0) {
        match_SimRecoHits = false;
        streamlog_out(WARNING) << " TrackerSimHitCollections or TrackerHitRelationCollections not available:"
                               << " no sim-reco matching!" << std::endl;
    } else {
        if (_inputTrackerSimHitCollections.size() != _inputTrackerHitCollections.size() ||
            _inputTrackerHitRelCollections.size() != _inputTrackerHitCollections.size()) {
            std::stringstream err_msg;
            err_msg << "Mismatch between the recohits, simhits, and relations input collections" << std::endl;
            throw EVENT::Exception(err_msg.str());
        }
    }

    // --- Get the run and event numbers
    _nRun = evt->getRunNumber();
    _nEvt = evt->getEventNumber();

    // --- Get the MC particles collection
    LCCollection* inputMCParticles = nullptr;

    try {
        inputMCParticles = evt->getCollection(_inputMCParticleCollection);
    } catch (lcio::DataNotAvailableException& e) {
        streamlog_out(WARNING) << _inputMCParticleCollection << " collection not available" << std::endl;
        return;
    }

    // --- Get the tracker hit collections
    const unsigned int         nTrackerHitCol = _inputTrackerHitCollections.size();
    std::vector<LCCollection*> inputHitColls(nTrackerHitCol);
    std::vector<LCCollection*> inputSimHitColls(nTrackerHitCol);
    std::vector<LCCollection*> inputHitRelColls(nTrackerHitCol);

    for (unsigned int icol = 0; icol < nTrackerHitCol; ++icol) {

        try {
            inputHitColls[icol] = evt->getCollection(_inputTrackerHitCollections[icol]);
        } catch (lcio::DataNotAvailableException& e) {
            streamlog_out(WARNING) << _inputTrackerHitCollections[icol] << " collection not available" << std::endl;
            return;
        }

        if (match_SimRecoHits) {

            try {
                inputSimHitColls[icol] = evt->getCollection(_inputTrackerSimHitCollections[icol]);
            } catch (lcio::DataNotAvailableException& e) {
                streamlog_out(WARNING) << _inputTrackerSimHitCollections[icol] << " collection not available"
                                       << std::endl;
                return;
            }

            try {
                inputHitRelColls[icol] = evt->getCollection(_inputTrackerHitRelCollections[icol]);
            } catch (lcio::DataNotAvailableException& e) {
                streamlog_out(WARNING) << _inputTrackerHitRelCollections[icol] << " collection not available"
                                       << std::endl;
                return;
            }

        } // if match_SimRecoHits

    } // icol loop

    // --- Loop over the MC particles
    for (int ipart = 0; ipart < inputMCParticles->getNumberOfElements(); ++ipart) {

        MCParticle* part = dynamic_cast<MCParticle*>(inputMCParticles->getElementAt(ipart));

        // Keep only the generator-level particles:
        if (part->getGeneratorStatus() != 1)
            continue;

        _part_pdg = part->getPDG();

        // Check the particle type
        if (std::find(_particleTypes.begin(), _particleTypes.end(), fabs(_part_pdg)) == _particleTypes.end())
            continue;

        _part_px = part->getMomentum()[0];
        _part_py = part->getMomentum()[1];
        _part_pz = part->getMomentum()[2];
        _part_e  = part->getEnergy();
        _part_vx = part->getVertex()[0];
        _part_vy = part->getVertex()[1];
        _part_vz = part->getVertex()[2];
        _part_t0 = part->getTime();
        _part_q  = part->getCharge();

    } // ipart loop

    // --- Loop over the tracker hits

    _n_hit = 0;
    for (unsigned int icol = 0; icol < inputHitColls.size(); ++icol) {

        LCCollection* hit_col = inputHitColls[icol];
        if (!hit_col)
            continue;

        UTIL::LCRelationNavigator* hit_rel = nullptr;
        if (match_SimRecoHits)
            hit_rel = new UTIL::LCRelationNavigator(inputHitRelColls[icol]);

        for (int ihit = 0; ihit < hit_col->getNumberOfElements(); ++ihit) {

            TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(hit_col->getElementAt(ihit));

            // --- Keep only hits in the given phi range
            float hit_phi = std::atan2(hit->getPosition()[1], hit->getPosition()[0]);
            if (hit_phi < _trackerHitPhiRange[0] || hit_phi > _trackerHitPhiRange[1])
                continue;

            dd4hep::rec::SurfaceMap::const_iterator si   = _surfMap->find(hit->getCellID0());
            dd4hep::rec::ISurface*                  surf = (si != _surfMap->end() ? si->second : nullptr);
            if (!surf)
                continue;

            dd4hep::rec::Vector3D globalPoint(hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2]);
            dd4hep::rec::Vector2D localPoint = surf->globalToLocal(dd4hep::mm * globalPoint);

            dd4hep::rec::Vector3D u = surf->u();
            dd4hep::rec::Vector3D v = surf->v();

            int hit_mother = 0;
            if (match_SimRecoHits && hit_rel != nullptr) {

                const LCObjectVec& simHitVector = hit_rel->getRelatedToObjects(hit);

                if (simHitVector.size() != 0) {
                    SimTrackerHit* simhit = dynamic_cast<SimTrackerHit*>(simHitVector.at(0));
                    hit_mother            = simhit->getMCParticle()->getPDG();
                    int part_status       = simhit->getMCParticle()->getGeneratorStatus();

                    // If _saveOnlyPartHits  skip hits from secondary particles
                    if (_saveOnlyPartHits &&
                        (part_status != 1 || std::find(_particleTypes.begin(), _particleTypes.end(),
                                                       fabs(hit_mother)) == _particleTypes.end()))
                        continue;

                } // if simHitVector.size()

            } // if match_SimRecoHits && hit_rel != nullptr

            _hit_index->push_back(ihit);
            _hit_mcp->push_back(hit_mother);
            _hit_id0->push_back(hit->getCellID0());
            _hit_x->push_back(globalPoint[0]);
            _hit_y->push_back(globalPoint[1]);
            _hit_z->push_back(globalPoint[2]);
            _hit_t->push_back(hit->getTime());
            _hit_u->push_back(localPoint[0] / dd4hep::mm);
            _hit_v->push_back(localPoint[1] / dd4hep::mm);
            _hit_du->push_back(hit->getdU());
            _hit_dv->push_back(hit->getdV());
            _hit_ux->push_back(u[0]);
            _hit_uy->push_back(u[1]);
            _hit_uz->push_back(u[2]);
            _hit_vx->push_back(v[0]);
            _hit_vy->push_back(v[1]);
            _hit_vz->push_back(v[2]);

            _n_hit++;

        } // ihit loop

        if (match_SimRecoHits)
            delete hit_rel;

    } // icol loop

    // --- Fill the tree
    _HTAtree->Fill();

    // --- Clear the hit vectors

    _hit_index->clear();
    _hit_mcp->clear();
    _hit_id0->clear();
    _hit_x->clear();
    _hit_y->clear();
    _hit_z->clear();
    _hit_t->clear();
    _hit_u->clear();
    _hit_v->clear();
    _hit_du->clear();
    _hit_dv->clear();
    _hit_ux->clear();
    _hit_uy->clear();
    _hit_uz->clear();
    _hit_vx->clear();
    _hit_vy->clear();
    _hit_vz->clear();
}

void HTATrainingTree::check(LCEvent* /*evt*/) {}

void HTATrainingTree::end() {

    // --- Clean up the heap

    delete _hit_index;
    delete _hit_mcp;
    delete _hit_id0;
    delete _hit_x;
    delete _hit_y;
    delete _hit_z;
    delete _hit_t;
    delete _hit_u;
    delete _hit_v;
    delete _hit_du;
    delete _hit_dv;
    delete _hit_ux;
    delete _hit_uy;
    delete _hit_uz;
    delete _hit_vx;
    delete _hit_vy;
    delete _hit_vz;
}
