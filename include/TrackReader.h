
#ifndef MDHT_TRACKREADER_H
#define MDHT_TRACKREADER_H

#include <iostream>

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"

#include "include/Track.h"
#include "include/Hit.h"
#include "include/CellIDtoLayer.h"
#include "include/GlobalConstants.h"

/*

 int n_hit;
  std::vector<unsigned int>* hit_index = nullptr;
  std::vector<int>* hit_mcp = nullptr;
  std::vector<unsigned int>* hit_id0 = nullptr;
  std::vector<float>* hit_x = nullptr;
  std::vector<float>* hit_y = nullptr;
  std::vector<float>* hit_z = nullptr;
  std::vector<float>* hit_t = nullptr;
  std::vector<float>* hit_u = nullptr;
  std::vector<float>* hit_v = nullptr;
  std::vector<float>* hit_du = nullptr;
  std::vector<float>* hit_dv = nullptr;
  std::vector<float>* hit_ux = nullptr;
  std::vector<float>* hit_uy = nullptr;
  std::vector<float>* hit_uz = nullptr;
  std::vector<float>* hit_vx = nullptr;
  std::vector<float>* hit_vy = nullptr;
  std::vector<float>* hit_vz = nullptr;

  tree->SetBranchAddress("n_hit", &n_hit);
  tree->SetBranchAddress("hit_index", &hit_index);
  tree->SetBranchAddress("hit_mcp", &hit_mcp);
  tree->SetBranchAddress("hit_id0", &hit_id0);
  tree->SetBranchAddress("hit_x", &hit_x);
  tree->SetBranchAddress("hit_y", &hit_y);
  tree->SetBranchAddress("hit_z", &hit_z);
  tree->SetBranchAddress("hit_t", &hit_t);
  tree->SetBranchAddress("hit_u", &hit_u);
  tree->SetBranchAddress("hit_v", &hit_v);
  tree->SetBranchAddress("hit_du", &hit_du);
  tree->SetBranchAddress("hit_dv", &hit_dv);
  tree->SetBranchAddress("hit_ux", &hit_ux);
  tree->SetBranchAddress("hit_uy", &hit_uy);
  tree->SetBranchAddress("hit_uz", &hit_uz);
  tree->SetBranchAddress("hit_vx", &hit_vx);
  tree->SetBranchAddress("hit_vy", &hit_vy);
  tree->SetBranchAddress("hit_vz", &hit_vz);
  
 */

class TrackReader {

  public:
    TFile*   input_file;
    TTree*   tree;
    Long64_t nEntries;
    Long64_t ientry;

    unsigned int min_system = 1e6;
    unsigned int min_side   = 1e6;
    unsigned int min_layer  = 1e6;
    unsigned int min_module = 1e6;
    unsigned int min_sensor = 1e6;

    unsigned int max_system = 0;
    unsigned int max_side   = 0;
    unsigned int max_layer  = 0;
    unsigned int max_module = 0;
    unsigned int max_sensor = 0;

    // MC particle

    int   part_pdg;
    float part_px;
    float part_py;
    float part_pz;
    float part_e;
    float part_vx;
    float part_vy;
    float part_vz;
    float part_t0;
    float part_q;

    // tracker hits

    int n_hit;

    std::vector<unsigned int>* hit_index = nullptr;
    std::vector<int>*          hit_mcp   = nullptr;
    std::vector<unsigned int>* hit_id0   = nullptr;
    std::vector<float>*        hit_x     = nullptr;
    std::vector<float>*        hit_y     = nullptr;
    std::vector<float>*        hit_z     = nullptr;
    std::vector<float>*        hit_t     = nullptr;
    std::vector<float>*        hit_u     = nullptr;
    std::vector<float>*        hit_v     = nullptr;
    std::vector<float>*        hit_du    = nullptr;
    std::vector<float>*        hit_dv    = nullptr;
    std::vector<float>*        hit_ux    = nullptr;
    std::vector<float>*        hit_uy    = nullptr;
    std::vector<float>*        hit_uz    = nullptr;
    std::vector<float>*        hit_vx    = nullptr;
    std::vector<float>*        hit_vy    = nullptr;
    std::vector<float>*        hit_vz    = nullptr;

    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    TrackReader(std::string fileName) { // constructor

        //std::cout << "CONSTRUCTOR" << std::endl;

        // --- open the root file (assign to member input_file)
        input_file = TFile::Open(fileName.c_str());
        if (!input_file) {
            throw std::runtime_error("File not found!");
            return;
        }

        std::cout << "FILE " << fileName << " opened successfully" << std::endl;

        // --- retrieve the tree
        TTree* xtree = (TTree*)input_file->Get("HTAtree");
        tree         = xtree;

        std::cout << "tree pointer: " << xtree << std::endl;

        // --- define branch variables and set the branch addresses

        tree->SetBranchAddress("part_pdg", &part_pdg);
        tree->SetBranchAddress("part_px", &part_px);
        tree->SetBranchAddress("part_py", &part_py);
        tree->SetBranchAddress("part_pz", &part_pz);
        tree->SetBranchAddress("part_e", &part_e);
        tree->SetBranchAddress("part_vx", &part_vx);
        tree->SetBranchAddress("part_vy", &part_vy);
        tree->SetBranchAddress("part_vz", &part_vz);
        tree->SetBranchAddress("part_t0", &part_t0);
        tree->SetBranchAddress("part_q", &part_q);

        tree->SetBranchAddress("n_hit", &n_hit);

        tree->SetBranchAddress("hit_index", &hit_index);
        tree->SetBranchAddress("hit_mcp", &hit_mcp);
        tree->SetBranchAddress("hit_id0", &hit_id0);
        tree->SetBranchAddress("hit_x", &hit_x);
        tree->SetBranchAddress("hit_y", &hit_y);
        tree->SetBranchAddress("hit_z", &hit_z);
        tree->SetBranchAddress("hit_t", &hit_t);
        tree->SetBranchAddress("hit_u", &hit_u);
        tree->SetBranchAddress("hit_v", &hit_v);
        tree->SetBranchAddress("hit_du", &hit_du);
        tree->SetBranchAddress("hit_dv", &hit_dv);
        tree->SetBranchAddress("hit_ux", &hit_ux);
        tree->SetBranchAddress("hit_uy", &hit_uy);
        tree->SetBranchAddress("hit_uz", &hit_uz);
        tree->SetBranchAddress("hit_vx", &hit_vx);
        tree->SetBranchAddress("hit_vy", &hit_vy);
        tree->SetBranchAddress("hit_vz", &hit_vz);

        nEntries = tree->GetEntries();
        ientry   = 0;

    } // end constructor

    bool read(Track& thisTrack) {

        if (ientry != nEntries) {

            tree->GetEntry(ientry++);

            // loop over the MC particles
            const float part_pt    = sqrt(part_px * part_px + part_py * part_py);
            const float part_phi   = atan2(part_py, part_px);
            const float part_theta = atan2(part_pt, part_pz);
            const float part_eta   = asinh(1. / tan(part_theta));

            thisTrack.init(0., part_vx, part_vy, part_vz, part_t0 * SpeedOfLight, part_q / part_pt, part_eta,
                           part_phi); // specific track constructor

            // loop over all hits

            for (int ihit = 0; ihit != n_hit; ++ihit) {

                // CellID encoding: "system:5,side:-2,layer:6,module:11,sensor:8"
                //  system:
                //    VXD barrel = 1
                //    VXD endcap = 2
                //    IT barrel  = 3
                //    IT endcap  = 4
                //    OT barrel  = 5
                //    OT endcap  = 6

                int CellID = (*hit_id0)[ihit];
                // decode fields but only use layerInd below; silence unused-variable warnings
                const unsigned int system = (unsigned)(CellID & 0x1f);
                (void)system;
                const unsigned int side = (unsigned)((CellID >> 5) & 0x3);
                (void)side;
                const unsigned int layer = (unsigned)((CellID >> 7) & 0x3f);
                (void)layer;
                const unsigned int module = (unsigned)((CellID >> 13) & 0x7ff);
                (void)module;
                const unsigned int sensor = (unsigned)((CellID >> 24) & 0xff);
                (void)sensor;

                int layerInd = CellIDtoLayer(CellID);

                int trackInd;
                (void)trackInd;
                double x1, x2, t, x, y, z;
                double u1x, u1y, u1z, u2x, u2y, u2z;
                double errx1, errx2;

                x = (*hit_x)[ihit];
                y = (*hit_y)[ihit];
                z = (*hit_z)[ihit];
                t = (*hit_t)[ihit] * SpeedOfLight;

                //if(	(*hit_du)[ihit] < (*hit_dv)[ihit]){
                if (true) {

                    x1 = (*hit_u)[ihit];
                    x2 = (*hit_v)[ihit];

                    u1x = (*hit_ux)[ihit];
                    u1y = (*hit_uy)[ihit];
                    u1z = (*hit_uz)[ihit];

                    u2x = (*hit_vx)[ihit];
                    u2y = (*hit_vy)[ihit];
                    u2z = (*hit_vz)[ihit];

                    errx1 = (*hit_du)[ihit];
                    errx2 = (*hit_dv)[ihit];
                } else {

                    x2 = (*hit_u)[ihit];
                    x1 = (*hit_v)[ihit];

                    u2x = (*hit_ux)[ihit];
                    u2y = (*hit_uy)[ihit];
                    u2z = (*hit_uz)[ihit];

                    u1x = (*hit_vx)[ihit];
                    u1y = (*hit_vy)[ihit];
                    u1z = (*hit_vz)[ihit];

                    errx2 = (*hit_du)[ihit];
                    errx1 = (*hit_dv)[ihit];
                }

                Hit thisHit('H', CellID, layerInd, 0, x1, x2, t, x, y, z, u1x, u1y, u1z, u2x, u2y, u2z, errx1, errx2);

                thisTrack.hitList.push_back(thisHit);

                //skip sorting for now

                //std::sort(thisTrack.hitList.begin(), thisTrack.hitList.end(), [](const Hit& a, const Hit& b) {
                //return ((a.CellID >> 7) & 0x3f) < ((b.CellID >> 7) & 0x3f);
                //});

            } // end loop on ihit

            return false;

        } // ientry good

        else {

            // EOF close the root file

            //input_file->Close();

            return true; //EOF
        }

    } // end read

}; //TrainingTracksReader

#endif // MDHT_TRACKREADER_H
