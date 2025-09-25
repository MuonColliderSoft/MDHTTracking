#ifndef MDHT_PARAMETERS_H
#define MDHT_PARAMETERS_H

//
//  Parameters.h
//  HTA_New_Training
//
//  This version created by Luciano Ristori on 3/21/25
//

#include <string>

class Parameters {

  public:
    ////////////////////////////////////////////////////////////////////////
    // TRAINING SECTION ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    //std::string train_inputTrackFileName = "ntu_muongun_pt3_theta10-170_phi0-30_20M.root";
    std::string train_inputTrackFileName = "ntu_muongun_pt3_theta10-170_phi0-30_training_30M.root";
    //std::string train_inputTrackFileName = "ntu_muongun_pt3_theta10-170_phi0-30_training_1M.root";

    long unsigned train_maxTracks   = 0;     // 0 == no limit
    bool          train_printTracks = false; // print all tracks while reading them

    // Hit coordinates distributions are diagonalized to find principal components
    // before finding minimum and maximum

    bool train_Diagonalize = false;

    // Write files with tracks in each cell at TrainingPhase == 3

    bool train_WriteFiles = false;

    // At the end of training, a summary of results is printed

    bool train_Summary = false;

    // Minimum number of tracks to be used to train each cell of the array

    const static unsigned train_nTracksPerCell = 100;

    // Histogram files for three training phases

    std::string train_histFileName1 = "AAATrainingHists1.root";
    std::string train_histFileName2 = "AAATrainingHists2.root";
    std::string train_histFileName3 = "AAATrainingHists3.root";

    // File containing all the data to init the Hough Transform array
    // This is where all the data from training are written to

    std::string train_dataFileName = "HTAdata.txt";

    // track parameters for Hough Transform Array

    const static unsigned HTA_NphiBins   = 15;
    double                HTA_t_phi      = 0.262; // track phi mean
    double                HTA_t_deltaPhi = 0.262; // track delta phi

    const static unsigned HTA_NetaBins   = 360;
    double                HTA_t_eta      = 0.;  // track eta mean
    double                HTA_t_deltaEta = 2.5; // track delta eta

    const static unsigned HTA_NinvptBins  = 6;
    double                HTA_t_invPt_max = 1. / 3.0;  // track invPt max Gev/c^(-1)
    double                HTA_t_invPt_min = -1. / 3.0; // track invPt min Gev/c^(-1)

    // coordinates of a single HTA cell and detector layer
    // whose hit coordinates we want to plot in 3-D

    unsigned HTA_plotBinX = 14;
    unsigned HTA_plotBinY = 189;
    unsigned HTA_plotBinZ = 4;
    unsigned HTA_plotLay  = 302;

    ////////////////////////////////////////////////////////////////////////
    //   RECONSTRUCTION SECTION  ///////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    // file with tracks to build events to be reconstructed

    //std::string reco_inputTrackFileName = "ntu_muongun_pt3_theta10-170_phi0-30_100k.root";

    //std::string reco_inputTrackFileName = "ntu_muongun_pt3_theta10-170_phi0-30_dz1p5-10k.root";

    std::string reco_inputTrackFileName = "ntu_muongun_pt3_theta10-170_phi0-30_dz1p5-10k-2.root";

    // Plot reconstructed tracks in 3D

    bool reco_plotTrack3D = false;

    // Number of events to be generated for simulation

    unsigned reco_nEvents = 100;

    // number of tracks to be generated for each event

    unsigned reco_nTracks = 1; // Number of tracks per event

    // Background events to be generated

    double reco_backGnd = -1.; // -1. == no BIB hit

    // perform track fitting of candidates

    bool     reco_TrackFit         = true;
    unsigned reco_minLayersForFit  = 5;
    unsigned reco_maxDroppedLayers = 1;
    double   reco_chi2Cut          = 1000.;
    bool     reco_massFit          = false;
    double   reco_massFitMaxP      = 3.; // GeV/c

    // print candidates for every event

    bool reco_printCandidates = true;

    // optimization mode for HTA fill = 0 (safe and slow) or 1 (faster) or 2 (fastest)

    int reco_fillMode = 0;

    // create a file with data to 3D plot track candidates
    // this file is the input for PlotTracks

    bool        reco_PlotTracks       = false;
    std::string reco_plotDataFileName = "PlotData.txt";

    // random generator seeds

    long int reco_randomSeed = 121348;

    // data file from where HT array is initialized from

    std::string reco_dataFileName = "HTAdata_ntu_muongun_pt3_theta10-170_phi0-30_training_30M.txt";

    // data file with all Bib hits

    std::string reco_bibFileName = "BIBdata.txt";

    // Reconstruction histogram file

    std::string reco_histFileName = "AAAReconstruction.root";

    bool reco_fillHitHistograms = false; // selective histogram filling

    // printing control

    bool reco_verbose = false;

    ////////////////////////////////////////////////////////////////////////
    // RANGE FOR TRACK PARAMETER HISTOGRAMS ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    double geo_gen_t_phi       = 0.262;     // track phi mean
    double geo_gen_t_deltaPhi  = 0.262;     // track delta phi
    double geo_gen_t_eta       = 0.;        // track eta mean
    double geo_gen_t_deltaEta  = 2.0;       // track delta eta
    double geo_gen_t_invPt_max = 1. / 1.5;  // track invPt max Gev/c^(-1)
    double geo_gen_t_invPt_min = -1. / 1.5; // track invPt min Gev/c^(-1)
    double geo_gen_t_x0        = 0.0;       // track mean x0
    double geo_gen_t_y0        = 0.0;       // track mean y0
    double geo_gen_t_z0        = 0.0;       // track mean z0
    double geo_gen_t_t0        = 0.0;       // track mean t0
    double geo_gen_t_deltaX0   = 0.0;       // track sigma x0
    double geo_gen_t_deltaY0   = 0.0;       // track sigma y0
    double geo_gen_t_deltaZ0   = 1.5;       // track sigma z0
    double geo_gen_t_deltaT0   = 1.5;       // track sigma t0 mm

    ////////////////////////////////////////////////////////////////////////
    // MISCELLANEA /////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    // ignore vertex

    bool ignoreVertex = false;

    // magnetic field in Tesla

    const double magneticField = 5.0; // Tesla

    ////////////////////////////////////////////////////////////////////////
    // CONSTRUCTOR//////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    Parameters() {}
};

#endif // MDHT_PARAMETERS_H
