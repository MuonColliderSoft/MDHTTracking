#include <iostream>

#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

static constexpr Double_t MagneticField = 5.;
static constexpr Double_t K = 0.000299792458*MagneticField;

TH1F h_muon_pt("muon_pt","muon p_{T};p_{T} [GeV]", 100., 0., 100);
TH1F h_muon_invpt("muon_invpt","muon p_{T}^{-1};p_{T}^{-1} [GeV^{-1}]", 100., -0.4, 0.4);
TH1F h_muon_phi("muon_phi","muon #phi;#phi [rad]", 100., -3.1415, 3.1415);
TH1F h_muon_theta("muon_theta","muon #theta;#theta [rad]",100., 0., 3.1415);
TH1F h_muon_eta("muon_eta","muon #eta;#eta",100., -2.5, 2.5);
TH1F h_muon_vz("muon_vz","muon vz;vz [mm]", 100, -10., 10.);
TH1F h_muon_t0("muon_t0","muon T_{0}; T_{0} [ns]", 100, -0.03, 0.03);

TH2F h2_yx("h2_yx","hit y vs x; x [mm]; y [mm]",500,-1600.,1600.,500., -1600.,1600.);
TH2F h2_rz("h2_rz","hit #rho vs z; z [mm]; #rho [mm]",500,-2500.,2500.,500.,0.,1600.);

TH2F h2_yx_loc[6][7];
const int nlayer[6] = {5, 4, 3, 7, 3, 4};
const TString system_name[6] = { "VXD barrel", "VXD endcaps", "IT barrel", "IT endcaps", "OT barrel", "OT endcaps"};
const float limits[6] = {50.,50.,50.,100.,50.,160.};

void check_HTAtree(){

  // --- book histograms
  for (int isystem=0; isystem<6; ++isystem){
    for (int ilayer=0; ilayer<nlayer[isystem]; ++ilayer){
      TString hname = Form("h2_yx_loc_%d_%d",isystem,ilayer);
      TString htitle = Form("%s - layer %d; x_{loc} [mm]; y_{loc} [mm]", system_name[isystem].Data(), ilayer);
      h2_yx_loc[isystem][ilayer] = TH2F(hname, htitle, 100, -limits[isystem], limits[isystem], 100, -limits[isystem], limits[isystem]);
    }
  }
  
  // --- open the root file
  TFile *input_file = TFile::Open("ntu/ntu_muongun_pt3_theta10-170_phi0-30_20M.root");

  // --- retrieve the tree
  TTree *tree = (TTree*) input_file->Get("HTAtree");

  // --- define branch variables and set the branch addresses

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

  tree->SetBranchAddress("part_pdg", &part_pdg);
  tree->SetBranchAddress("part_px",  &part_px);
  tree->SetBranchAddress("part_py",  &part_py);
  tree->SetBranchAddress("part_pz",  &part_pz);
  tree->SetBranchAddress("part_e",   &part_e);
  tree->SetBranchAddress("part_vx",  &part_vx);
  tree->SetBranchAddress("part_vy",  &part_vy);
  tree->SetBranchAddress("part_vz",  &part_vz);
  tree->SetBranchAddress("part_t0",  &part_t0);
  tree->SetBranchAddress("part_q",   &part_q);

  // tracker hits
  static const int MAX_HITS = 1000;

  int n_hit;
  int hit_index[MAX_HITS];
  int hit_id0[MAX_HITS];
  float hit_x[MAX_HITS];
  float hit_y[MAX_HITS];
  float hit_z[MAX_HITS];
  float hit_t[MAX_HITS];
  float hit_xloc[MAX_HITS];
  float hit_yloc[MAX_HITS];

  tree->SetBranchAddress("n_hit", &n_hit);
  tree->SetBranchAddress("hit_index", hit_index);
  tree->SetBranchAddress("hit_id0", hit_id0);
  tree->SetBranchAddress("hit_x", hit_x);
  tree->SetBranchAddress("hit_y", hit_y);
  tree->SetBranchAddress("hit_z", hit_z);
  tree->SetBranchAddress("hit_t", hit_t);
  tree->SetBranchAddress("hit_xloc", hit_xloc);
  tree->SetBranchAddress("hit_yloc", hit_yloc);

  // --- loop over all entries
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t ientry = 0; ientry < nEntries; ++ientry) {

    if ( ientry % 100000 == 0 )
      std::cout << ientry << " / " << nEntries << std::endl;
    
    tree->GetEntry(ientry);
    
    // loop over the MC partiles
    const float part_pt = sqrt(part_px*part_px + part_py*part_py);
    const float part_invpt = part_q/part_pt;
    const float part_phi = atan2(part_py, part_px);
    const float part_theta = atan2(part_pt, part_pz);
    const float part_eta = -log(tan(0.5*part_theta));
    
    
    h_muon_pt.Fill(part_pt);
    h_muon_invpt.Fill(part_invpt);
    h_muon_phi.Fill(part_phi);
    h_muon_theta.Fill(part_theta);
    h_muon_eta.Fill(part_eta);
    h_muon_vz.Fill(part_vz);
    h_muon_t0.Fill(part_t0);
      

    // loop over the tracker hits
    for (int ihit=0; ihit<n_hit; ++ihit){

      // CellID encoding: "system:5,side:-2,layer:6,module:11,sensor:8"
      //  system:
      //    VXD barrel = 1
      //    VXD endcap = 2
      //    IT barrel  = 3
      //    IT endcap  = 4
      //    OT barrel  = 5
      //    OT endcap  = 6
      
      const unsigned int system = (unsigned) ( hit_id0[ihit] & 0x1f );
      const int side = (int) ( (hit_id0[ihit] >> 5) & 0x3 );
      unsigned int layer = (unsigned) ( (hit_id0[ihit] >> 7) & 0x3f );
      if (system == 1 || system == 2 ) layer /= 2;
      
      const float hit_rho = sqrt(hit_x[ihit]*hit_x[ihit]+hit_y[ihit]*hit_y[ihit]);

      h2_yx.Fill(hit_x[ihit], hit_y[ihit]);
      h2_rz.Fill(hit_z[ihit], hit_rho);

      h2_yx_loc[system-1][layer].Fill(hit_xloc[ihit], hit_yloc[ihit]);

      //cout << "    " << ihit << " " <<  system << " " << side << " " << layer << " "
      //	   << hit_x[ihit] << " " <<  hit_y[ihit] << " " <<  hit_z[ihit] << " " <<  (*hit_t)[ihit] << endl;
      
    } // ihit

    
  } // ientry loop
  
  
  // --- close the root file
  input_file->Close();

}
