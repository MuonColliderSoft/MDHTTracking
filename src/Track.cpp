
#ifndef MDHT_TRACK_CPP
#define MDHT_TRACK_CPP

//
//  Track.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 6/25/21
//  Modified for use with Muon Collider code April 3 2025
//

#include <cmath>

#include "include/Track.h"
#include "include/Hit.h"
#include "include/HTAmapper.h"

// constructor

Track::Track(double mass_, double x0_, double y0_, double z0_, double t0_, double invPt_, double eta_, double phi_) {

    init(mass_, x0_, y0_, z0_, t0_, invPt_, eta_, phi_); // initialize track
};

void Track::print(std::ostream& out, int mode = 0) {

    // mode = 0: track parameters only
    // mode = 1: list all hits

    int i, j, k;
    trackToCell(i, j, k, *this);

    out << "Track ID: " << ID << " mass: " << mass << " beta: " << beta << " vertex: x: " << x0 << " y: " << y0
        << " z: " << z0 << " t: " << t0 << " charge: " << charge << " invPt: " << invPt << " eta: " << eta
        << " phi: " << phi << " [" << i << "," << j << "," << k << "]" << std::endl;

    if (mode == 0)
        return;

    out << hitList.size() << " hits" << std::endl;
    for (unsigned iH = 0; iH != hitList.size(); ++iH) {
        out << "hit " << iH;
        (hitList[iH]).print(out);
    }

}; // end print

void Track::write(std::ostream& out) { // write track to a file

    out << ID << " " << mass << " " << beta << " " << x0 << " " << y0 << " " << z0 << " " << t0 << " " << charge << " "
        << Pt << " " << eta << " " << phi << std::endl;

    for (unsigned iH = 0; iH != hitList.size(); ++iH) {
        out << ID << " " << iH << " ";
        (hitList[iH]).write(out);
    }

}; // end write

void Track::init(double mass_, double x0_, double y0_, double z0_, double t0_, double invPt_, double eta_,
                 double phi_) {

    hitList.clear();

    // coordinates of primary vertex

    x0 = x0_;
    y0 = y0_;
    z0 = z0_;
    t0 = t0_;

    // primary parameters from constructor arguments

    mass  = mass_;
    invPt = invPt_;
    eta   = eta_;
    phi   = phi_;

    // derived parameters

    //c = 3.e-4*invPt*par.magneticField; // signed curvature in mmm^(-1)

    charge = 1.;
    if (invPt < 0.)
        charge = -1.;

    cotTheta = sinh(eta); // cotTheta = pz/pt
    tgTheta  = 1. / cotTheta;

    Pt = 1. / invPt; // transverse momentum
    if (Pt < 0.)
        Pt = -Pt;
    Pz = Pt * cotTheta; // longitudinal momentum

    double P2 = Pt * Pt + Pz * Pz;
    double E2 = P2 + mass * mass;
    beta      = sqrt(P2 / E2); // velocity

    P = sqrt(P2);

    cosTheta = Pz / sqrt(P2);

    phi0 = phi - c * Pt * z0 / 2. / Pz;

}; // end init

//void Track::getIJK(unsigned &I, unsigned &J, unsigned &K){

//I = phi_to_xi(phi);
//J = eta_to_xj(eta);
//K = invPt_to_xk(invPt);

//return;

//}; // returns indices in HTM array

#endif // MDHT_TRACK_CPP
