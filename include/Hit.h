
#ifndef MDHT_HIT_H
#define MDHT_HIT_H

//
//  Hit.h
//  MuonColliderToy
//
//  Created by Luciano Ristori on 2/28/23
//  Modified for use with Muon Collider code April 3 2025
//

#include <iostream>

class Hit {

  public:
    long int ID;       // sequence number in hitList contained in event
    char     hitType;  // 'B' for barrel, 'D' for disc, 'G' for ghost (other codes possible)
    unsigned CellID;   // hit ID in the Muon Collider framework
    unsigned collInd;  // index of the hit collection in the Muon Collider framework
    unsigned hitInd;   // hit index in the collection of the Muon Collider framework
    unsigned layerInd; // single index going through all barrels and discs
    unsigned trackInd; // index of parent track in track list (member of Event)

    double x1; // primary local coordinate [mm]
    double x2; // secondary local coordinate [mm]
    double t;  // time coordinate [mm]

    double x; // global coordinate [mm]
    double y; // global coordinate [mm]
    double z; // global coordinate [mm]

    double u1x; // components of x1 unit vector
    double u1y;
    double u1z;

    double u2x; // components of x2 unit vector
    double u2y;
    double u2z;

    double errx1; // error on u local coortinate
    double errx2; // errors on v local coordinate

    double u0, u1, u2; // diagonalized local coordinates

    Hit(char hitType, unsigned CellID, unsigned collInd, unsigned hitInd, unsigned layerInd, unsigned trackInd,
        double x1, double x2, double t, double x, double y, double z, double u1x, double u1y, double u1z, double u2x,
        double u2y, double u2z, double errx1, double errx2);

    Hit(char hitType, unsigned CellID, unsigned layerInd, unsigned trackInd, double x1, double x2, double t, double x,
        double y, double z, double u1x, double u1y, double u1z, double u2x, double u2y, double u2z, double errx1,
        double errx2);

    Hit(char hitType);

    void print(std::ostream& out);
    void write(std::ostream& out);

    bool isSeed();
};

#endif //MDHT_HIT_H
