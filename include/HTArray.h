
#ifndef MDHT_HTARRAY_CPP
#define MDHT_HTARRAY_CPP

//
//  HTArray.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 5/23/22
//

#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "TH1D.h"
#include "TH2I.h"
#include "TH3I.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "include/TrackReader.h"
#include "include/Event.h"
#include "include/Hit.h"
#include "include/Track.h"

#include "include/Parameters.h"
#include <limits>

/// global variables

extern int  TrainingPhase;
extern bool Special;
extern bool verbose;

extern Parameters par;

double xi_to_phi(double xi) {
    return par.HTA_t_phi - par.HTA_t_deltaPhi + (2 * par.HTA_t_deltaPhi / (double)par.HTA_NphiBins) * xi;
}
double xj_to_eta(double xj) {
    return par.HTA_t_eta - par.HTA_t_deltaEta + 2 * (par.HTA_t_deltaEta / (double)par.HTA_NetaBins) * xj;
}
double xk_to_invPt(double xk) {
    return par.HTA_t_invPt_min + (par.HTA_t_invPt_max - par.HTA_t_invPt_min) / (double)par.HTA_NinvptBins * xk;
}

double phi_to_xi(double phi) {
    return (phi - par.HTA_t_phi + par.HTA_t_deltaPhi) / par.HTA_t_deltaPhi * (double)par.HTA_NphiBins / 2.;
}
double eta_to_xj(double eta) {
    return (eta - par.HTA_t_eta + par.HTA_t_deltaEta) / par.HTA_t_deltaEta * (double)par.HTA_NetaBins / 2.;
}
double invPt_to_xk(double invPt) {
    return (invPt - par.HTA_t_invPt_min) / (par.HTA_t_invPt_max - par.HTA_t_invPt_min) * (double)par.HTA_NinvptBins;
}

struct HelixPars {
    double phi, theta, R, z0, t0, beta;
};

struct PointCoordinates {
    double X, Y, Z, T;
};

// Helix function. time is time from crossing, result is X,Y,Z,
//					track origin is assumed to be at (x=0, y=0, z=z0),
//					t0 is time of track at origin

PointCoordinates Helix(double time, HelixPars pars) {

    double s     = (time - pars.t0) * pars.beta;
    double phi   = pars.phi;
    double theta = pars.theta;
    double R     = pars.R;
    double z0    = pars.z0;
    // pars.t0 and pars.beta are used in s; silence unused warnings if any
    (void)pars.t0;
    (void)pars.beta;

    double sinTheta = sin(theta);
    double cosTheta = cos(theta);
    double cosPhi   = cos(phi);
    double sinPhi   = sin(phi);

    PointCoordinates result;

    result.X = R * (1 - cos((s * sinTheta) / R)) * sinPhi + R * cosPhi * sin((s * sinTheta) / R);
    result.Y = -(R * cosPhi * (1 - cos((s * sinTheta) / R))) + R * sinPhi * sin((s * sinTheta) / R);
    result.Z = z0 + s * cosTheta;
    result.T = time;

    return result;
};

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// global variables for track fitting

std::vector<Hit>                   fitHitList;   // list of all hits to be used in one fit
std::vector<std::vector<long int>> combinations; // all hit combinations to be fitted in one candidate

////////////////////////////////////////////////////////////////////////////////////////
// global function for track fitting
// Gets track parameters in input and returns chi2 using fitHitList as the list of hits

double chi2Function(const double* x) {

    double invPt = xk_to_invPt(x[0]);
    double eta   = xj_to_eta(x[1]);
    double phi   = xi_to_phi(x[2]);
    double z0    = x[3];
    double t0    = x[4];
    double beta  = x[5];

    // need some parameter transformation

    double       theta       = 2. * atan(exp(-eta));               // eta = -log(tan(theta/2))
    const double RConvFactor = 0.000299792458 * par.magneticField; // conversion factor from invPt to R
    double       R           = 1. / RConvFactor / invPt;

    //const float MagneticField = 5.; // [Tesla]
    //const float K = 0.000299792458*MagneticField;

    HelixPars hPars;

    hPars.phi   = phi;
    hPars.theta = theta;
    hPars.R     = R;
    hPars.z0    = z0;
    hPars.t0    = t0;
    hPars.beta  = beta;

    std::vector<double> fitTimes;
    unsigned            nPoints = fitHitList.size();
    for (unsigned i = 0; i < nPoints; ++i)
        fitTimes.push_back(x[6 + i]);

    PointCoordinates pointXYZ;

    double chi2 = 0.;

    for (unsigned iPoint = 0; iPoint != nPoints; ++iPoint) {

        pointXYZ = Helix(fitTimes[iPoint], hPars);

        Hit thisHit = fitHitList[iPoint];

        double deltax = pointXYZ.X - thisHit.x;
        double deltay = pointXYZ.Y - thisHit.y;
        double deltaz = pointXYZ.Z - thisHit.z;
        double deltat = pointXYZ.T - thisHit.t;

        double u1x = thisHit.u1x;
        double u1y = thisHit.u1y;
        double u1z = thisHit.u1z;

        double u2x = thisHit.u2x;
        double u2y = thisHit.u2y;
        double u2z = thisHit.u2z;

        double errx1 = thisHit.errx1;
        double errx2 = thisHit.errx2;

        const double errt  = 18.;   // 18 mm = 60 ps
        const double errx3 = 0.001; // 1µm error perpendicular to the sensor plane

        // find u3 by vector product of u1 and u2

        double u3x = u1y * u2z - u1z * u2y;
        double u3y = u1z * u2x - u1x * u2z;
        double u3z = u1x * u2y - u1y * u2x;

        double deltax1 = deltax * u1x + deltay * u1y + deltaz * u1z; // deltax * u1 scalar product
        double deltax2 = deltax * u2x + deltay * u2y + deltaz * u2z; // deltax * u2 scalar product
        double deltax3 = deltax * u3x + deltay * u3y + deltaz * u3z; // deltax * u3 scalar product

        chi2 += (deltax1 / errx1) * (deltax1 / errx1);
        chi2 += (deltax2 / errx2) * (deltax2 / errx2);
        chi2 += (deltax3 / errx3) * (deltax3 / errx3);
        chi2 += (deltat / errt) * (deltat / errt);
        /*	
		std::cout << " " << (deltax1/errx1)*(deltax1/errx1);
		std::cout << " " << (deltax2/errx2)*(deltax2/errx2);	
		std::cout << " " << (deltax3/errx3)*(deltax3/errx3);
		std::cout << " " << (deltat/errt)*(deltat/errt);
		std::cout << std::endl;
	*/
    }

    //std::cout << "chi2 = " << chi2 << std::endl;

    return chi2;

} // end chi2Functions

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

class Hitstat {

    // This is the class that contains all that is needed for each layer in each HT cell

  public:
    unsigned nHits;

    long int nEntries;
    double   mean[3];                                  // mean values of x1, x2 and t
    double   low[3];                                   // min value
    double   high[3];                                  // max value
    double   c00, c11, c22, c01, c02, c12;             //  covariance matrix c[i,j]
    double   rotAngle;                                 // rotation angle for x2-t diagonalization
    double   sinAngle;                                 // sin of rotAngle
    double   cosAngle;                                 // cos of rotAngle
    double   u0, u1, u2;                               // diagonal coordinates
    double   u0m, u1m, u2m;                            // mean of diagonal coordinates
    double   u0v, u1v, u2v;                            // variance of diagonal coordinates
    double   u0Max, u0Min, u1Max, u1Min, u2Max, u2Min; // min-max values

    Hitstat() { // constructor
        nEntries = 0;
        nHits    = 0;
        //	hitIDList.clear();
        c00 = 0.;
        c11 = 0.;
        c22 = 0.;
        c01 = 0.;
        c02 = 0.;
        c12 = 0.;
        for (int k = 0; k != 3; ++k) {
            mean[k] = 0.;
            low[k]  = 0.;
            high[k] = 0.;
        }
        rotAngle = 0.;
        sinAngle = sin(rotAngle);
        cosAngle = cos(rotAngle);
        u0       = 0.;
        u1       = 0.;
        u2       = 0.;
        u0m      = 0.;
        u1m      = 0.;
        u2m      = 0.;
        u0v      = 0.;
        u1v      = 0.;
        u2v      = 0.;
        u0Max    = 0.;
        u0Min    = 0.;
        u1Max    = 0.;
        u1Min    = 0.;
        u2Max    = 0.;
        u2Min    = 0.;
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    void train(Hit& h) {

        // incremental calculation of means, variances and covariance.
        // classic algorithm

        ++nEntries;

        if (TrainingPhase == 1) {

            mean[0] += (h.x1 - mean[0]) / nEntries;
            mean[1] += (h.x2 - mean[1]) / nEntries;
            mean[2] += (h.t - mean[2]) / nEntries;

            if (nEntries > 1) {
                c00 += (h.x1 - mean[0]) * (h.x1 - mean[0]) / (nEntries - 1) - c00 / nEntries;
                c11 += (h.x2 - mean[1]) * (h.x2 - mean[1]) / (nEntries - 1) - c11 / nEntries;
                c22 += (h.t - mean[2]) * (h.t - mean[2]) / (nEntries - 1) - c22 / nEntries;
                c01 += (h.x1 - mean[0]) * (h.x2 - mean[1]) / (nEntries - 1) - c01 / nEntries;
                c02 += (h.x1 - mean[0]) * (h.t - mean[2]) / (nEntries - 1) - c02 / nEntries;
                c12 += (h.x2 - mean[1]) * (h.t - mean[2]) / (nEntries - 1) - c12 / nEntries;
            }

            if (nEntries == 1) {
                low[0]  = h.x1;
                high[0] = h.x1;
                low[1]  = h.x2;
                high[1] = h.x2;
                low[2]  = h.t;
                high[2] = h.t;
            } else {
                if (h.x1 > high[0])
                    high[0] = h.x1;
                if (h.x1 < low[0])
                    low[0] = h.x1;
                if (h.x2 > high[1])
                    high[1] = h.x2;
                if (h.x2 < low[1])
                    low[1] = h.x2;
                if (h.t > high[2])
                    high[2] = h.t;
                if (h.t < low[2])
                    low[2] = h.t;
            }

        } // end if(TrainingPhase == 1)

        if (TrainingPhase >= 2) {
            // rotation into diagonal space u0,u1,u2

            u0 = h.x1 - mean[0];
            u1 = (h.x2 - mean[1]) * cosAngle - (h.t - mean[2]) * sinAngle;
            u2 = (h.x2 - mean[1]) * sinAngle + (h.t - mean[2]) * cosAngle;

            // copy vector u into hit

            h.u0 = u0;
            h.u1 = u1;
            h.u2 = u2;

        } // end if(TrainingPhase >= 2)

        if (TrainingPhase == 2) {

            // find min and max values of u's

            if (nEntries == 1) {
                u0Max = u0;
                u0Min = u0;
                u1Max = u1;
                u1Min = u1;
                u2Max = u2;
                u2Min = u2;
            } else {
                if (u0 > u0Max)
                    u0Max = u0;
                if (u0 < u0Min)
                    u0Min = u0;
                if (u1 > u1Max)
                    u1Max = u1;
                if (u1 < u1Min)
                    u1Min = u1;
                if (u2 > u2Max)
                    u2Max = u2;
                if (u2 < u2Min)
                    u2Min = u2;
            }

            // incremental calculation of mean and variance of diagonal coordinates u0,u1,u2

            u0m += (u0 - u0m) / nEntries;
            u1m += (u1 - u1m) / nEntries;
            u2m += (u2 - u2m) / nEntries;

            if (nEntries > 1) {
                u0v += (u0 - u0m) * (u0 - u0m) / (nEntries - 1) - u0v / nEntries;
                u1v += (u1 - u1m) * (u1 - u1m) / (nEntries - 1) - u1v / nEntries;
                u2v += (u2 - u2m) * (u2 - u2m) / (nEntries - 1) - u2v / nEntries;
            }
        } // end if(TrainingPhase == 2)

    } // end train

    ////////////////////////////////////////////////////////////////////////////////////////

    void diagonalize() {

        // finds diagonalizing rotation angle
        // for this cell, this layer

        rotAngle = -0.5 * atan2(2 * c12, c11 - c22);
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    long int GetEntries() { return nEntries; }

    ////////////////////////////////////////////////////////////////////////////////////////

    double GetMean(int iPar) { return mean[iPar]; }

    ////////////////////////////////////////////////////////////////////////////////////////

    // not clear where this was needed (LR)

    //void reset(){

    //	hitLayer = false;
    //	hitIDList.clear();
    //	nHits = 0;
    //}

    ////////////////////////////////////////////////////////////////////////////////////////

    void print(std::ostream& out) {

        out << "N: " << nEntries << " x1: " << mean[0] << " x2: " << mean[1] << " t: " << mean[2] << " ";
        //out << "c00: " << c00 << " c11: " << c11 << " c22: " << c22 << " c01: " << c01 << " c02: " << c02 << " c12: " << c12;
        out << " rotAngle: " << rotAngle;
        out << "     u0:[" << u0Min << "," << u0Max << "]; u1:[" << u1Min << "," << u1Max << "]; u2:[" << u2Min << ","
            << u2Max << "]";
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    // are these needed? (LR)

    /* 	void printHits(std::ostream &out){
    	
    		int nHits = hitIDList.size();
    		out << "nHits = " << nHits << ": ";
    		for(int iH = 0; iH != nHits; ++ iH) out << hitIDList[iH] << " ";
    		out << std::endl;
    	}
   */

    /* 	void writeHits(std::ostream &out, std::vector<Hit> & hitList){
    	
    		int nHits = hitIDList.size();
    		if(nHits) out << nHits << std::endl;
    		for(int iH = 0; iH != nHits; ++ iH) {
    			hitList[hitIDList[iH]].write(out);
    		}
    */

    ////////////////////////////////////////////////////////////////////////////////////////

    void write(std::ostream& out) {

        out << mean[0] << " " << mean[1] << " " << mean[2] << " ";
        out << low[0] << " " << low[1] << " " << low[2] << " ";
        out << high[0] << " " << high[1] << " " << high[2] << " ";
        out << rotAngle << " " << u0Min << " " << u1Min << " " << u2Min << " " << u0Max << " " << u1Max << " " << u2Max
            << " ";
    }

}; // end Hitstat

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

int decodeLayer(int CellID) {

    const unsigned int system = (unsigned)(CellID & 0x1f);
    const unsigned int side   = (int)((CellID >> 5) & 0x3);
    const unsigned int layer  = (unsigned)((CellID >> 7) & 0x3f);

    return 100 * system + 10 * side + layer;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

class HTArrayElement {

    // This is one element of the HT Array with parameters, attached hits and accessories

  public:
    // pointer to the vector holding all hits for this event
    // all hit manipulation done in this class is done using indices pointing
    // in this vector.
    // indices are substituted by hits before returning the results
    // (e.g. fitCandidate)

    //std::vector<Hit> * allHits;

    // array indices and mean values of parameters: initialized by constructor

    int i_index;
    int j_index;
    int k_index;

    double imeanPhi;
    double jmeanEta;
    double kmeanInvPt;

    bool thisCellDone = false; //used by algorithm to avoid duplicates

    // list of sensors for this element implemented as a map
    // each sensor (encoded as an int), is linked to its own   object with training hit statistics and derived cuts
    // coordinates here are local

    std::map<int, Hitstat> layerIndHitStat;

    // This map links each layer (layerInd) to a vector of hit IDs
    // Those are the hits that have been attached to this HTA element.
    // The ID's point into hitList member of Event.

    std::map<int, std::vector<long int>> layerIndHitIDList;

    unsigned nHitLayers = 0; // number of different layers hit for current event in this element

    // next two params calculated in training phase and loaded at initialization
    unsigned minLayers = 10000; // min number of layers hit by a track from this param space
    unsigned maxLayers = 0;     // max number of layers hit by a track from this param space

    ////////////////////////////////////////////////////////////////////////////////////////

    void train(Hit& h) {

        if (TrainingPhase == 1)
            layerIndHitStat[h.CellID].train(h);
        else if (layerIndHitStat.find(h.CellID) != layerIndHitStat.end())
            layerIndHitStat[h.CellID].train(h);
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    // diagonalize this element

    void diagonalize() {

        for (auto it = layerIndHitStat.begin(); it != layerIndHitStat.end(); ++it) {
            (it->second).diagonalize();
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    // reset all hit lists in all layers of this element

    void reset() {

        thisCellDone = false;
        nHitLayers   = 0;
        for (auto it = layerIndHitIDList.begin(); it != layerIndHitIDList.end(); ++it) {
            (it->second).clear();
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    void print(std::ostream& out) {
        std::vector<std::map<int, Hitstat>::iterator> itList;
        out << "min,max Layers: " << minLayers << "," << maxLayers << std::endl;
        for (auto iter = layerIndHitStat.begin(); iter != layerIndHitStat.end(); ++iter) {
            itList.push_back(iter);
        }

        // Sort the vector in descending order using a lambda
        std::sort(itList.begin(), itList.end(),
                  [](std::map<int, Hitstat>::iterator a, std::map<int, Hitstat>::iterator b) {
                      return decodeLayer(a->first) < decodeLayer(b->first);
                  });

        std::cout << std::endl;
        for (size_t i = 0; i < itList.size(); ++i) {
            auto it2 = itList[i];
            out << "   sensor: " << it2->first << " layer: " << decodeLayer(it2->first) << " ";
            (it2->second).print(out);
            out << std::endl;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    /*     
        void printHits(std::ostream &out){
        	out << "min,max Layers: " << minLayers << "," << maxLayers << std::endl;
        	for(it = layerIndHitStat.begin(); it != layerIndHitStat.end(); ++it){
        		out << "lay: " << it->first; out << " ";
        			(it->second).printHits(out); out << std::endl;
        	}
        }
        
   */

    ////////////////////////////////////////////////////////////////////////////////////////

    // prints all the hits sorted by layer

    void printHits(std::ostream& out) {
        for (auto it = layerIndHitIDList.begin(); it != layerIndHitIDList.end(); ++it) {
            out << "layerInd: " << it->first << " hit IDs : ";
            const std::vector<long int>& hitIDList = it->second;
            for (size_t iH = 0; iH < hitIDList.size(); ++iH) {
                long int hitID = hitIDList[iH];
                out << hitID << ", ";
            }
            out << std::endl;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    /*        
        
          void writeHits(std::ostream &out, std::vector<Hit> & hitList){
        	for(it = layerIndHitStat.begin(); it != layerIndHitStat.end(); ++it){
        			(it->second).writeHits(out, hitList);// out << std::endl;
        	}
        }
        
    */

    ////////////////////////////////////////////////////////////////////////////////////////

    int fill(Hit& h) {

        auto it_loc = layerIndHitStat.find(h.CellID);
        if (it_loc != layerIndHitStat.end()) {
            Hitstat stat = it_loc->second;

            double u0 = h.x1 - stat.mean[0];
            if (Special)
                std::cout << u0 << " " << stat.u0Min << " " << stat.u0Max << std::endl;
            if (u0 < stat.u0Min || u0 > stat.u0Max) {
                Special = false;
                return -1;
            }
            if (Special)
                std::cout << "good u0" << std::endl;
            double u1 = (h.x2 - stat.mean[1]) * stat.cosAngle - (h.t - stat.mean[2]) * stat.sinAngle;
            if (Special)
                std::cout << u1 << " " << stat.u1Min << " " << stat.u2Max << std::endl;
            if (u1 < stat.u1Min || u1 > stat.u1Max) {
                Special = false;
                return -2;
            }
            if (Special)
                std::cout << "good u1" << std::endl;
            double u2 = (h.x2 - stat.mean[1]) * stat.sinAngle + (h.t - stat.mean[2]) * stat.cosAngle;
            if (Special)
                std::cout << u2 << " " << stat.u2Min << " " << stat.u2Max << std::endl;
            if (u2 < stat.u2Min || u2 > stat.u2Max) {
                Special = false;
                return -3;
            }
            if (Special)
                std::cout << "good u2" << std::endl;

            const unsigned int layerInd = CellIDtoLayer(h.CellID);

            //if(hitList[layer].size() == 0) ++nHitLayers;
            //hitIDList[layer].push_back(h.ID);

            layerIndHitIDList[layerInd].push_back(h.ID);
            if ((int)layerIndHitIDList[layerInd].size() == 1)
                ++nHitLayers;

            Special = false;
            return 0;
        } else
            return -4;

    } // end fill

    ////////////////////////////////////////////////////////////////////////////////////////

    static bool compareComb(std::vector<long int> v1, std::vector<long int> v2) { return v1.size() > v2.size(); }

    ////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////
    // This array element has been selected as a candidate and all the hit combinations
    // must be formed to find the best track fit. A maximum of one track is chosen
    // on the basis of best chi square and number of hits. Combinations of hits
    // where hits are discarded are allowed provided the minimum number of hits
    // required for a fit is satisfied. A combination of lesser number of hits
    // is considered only if no combination with more hits has an acceptable chi
    // square.                           Luciano Feb 12, 2024
    //////////////////////////////////////////////////////////////////////////////////

    int fitCandidate(std::vector<Hit>* allHits, double& chi2, double& phi, double& eta, double& invPt, double& z0,
                     double& t0, double& beta, unsigned& nLayers, std::vector<long int>& goodFitHitList) {

        if (verbose)
            std::cout << "fitCandidate " << i_index << " " << j_index << " " << k_index << std::endl;

        struct Result {
            double   chi2;
            double   phi;
            double   eta;
            double   invPt;
            double   z0;
            double   t0;
            double   beta;
            unsigned nLayers;
        };

    Result bestFitRes;
    // initialize to safe defaults
    bestFitRes.chi2  = std::numeric_limits<double>::infinity();
    bestFitRes.phi   = 0.;
    bestFitRes.eta   = 0.;
    bestFitRes.invPt = 0.;
    bestFitRes.z0    = 0.;
    bestFitRes.t0    = 0.;
    bestFitRes.beta  = 0.;
    bestFitRes.nLayers = 0;
        int    bestFitComb = 0;

        // create list of combination of hits to be fitted

        std::vector<std::vector<long int>>  combinations2; // all hit combinations to be fitted in one candidate
        std::vector<std::vector<long int>>* comb1 = &combinations;
        std::vector<std::vector<long int>>* comb2 = &combinations2;
        std::vector<std::vector<long int>>* temp;
        std::vector<long int>               hitList;

        // create one single empty combination

        hitList.clear();
        comb1->clear();
        comb2->clear();
        comb1->push_back(hitList);

        // iterate on all layers

        //std::map<int, std::vector<long int>>::iterator it;
        for (auto it = layerIndHitIDList.begin(); it != layerIndHitIDList.end(); ++it) {

            std::vector<long int>* hitIDList = &(it->second);

            int nHits = hitIDList->size();
            if (nHits == 0)
                continue; // skip possible empty layers

            //define ghost hit ID

            long int ghostHitID = -1;

            // iterate on all previous combinations
            int nComb = comb1->size();
            for (int iC = 0; iC != nComb; ++iC) {

                // iterate on all hits in this layer

                for (int ih = 0; ih != nHits; ++ih) {
                    long int hitID = (*hitIDList)[ih];
                    hitList        = (*comb1)[iC];
                    hitList.push_back(hitID);
                    comb2->push_back(hitList);
                }
                // add ghost hit

                hitList = (*comb1)[iC];
                hitList.push_back(ghostHitID);
                comb2->push_back(hitList);

                (*comb1)[iC].clear();
            }

            comb1->clear();
            temp  = comb1;
            comb1 = comb2;
            comb2 = temp;
        }

        combinations = *comb1; // including ghost hits

        int nCombs = combinations.size();

        std::vector<std::vector<long int>> finalCombinations;
        std::vector<long int>              tempHitList;

        for (int iComb = 0; iComb != nCombs; ++iComb) {

            tempHitList.clear();

            // discard ghost hits from combination
            for (int i = 0; i != (int)combinations[iComb].size(); ++i) {
                long int thisHitID = combinations[iComb][i];
                if (thisHitID != -1)
                    tempHitList.push_back(thisHitID);
            }

            // check for minimum number of hits	in this combination

            unsigned minHits = std::max(nHitLayers - par.reco_maxDroppedLayers, par.reco_minLayersForFit);

            if (tempHitList.size() >= minHits)
                finalCombinations.push_back(tempHitList);
        }

        nCombs = finalCombinations.size(); // final unsorted list

        sort(finalCombinations.begin(), finalCombinations.end(), compareComb);

        if (verbose) {
            std::cout << nCombs << " final combinations:" << std::endl;
            for (int i = 0; i != nCombs; ++i) {
                std::cout << finalCombinations[i].size() << " ";
            }
            std::cout << std::endl;
        }

        // Minimize with Root Minimizer

        // create minimizer giving a name and a name (optionally) for the specific
        // algorithm
        // possible choices are:
        //     minName                  algoName
        // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
        //  Minuit2                     Fumili2
        //  Fumili
        //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
        //                              BFGS2, SteepestDescent
        //  GSLMultiFit
        //   GSLSimAn xxx problems
        //   Genetic

        const char* minName  = "Minuit2";
        const char* algoName = "";

        ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

        ////////////////////////////////////////////////////////////////
        // loop on all combinations for this candidate

        bool oneGoodFit   = false; // did we find at least one good hit combination?
        bool firstGoodFit = true;  // is this the first good fit for min chi2?

        for (int iComb = 0; iComb != nCombs; ++iComb) {

            // "fitHitList" is the list of hits that will be fitted for this combination
            // retrive hits from allHits using hit index

            fitHitList.clear();
            for (int iH = 0; iH != (int)finalCombinations[iComb].size(); ++iH) {
                fitHitList.push_back((*allHits)[finalCombinations[iComb][iH]]);
            }

            // If we already have a good fit from a previous combination
            // and the number of hits in this combination is less that
            // the number of hits in the combination with the best fit,
            // we end the loop and ignore all the remaining combinations.
            // Note that combinations are sorted by number of hits.

            if (oneGoodFit && (fitHitList.size() < bestFitRes.nLayers))
                break;

            if (verbose) {
                std::cout << std::endl;
                std::cout << "Combination " << iComb << " ";
                std::cout << fitHitList.size() << " Hits" << std::endl;
                //for(unsigned iH = 0; iH != fitHitList.size(); ++iH){
                //std::cout << fitHitList[iH].trackInd << " ";
                //}
                std::cout << "Fit hit list: " << std::endl;
                for (unsigned iH = 0; iH != fitHitList.size(); ++iH) {
                    fitHitList[iH].print(std::cout);
                }
            }

            // set tolerance , etc...
            min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
            min->SetMaxIterations(1000);       // for GSL
            min->SetTolerance(0.001);
            min->SetPrintLevel(0);

            // create function wrapper for minimizer
            // a IMultiGenFunction type

            const unsigned      maxPoints   = 20;
            const unsigned      nMaxFitPars = 6 + maxPoints;
            ROOT::Math::Functor f(&chi2Function, nMaxFitPars);
            double              step[nMaxFitPars];
            double              variable[nMaxFitPars];

            // 6 TRack parameters

            step[0] = 0.1; // invPt
            step[1] = 0.1; // eta
            step[2] = 0.1; // phi
            step[3] = 0.1; // z0
            step[4] = 0.1; // t0
            step[5] = 0.1; // beta

            variable[0] = kmeanInvPt; // invPt
            variable[1] = jmeanEta;   // eta
            variable[2] = imeanPhi;   // phi
            variable[3] = 0.;         // z0
            variable[4] = 0.;         // t0
            variable[5] = 1.;         // beta

            min->SetFunction(f);

            // Set the free variables to be minimized!

            min->SetVariable(0, "invPt", variable[0], step[0]);
            min->SetVariable(1, "Eta", variable[1], step[1]);
            min->SetVariable(2, "Phi", variable[2], step[2]);
            min->SetVariable(3, "z0", variable[3], step[3]);
            min->SetVariable(4, "t0", variable[4], step[4]);
            min->SetVariable(5, "beta", variable[5], step[5]);

            unsigned nHits = fitHitList.size();
            if (nHits > maxPoints) {
                std::cout << "Too many Hits" << std::endl;
                return -2;
            }

            unsigned totFitPars = 6 + nHits; // number of paramers that need to be fitted

            // loop on time parameters

            for (int iFitPar = 6; iFitPar != (int)totFitPars; ++iFitPar) {

                min->SetVariable(iFitPar, ("TOF" + std::to_string(iFitPar - 6)).c_str(), fitHitList[iFitPar - 6].t,
                                 0.1);
            };

            // loop on unused parameters and fix them

            for (int iFitPar = totFitPars; iFitPar != nMaxFitPars; ++iFitPar) {

                variable[iFitPar] = 0.;
                min->SetVariable(iFitPar, "UNUSED", 0., 0.1);
                min->FixVariable(iFitPar);
            };

            min->Minimize(); // calls Minuit

            // extract results returned by the fit

            const double* res = min->X();

            chi2  = min->MinValue();
            eta   = xj_to_eta(res[1]);
            phi   = xi_to_phi(res[2]);
            invPt = xk_to_invPt(res[0]);
            z0    = res[3];
            t0    = res[4];
            beta  = res[5];

            if (chi2 > par.reco_chi2Cut)
                continue; // failed chi2

            oneGoodFit = true;

            if (firstGoodFit || bestFitRes.chi2 > chi2) {
                firstGoodFit       = false;
                bestFitRes.chi2    = chi2;
                bestFitRes.phi     = phi;
                bestFitRes.eta     = eta;
                bestFitRes.invPt   = invPt;
                bestFitRes.z0      = z0;
                bestFitRes.t0      = t0;
                bestFitRes.beta    = beta;
                bestFitRes.nLayers = fitHitList.size();
                bestFitComb        = iComb;
            }

        } // end loop on final combinations

        if (oneGoodFit) {
            nLayers = finalCombinations[bestFitComb].size();

            chi2           = bestFitRes.chi2;
            phi            = bestFitRes.phi;
            eta            = bestFitRes.eta;
            invPt          = bestFitRes.invPt;
            z0             = bestFitRes.z0;
            t0             = bestFitRes.t0;
            beta           = bestFitRes.beta;
            goodFitHitList = finalCombinations[bestFitComb];

            if (verbose)
                std::cout << "bestFitComb: " << bestFitComb << " chi2: " << chi2 << std::endl;

            return 0;
        } else
            return -1;

    } // end fitCandidate

    ////////////////////////////////////////////////////////////////////////////////////////

}; // end HTArrayElement

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

class HTArray {

    // This is the whole Hough Transform Array with parameters and accessories

  public:
    // This vector holds all  hits for this event
    // It is filled by the "fill" method of this class
    // It is cleared by the "reset" method of this class
    // Hits in this vector are referred by index

    std::vector<Hit> allHits;

    // some histograms of inside quantities

    TH1D* HDeltaEta;
    TH1D* HDeltaEtaCorrected;
    TH1D* HDeltaPhi;
    TH1D* HDeltaPhi2;
    TH1D* HInvptTimesR;
    TH2I* HDeltaPhiVsInvptTimesR;

    // histograms of a single HT cell content as an example

    TH3I* H3D_HTcellx;
    TH3I* H3D_HTcellu;
    TH3I* H3D_HTAtrainStat;
    TH1D* HstatPhi;
    TH1D* HstatEta;
    TH1D* HstatInvpt;

    TH1D* Hcellx1;
    TH1D* Hcellx2;
    TH1D* Hcellt;
    TH1D* Hcellu0;
    TH1D* Hcellu1;
    TH1D* Hcellu2;

    // dimensions of the array

    const static unsigned NphiBins   = par.HTA_NphiBins;
    const static unsigned NetaBins   = par.HTA_NetaBins;
    const static unsigned NinvptBins = par.HTA_NinvptBins;

    // origin and steps of the array

    double phiMin;
    double phiStep;
    double etaMin;
    double etaStep;
    double invPtMin;
    double invPtStep;

    // HT Array bin and layer to plot content in 3D

    unsigned plotBinX = par.HTA_plotBinX;
    unsigned plotBinY = par.HTA_plotBinY;
    unsigned plotBinZ = par.HTA_plotBinZ;
    unsigned plotLay  = par.HTA_plotLay;

    // HT Array proper

    HTArrayElement ArrElem[NphiBins][NetaBins][NinvptBins];

    double signedCurvaturePerInvPt;

    ////////////////////////////////////////////////////////////////////////////////////////

    double getCell_Phi(int i) { return phiMin + ((double)i + 0.5) * phiStep; }

    ////////////////////////////////////////////////////////////////////////////////////////

    double getCell_Eta(int j) { return etaMin + ((double)j + 0.5) * etaStep; }

    ////////////////////////////////////////////////////////////////////////////////////////

    double getCell_InvPt(int k) { return invPtMin + ((double)k + 0.5) * invPtStep; }

    ////////////////////////////////////////////////////////////////////////////////////////

    double get_Phi(double xi) { return phiMin + xi * phiStep; }

    ////////////////////////////////////////////////////////////////////////////////////////

    double get_Eta(double xj) { return etaMin + xj * etaStep; }

    ////////////////////////////////////////////////////////////////////////////////////////

    double get_InvPt(double xk) { return invPtMin + xk * invPtStep; }

    ////////////////////////////////////////////////////////////////////////////////////////

    HTArray() { // constructor

        phiMin    = par.HTA_t_phi - par.HTA_t_deltaPhi;
        phiStep   = 2 * par.HTA_t_deltaPhi / double(NphiBins);
        etaMin    = par.HTA_t_eta - par.HTA_t_deltaEta;
        etaStep   = 2 * par.HTA_t_deltaEta / double(NetaBins);
        invPtMin  = par.HTA_t_invPt_min;
        invPtStep = (par.HTA_t_invPt_max - par.HTA_t_invPt_min) / double(NinvptBins);

        // constant for phi cell prediction

        signedCurvaturePerInvPt = 3.e-4 * par.magneticField;

        for (int i = 0; i != NphiBins; ++i)
            for (int j = 0; j != NetaBins; ++j)
                for (int k = 0; k != NinvptBins; ++k) {
                    ArrElem[i][j][k].imeanPhi   = i + 0.5;
                    ArrElem[i][j][k].jmeanEta   = j + 0.5;
                    ArrElem[i][j][k].kmeanInvPt = k + 0.5;
                    ArrElem[i][j][k].i_index    = i;
                    ArrElem[i][j][k].j_index    = j;
                    ArrElem[i][j][k].k_index    = k;
                }

        /*		   
			for(int i = 0; i != NphiBins; ++i)
				for(int j = 0; j != NetaBins; ++j)
					for(int k = 0; k != NinvptBins; ++k)
						ArrElem[i][j][k].allHits = &allHits;	
		*/
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    void initHists() {

        HDeltaEta          = new TH1D("DeltaEta", "DeltaEta", 21, -10.5, +10.5);
        HDeltaEtaCorrected = new TH1D("DeltaEtaCorrected", "DeltaEtaCorrected", 21, -10.5, +10.5);
        HDeltaPhi          = new TH1D("DeltaPhi", "DeltaPhi", 61, -30.5, +30.5);
        HDeltaPhi2         = new TH1D("DeltaPhi2", "DeltaPhi2", 21, -10.5, +10.5);
        HInvptTimesR       = new TH1D("InvptTimesR", "InvptTimesR", 1000, 0., 0.3);
        HDeltaPhiVsInvptTimesR =
            new TH2I("DeltaPhiVsInvptTimesR", "DeltaPhiVsInvptTimesR", 100, 0., 0.3, 510, -10.5, +40.5);
        H3D_HTcellx      = new TH3I("H3D_HTcellx", "H3D_HTcellx; x1; x2; t", 40, 0., 1., 40, 0., 1., 40, 0., 1.);
        H3D_HTcellu      = new TH3I("H3D_HTcellu", "H3D_HTcellu; x1; u1; u2", 40, 0., 1., 40, 0., 1., 40, 0., 1.);
        H3D_HTAtrainStat = new TH3I("H3D_HTAtrainStat", "H3D_HTAtrainStat; i; j; k", NphiBins, -0.5, NphiBins - 0.5,
                                    NetaBins, -0.5, NetaBins - 0.5, NinvptBins, -0.5, NinvptBins - 0.5);
        HstatPhi         = new TH1D("HstatPhi", "HstatPhi", NphiBins, -0.5, NphiBins - 0.5);
        HstatEta         = new TH1D("HstatEta", "HstatEta", NetaBins, -0.5, NetaBins - 0.5);
        HstatInvpt       = new TH1D("HstatInvpt", "HstatInvpt", NinvptBins, -0.5, NinvptBins - 0.5);

        Hcellx1 = new TH1D("Hcellx1", "Hcellx1", 100, 0., 1.);
        Hcellx2 = new TH1D("Hcellx2", "Hcellx2", 100, 0., 1.);
        Hcellt  = new TH1D("Hcellt", "Hcellt", 100, 0., 1.);
        Hcellu0 = new TH1D("Hcellu0", "Hcellu0", 100, 0., 1.);
        Hcellu1 = new TH1D("Hcellu1", "Hcellu1", 100, 0., 1.);
        Hcellu2 = new TH1D("Hcellu2", "Hcellu2", 100, 0., 1.);
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    void print(std::ostream& out, int level = 0) {

        long int totalNbins = NphiBins * NetaBins * NinvptBins;
        out << "Total number of bins in the array: " << totalNbins << std::endl;

        out << "phi bins: " << NphiBins << "; phi min: " << phiMin << "; phi max: " << phiMin + NphiBins * phiStep
            << "; phi step: " << phiStep << std::endl;
        out << "eta bins: " << NetaBins << "; eta min: " << etaMin << "; eta max: " << etaMin + NetaBins * etaStep
            << "; eta step: " << etaStep << std::endl;
        out << "invPT bins: " << NinvptBins << "; invPt min: " << invPtMin
            << "; invPt max: " << invPtMin + NinvptBins * invPtStep << "; invPt step: " << invPtStep << std::endl;

        //out << "ETA " << NetaBins << " " << etaMin << " " << etaStep << std::endl;
        //out << "INVPT " << NinvptBins << " " << invPtMin << " " << invPtStep << std::endl;

        if (level > 0)
            for (int i = 0; i != NphiBins; ++i)
                for (int j = 0; j != NetaBins; ++j)
                    for (int k = 0; k != NinvptBins; ++k) {

                        out << std::endl
                            << "cell ijk: " << i << " " << j << " " << k << "      "
                            << " phi: [" << phiMin + i * phiStep << "," << phiMin + (i + 1) * phiStep << "]"
                            << " eta: [" << etaMin + j * etaStep << "," << etaMin + (j + 1) * etaStep << "]"
                            << " invPt: [" << invPtMin + k * invPtStep << "," << invPtMin + (k + 1) * invPtStep << "]"
                            << std::endl;

                        ArrElem[i][j][k].print(out);
                    }
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    void write(std::ostream& out) {

        out << TrainingPhase << std::endl;

        // parameters to check consistency of file with current HTM

        out << NphiBins << " " << phiMin << " " << phiStep << std::endl;
        out << NetaBins << " " << etaMin << " " << etaStep << std::endl;
        out << NinvptBins << " " << invPtMin << " " << invPtStep << std::endl;

        for (int i = 0; i != NphiBins; ++i)
            for (int j = 0; j != NetaBins; ++j)
                for (int k = 0; k != NinvptBins; ++k) {
                    for (auto it = ArrElem[i][j][k].layerIndHitStat.begin();
                         it != ArrElem[i][j][k].layerIndHitStat.end(); ++it) {
                        out << i << " " << j << " " << k << "  " << it->first << " "; // cell coordinates and layer
                        out << ArrElem[i][j][k].minLayers << " " << ArrElem[i][j][k].maxLayers << " ";
                        (it->second).write(out); // umin-umax
                        out << std::endl;
                    }
                }
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    void clear() {

        for (int i = 0; i != NphiBins; ++i)
            for (int j = 0; j != NetaBins; ++j)
                for (int k = 0; k != NinvptBins; ++k)
                    ArrElem[i][j][k].layerIndHitStat.clear();
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    int read(std::istream& in) {

        int      i, j, k, iLay;
        unsigned minLayers, maxLayers;
        double   angle, mean[3], low[3], high[3], u0L, u0H, u1L, u1H, u2L, u2H;
        unsigned _NphiBins, _NetaBins, _NinvptBins;
        double   _phiMin, _phiStep, _etaMin, _etaStep, _invPtMin, _invPtStep;
        int      retcode;

        in >> retcode;

        // check consistency of file with current HTM

        double toll = 1.e-5; // relative tolerance

        in >> _NphiBins >> _phiMin >> _phiStep;
        in >> _NetaBins >> _etaMin >> _etaStep;
        in >> _NinvptBins >> _invPtMin >> _invPtStep;

        if ((_NphiBins - NphiBins) != 0)
            return -1;
        if (fabs((_phiMin - phiMin) / (_phiMin + phiMin)) > toll)
            return -2;
        if (fabs((_phiStep - phiStep) / (_phiStep + phiStep)) > toll)
            return -3;

        if ((_NetaBins - NetaBins) != 0)
            return -4;
        if (fabs((_etaMin - etaMin) / (_etaMin + etaMin)) > toll)
            return -5;
        if (fabs((_etaStep - etaStep) / (_etaStep + etaStep)) > toll)
            return -6;

        if ((_NinvptBins - NinvptBins) != 0)
            return -7;
        if (fabs((_invPtMin - invPtMin) / (_invPtMin + invPtMin)) > toll)
            return -8;
        if (fabs((_invPtStep - invPtStep) / (_invPtStep + invPtStep)) > toll)
            return -9;

        // read file and fill HT array

        for (; in;) {
            in >> i >> j >> k >> iLay >> minLayers >> maxLayers >> mean[0] >> mean[1] >> mean[2] >> low[0] >> low[1] >>
                low[2] >> high[0] >> high[1] >> high[2] >> angle >> u0L >> u1L >> u2L >> u0H >> u1H >> u2H;
            //std::cout << "read "<< i << " " << j << " " << k << " " << iLay << " "
            //	<< mean[0] << " " << mean[1] << " " << mean[2] << " " << angle << std::endl;
            if (i < 0 || i >= (int)NphiBins)
                return -1;
            if (j < 0 || j >= (int)NetaBins)
                return -1;
            if (k < 0 || k >= (int)NinvptBins)
                return -1;
            Hitstat hs;
            hs.mean[0]  = mean[0];
            hs.mean[1]  = mean[1];
            hs.mean[2]  = mean[2];
            hs.low[0]   = low[0];
            hs.low[1]   = low[1];
            hs.low[2]   = low[2];
            hs.high[0]  = high[0];
            hs.high[1]  = high[1];
            hs.high[2]  = high[2];
            hs.rotAngle = angle;
            hs.cosAngle = cos(angle);
            hs.sinAngle = sin(angle);
            hs.u0Min    = u0L * 1.01;
            hs.u0Max    = u0H * 1.01;
            hs.u1Min    = u1L * 1.01;
            hs.u1Max    = u1H * 1.01;
            hs.u2Min    = u2L * 1.01;
            hs.u2Max    = u2H * 1.01;

            ArrElem[i][j][k].layerIndHitStat[iLay] = hs;
            ArrElem[i][j][k].minLayers             = minLayers;
            ArrElem[i][j][k].maxLayers             = maxLayers;
        }
        if (retcode < 3)
            ++retcode;
        return retcode;

    } // end read

    ////////////////////////////////////////////////////////////////////////////////////////

    // diagonalize each element of the array

    void diagonalize() {

        for (int i = 0; i != NphiBins; ++i)
            for (int j = 0; j != NetaBins; ++j)
                for (int k = 0; k != NinvptBins; ++k)
                    ArrElem[i][j][k].diagonalize();
    }

    ////////////////////////////////////////////////////////////////////////////////////////
    // find cell i,j,k in HT array from track parameters

    int getCell(Track& t, int& i, int& j, int& k) {

        i = (t.phi - phiMin) / phiStep;
        j = (t.eta - etaMin) / etaStep;
        k = (t.invPt - invPtMin) / invPtStep;

        // check boundaries - return -1 if track parameters are out of bounds

        if (i < 0 || i >= (int)NphiBins)
            return -1;
        if (j < 0 || j >= (int)NetaBins)
            return -1;
        if (k < 0 || k >= (int)NinvptBins)
            return -1;

        return 0;
    }

    ////////////////////////////////////////////////////////////////////////////////////////
    int getCell_i(double xPhi, int& i) {
        i = (xPhi - phiMin) / phiStep;
        if (i < 0 || i >= (int)NphiBins)
            return -1;
        else
            return 0;
    }

    ////////////////////////////////////////////////////////////////////////////////////////
    int getCell_j(double xEta, int& j) {
        j = (xEta - etaMin) / etaStep;
        if (j < 0 || j >= (int)NetaBins)
            return -1;
        else
            return 0;
    }

    ////////////////////////////////////////////////////////////////////////////////////////
    int getCell_k(double xInvPt, int& k) {
        k = (xInvPt - invPtMin) / invPtStep;
        if (k < 0 || k >= (int)NinvptBins)
            return -1;
        else
            return 0;
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    // Train this Array with hit h coming from track t

    int train(Hit& h, Track t, unsigned i, unsigned j, unsigned k) {

        // train the appropriate cell

        ArrElem[i][j][k].train(h);

        // fill max and min for number of layers traversed by this track

        unsigned nL = t.hitList.size();
        if (nL < ArrElem[i][j][k].minLayers)
            ArrElem[i][j][k].minLayers = nL;
        if (nL > ArrElem[i][j][k].maxLayers)
            ArrElem[i][j][k].maxLayers = nL;

        // Fill the 3D plot just for the one special test cell

        if (i == plotBinX && j == plotBinY && k == plotBinZ) {
            int lay = h.CellID;
            if (lay == (int)plotLay) {
                // This is the special cell we want to plot

                if (H3D_HTcellx->GetEntries() == 0) { // init hist limits (first fill only)
                    // Set histogram boundaries
                    double xLow  = ArrElem[i][j][k].layerIndHitStat[lay].low[0];
                    double xHigh = ArrElem[i][j][k].layerIndHitStat[lay].high[0];
                    double yLow  = ArrElem[i][j][k].layerIndHitStat[lay].low[1];
                    double yHigh = ArrElem[i][j][k].layerIndHitStat[lay].high[1];
                    double zLow  = ArrElem[i][j][k].layerIndHitStat[lay].low[2];
                    double zHigh = ArrElem[i][j][k].layerIndHitStat[lay].high[2];
                    H3D_HTcellx->GetXaxis()->SetLimits(xLow, xHigh);
                    H3D_HTcellx->GetYaxis()->SetLimits(yLow, yHigh);
                    H3D_HTcellx->GetZaxis()->SetLimits(zLow, zHigh);

                    Hcellx1->GetXaxis()->SetLimits(xLow, xHigh);
                    Hcellx2->GetXaxis()->SetLimits(yLow, yHigh);
                    Hcellt->GetXaxis()->SetLimits(zLow, zHigh);

                    xLow  = ArrElem[i][j][k].layerIndHitStat[lay].u0Min;
                    xHigh = ArrElem[i][j][k].layerIndHitStat[lay].u0Max;
                    yLow  = ArrElem[i][j][k].layerIndHitStat[lay].u1Min;
                    yHigh = ArrElem[i][j][k].layerIndHitStat[lay].u1Max;
                    zLow  = ArrElem[i][j][k].layerIndHitStat[lay].u2Min;
                    zHigh = ArrElem[i][j][k].layerIndHitStat[lay].u2Max;
                    H3D_HTcellu->GetXaxis()->SetLimits(xLow, xHigh);
                    H3D_HTcellu->GetYaxis()->SetLimits(yLow, yHigh);
                    H3D_HTcellu->GetZaxis()->SetLimits(zLow, zHigh);

                    Hcellu0->GetXaxis()->SetLimits(xLow, xHigh);
                    Hcellu1->GetXaxis()->SetLimits(yLow, yHigh);
                    Hcellu2->GetXaxis()->SetLimits(zLow, zHigh);

                } // end init hist limits

                H3D_HTcellx->Fill(h.x1, h.x2, h.t);
                H3D_HTcellu->Fill(h.u0, h.u1, h.u2);
                Hcellx1->Fill(h.x1);
                Hcellx2->Fill(h.x2);
                Hcellt->Fill(h.t);
                Hcellu0->Fill(h.u0);
                Hcellu1->Fill(h.u1);
                Hcellu2->Fill(h.u2);
            }
        }

        return 0;

    } // end train

    ////////////////////////////////////////////////////////////////////////////////////////

    void reset() {

        for (unsigned iPhi = 0; iPhi != NphiBins; ++iPhi)
            for (unsigned iEta = 0; iEta != NetaBins; ++iEta)
                for (unsigned iInvpt = 0; iInvpt != NinvptBins; ++iInvpt)
                    ArrElem[iPhi][iEta][iInvpt].reset();

        allHits.clear();

        return;
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    int fill(Hit& h) {

        // add this hit to the pool of hits for this event;

        h.ID = allHits.size();
        allHits.push_back(h);

        int iPhi1   = 0;
        int iPhi2   = NphiBins;
        int iEta1   = 0;
        int iEta2   = NetaBins;
        int iInvpt1 = 0;
        int iInvpt2 = NinvptBins;

        for (int iInvpt = iInvpt1; iInvpt != iInvpt2; ++iInvpt) {

            // indexes of the hit coordinates

            int iPhiMap  = 0;
            int iPhiMap2 = 0;

            int iEtaMap     = 0;
            int iEtaMapCorr = 0;
            for (int iEta = iEta1; iEta != iEta2; ++iEta) {

                for (int iPhi = iPhi1; iPhi != iPhi2; ++iPhi) {

                    int strike = ArrElem[iPhi][iEta][iInvpt].fill(h);

                    if (strike == 0) {

                        // insert here code to build candidate list on the fly ////////////

                        int deltaEta     = iEtaMap - iEta;
                        int deltaEtaCorr = iEtaMapCorr - iEta;

                        HDeltaEta->Fill(deltaEta);
                        HDeltaEtaCorrected->Fill(deltaEtaCorr);

                        int deltaPhi = iPhiMap - iPhi;
                        HDeltaPhi->Fill(deltaPhi);

                        int deltaPhi2 = iPhiMap2 - iPhi;
                        HDeltaPhi2->Fill(deltaPhi2);

                    } // end if strike == 0

                } // end loop on phi
            }     // end loop on eta
        }         // end loop on pt

        return 0;

    } // end fill

    ////////////////////////////////////////////////////////////////////////////////////////

    unsigned getBestCell(unsigned& iPhi_, unsigned& iEta_, unsigned& iInvpt_) {

        unsigned maxHits = 0;
        for (unsigned iPhi = 0; iPhi != NphiBins; ++iPhi)
            for (unsigned iEta = 0; iEta != NetaBins; ++iEta)
                for (unsigned iInvpt = 0; iInvpt != NinvptBins; ++iInvpt) {
                    unsigned nHits = ArrElem[iPhi][iEta][iInvpt].nHitLayers;

                    if (nHits > maxHits) {
                        maxHits = nHits;
                        iPhi_   = iPhi;
                        iEta_   = iEta;
                        iInvpt_ = iInvpt;
                    }
                }

        return maxHits;

    } // end getBestCell

    ////////////////////////////////////////////////////////////////////////////////////////

    struct Pars {
        unsigned iPhi;
        unsigned iEta;
        unsigned iInvpt;
        unsigned nLayers;
    };
    std::vector<Pars> cellCandidateList;

    unsigned getCellCandidates() {

        cellCandidateList.clear();
        for (unsigned iPhi = 0; iPhi != NphiBins; ++iPhi)
            for (unsigned iEta = 0; iEta != NetaBins; ++iEta)
                for (unsigned iInvpt = 0; iInvpt != NinvptBins; ++iInvpt) {
                    unsigned nLayers = ArrElem[iPhi][iEta][iInvpt].nHitLayers;
                    //if(nLayers >= (ArrElem[iPhi][iEta][iInvpt].minLayers))
                    //if(nLayers >= (ArrElem[iPhi][iEta][iInvpt].minLayers-1))
                    if (nLayers >= par.reco_minLayersForFit) {
                        Pars p;
                        p.iPhi    = iPhi;
                        p.iEta    = iEta;
                        p.iInvpt  = iInvpt;
                        p.nLayers = nLayers;
                        cellCandidateList.push_back(p);
                    }
                }
        return (unsigned)cellCandidateList.size();

    } // end getCellCandidates

    ////////////////////////////////////////////////////////////////////////////////////////

    void printCellCandidateList(std::ostream& out) {

        for (unsigned i = 0; i != cellCandidateList.size(); ++i) {
            unsigned iPhi    = cellCandidateList[i].iPhi;
            unsigned iEta    = cellCandidateList[i].iEta;
            unsigned iInvpt  = cellCandidateList[i].iInvpt;
            unsigned nLayers = cellCandidateList[i].nLayers;
            out << "Candidate " << i << ": " << iPhi << " " << iEta << " " << iInvpt << " nLayers: " << nLayers
                << std::endl;
            HTArrayElement elem = ArrElem[iPhi][iEta][iInvpt];
            elem.printHits(out);
        }

    } // end printCellCandidateList

    ////////////////////////////////////////////////////////////////////////////////////////

    void writeCellCandidateList(std::ostream& out) {
        out << cellCandidateList.size() << std::endl;
        for (unsigned i = 0; i != cellCandidateList.size(); ++i) {
            unsigned iPhi    = cellCandidateList[i].iPhi;
            unsigned iEta    = cellCandidateList[i].iEta;
            unsigned iInvpt  = cellCandidateList[i].iInvpt;
            unsigned nLayers = cellCandidateList[i].nLayers;
            out << iPhi << " " << iEta << " " << iInvpt << " " << nLayers << std::endl;
            HTArrayElement elem = ArrElem[iPhi][iEta][iInvpt];
            elem.printHits(out);
        }

    } // end writeCellCandidateList

    ////////////////////////////////////////////////////////////////////////////////////////

    void trainStat(int i, int j, int k) {

        H3D_HTAtrainStat->Fill(i, j, k); // record track stat per cell
        HstatPhi->Fill(i);
        HstatEta->Fill(j);
        HstatInvpt->Fill(k);
    }

    //////////////////////////////////////////////////////////////////////////////////////

    std::string makeFileName(std::string prefix, int i, int j, int k) {

        std::stringstream fileName;

        fileName << prefix << "I" << i << "J" << j << "K" << k << ".txt";

        return fileName.str();
    }

    //////////////////////////////////////////////////////////////////////////////////////

}; // end HTArray

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

#endif // MDHT_HTARRAY_CPP
