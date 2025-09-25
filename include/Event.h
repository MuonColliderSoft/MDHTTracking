#ifndef MDHT_EVENT_CPP
#define MDHT_EVENT_CPP

//
//  Event.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 6/30/21
//  Revised for Muon Collider software integration 4/14/25
//

#include <iostream>
#include <vector>
#include <random>

#include "include/Track.h"
#include "include/Hit.h"

extern bool         verbose;
extern std::mt19937 generator; // Mersenne Twister

#include "include/BibFileReader.h"

class Event {

  public:
    std::vector<Track> trackList;
    std::vector<Hit>   hitList;

    //////////////////////////////////////////////////////////////////////////////////////

    // Constructor with background

    Event(TrackReader* newRecoTrack, int nTracks, BibFileReader& bibRead, long int backGnd = -1) {

        trackList.clear();
        hitList.clear();

        if (verbose)
            std::cout << "Event" << std::endl;

        // Create all tracks
        for (unsigned iT = 1; iT != (unsigned)nTracks + 1; ++iT) { // iT=0 is BIB

            Track thisTrack; //  create track

            //	if(iT == 0){ // put a dummy track at position 0;
            //		trackList.push_back(thisTrack);
            //		continue;
            //	}

            bool eof = newRecoTrack->read(thisTrack);
            if (eof)
                break;

            thisTrack.ID = iT; // assign track iD = track index

            // loop on hits of this track

            for (int iH = 0; iH != (int)thisTrack.hitList.size(); ++iH) {
                //thisTrack.hitList[iH].ID = hitList.size(); // set hit ID
                thisTrack.hitList[iH].trackInd = iT;      // set hit trackInd
                hitList.push_back(thisTrack.hitList[iH]); // add this hit to hitList

            } // end loop on hits of this track

            thisTrack.hitList.clear();      // clear hitList for this track
            trackList.push_back(thisTrack); // add this track to trackList

            if (verbose)
                std::cout << "Created track " << iT << std::endl;

        } // end loop on tracks

        if (backGnd > 0) {

            if (verbose)
                std::cout << "generate background" << std::endl;

            // Adds BIB hits to hitList for this event

            for (unsigned iH = 0; iH != backGnd; ++iH) {
                Hit bibHit = bibRead.randomBibHit();
                bibHit.ID  = hitList.size();
                hitList.push_back(bibHit);
            }
        }

        // shuffle hits

        std::shuffle(std::begin(hitList), std::end(hitList), generator);

        if (verbose)
            std::cout << "done shuffling hits" << std::endl;

        // define ID's of all hits and add hits to hitList of parent track

        for (unsigned iH = 0; iH != hitList.size(); ++iH) {
            Hit* thisHit      = &(hitList[iH]);
            thisHit->ID       = iH;
            unsigned trackInd = thisHit->trackInd;
            if (trackInd > 0)
                ((trackList[trackInd - 1]).hitList).push_back(*thisHit);
        }
        if (verbose)
            std::cout << "end of Event constructor" << std::endl;

    } // end event constructor

    //////////////////////////////////////////////////////////////////////////////////////

    void print(std::ostream& out, int mode = 0) {

        // mode = 0: print only track parameters
        // mode = 1: print all hits

        for (unsigned iT = 0; iT != trackList.size(); ++iT) {
            trackList[iT].print(out, mode);
        }

    } // end print

    //////////////////////////////////////////////////////////////////////////////////////

}; // end class Event

#endif //MDHT_EVENT_CPP
