
#ifndef MDHT_BIBFILEREADER_CPP
#define MDHT_BIBFILEREADER_CPP

//////////////////////////////////////////////////////
//
//  BibFileReader.cpp
//  HT_TrackFinding
//
//  Created by Luciano Ristori on 1/31/23
//
//////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include "Hit.h"

class BibFileReader {

  public:
    //-------------------------------------------------
    int readFile(std::string BibFileName) {
        std::cout << "BibFileReader  file " << BibFileName << std::endl;
        return 0;
    }

    //-------------------------------------------------
    unsigned size() { return 0; }

    //-------------------------------------------------
    Hit randomBibHit() {
        Hit h('G');
        return h;
    }

}; // end class BibFileReader

#endif // MDHT_BIBFILEREADER_CPP
