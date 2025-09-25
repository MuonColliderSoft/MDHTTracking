
#ifndef MDHT_CELLIDTOLAYER_H
#define MDHT_CELLIDTOLAYER_H

//
//  CellIDtoLayer.h
//
//  Created by Luciano Ristori on 4/23/25
//
//

int CellIDtoLayer(int CellID) {

    const unsigned int system = (unsigned)(CellID & 0x1f);
    const unsigned int side   = (int)((CellID >> 5) & 0x3);
    const unsigned int layer  = (unsigned)((CellID >> 7) & 0x3f);

    return 100 * system + 10 * side + layer;
};

#endif //MDHT_CELLIDTOLAYER_H
