#ifndef MDHT_HTAMAPPER_H
#define MDHT_HTAMAPPER_H

#include "Parameters.h"

extern Parameters par;

int trackToCell(int& i, int& j, int& k, Track track) {

    int NphiBins   = par.HTA_NphiBins;
    int NetaBins   = par.HTA_NetaBins;
    int NinvptBins = par.HTA_NinvptBins;

    i = (track.phi - par.HTA_t_phi + par.HTA_t_deltaPhi) / par.HTA_t_deltaPhi * (double)NphiBins / 2.;
    j = (track.eta - par.HTA_t_eta + par.HTA_t_deltaEta) / par.HTA_t_deltaEta * (double)NetaBins / 2.;
    k = (track.invPt - par.HTA_t_invPt_min) / (par.HTA_t_invPt_max - par.HTA_t_invPt_min) * (double)NinvptBins;

    if (i < 0)
        return -1;
    if (j < 0)
        return -2;
    if (k < 0)
        return -3;
    if (i >= NphiBins)
        return -4;
    if (j >= NetaBins)
        return -5;
    if (k >= NinvptBins)
        return -6;

    return 0;
}

#endif // MDHT_HTAMAPPER_H