import numpy as np


def approxImpVol(
    inStrike: float,
    inInitPrice: float,
    inInitVol: float,
    inCorr: float,
    inExponent: float,
    inShift: float,
    inStandard: float,
    inVolvol: float,
    inTimeMaturity: float
):
    lShiftExp = inShift ** inExponent

    def funcC(f):
        lTmp = f ** inExponent
        return inStandard * lTmp / (lTmp + lShiftExp)

    lFav = np.sqrt(inInitPrice * inStrike)
    lFavExp = lFav ** inExponent
    lFK = lFav ** (1.0 - inExponent)

    lLogFK = np.log(inInitPrice / inStrike)
    lZeta = (inVolvol / inInitVol) * (inInitPrice - inStrike) / funcC(lFav)

    lGamma = inExponent * lShiftExp / (lFav * (lFavExp + lShiftExp))

    lPhi = inStandard * inStandard
    lPhi *= inExponent * lShiftExp / (lFK * lFK)
    lPhi *= (inExponent - 2.0) * lShiftExp - 2.0*(inExponent + 1) * lFavExp
    lPhi /= (lShiftExp + lFavExp) ** 4.0

    lSigma = inInitVol * inStandard / (lFav + lShiftExp * lFK)

    lX = np.log((np.sqrt(1.0 - 2.0 * inCorr * lZeta + lZeta **
                2) - inCorr + lZeta) / (1.0 - inCorr))

    def funcFac(inExp):
        lResult = lFav**inExp
        lResult *= 1.0 + (inExp * lLogFK)**2 / 24.0 + \
            (inExp * lLogFK)**4 / 1920.0
        return lResult

    lResult = (inInitVol * inStandard)
    if (lX != 0.0):
        lResult *= (lZeta / lX)
    lResult /= funcFac(1.0) + lShiftExp * funcFac(1.0-inExponent)

    lFactor = 1.0
    lFactor += inTimeMaturity * \
        (inInitVol * inInitVol * (lPhi + lSigma * lSigma)) / 24.0
    lFactor += 0.25 * inTimeMaturity * inCorr * \
        inInitVol * inVolvol * lGamma * funcC(lFav)
    lFactor += inTimeMaturity * \
        (2.0 - 3.0 * inCorr * inCorr) * inVolvol * inVolvol / 24.0
    return lResult * lFactor
