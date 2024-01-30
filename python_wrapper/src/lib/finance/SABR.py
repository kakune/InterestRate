import numpy as np


def approxImpVol(inVolvol, inExponent, inCorr, inInitPrice, inInitVol, inStrike, inTimeMaturity):
    lFK = (inInitPrice * inStrike) ** ((1.0 - inExponent) * 0.5)
    lLogFK = np.log(inInitPrice / inStrike)
    lZeta = inVolvol * lFK * lLogFK / inInitVol
    lX = np.log((np.sqrt(1.0 - 2.0 * inCorr * lZeta + lZeta **
                2) - inCorr + lZeta) / (1.0 - inCorr))
    lBetaLog = ((1.0 - inExponent) * lLogFK)**2

    lResult = (inInitVol / lFK)
    if (lX != 0.0):
        lResult *= (lZeta / lX)
    lResult /= 1.0 + lBetaLog / 24.0 + lBetaLog ** 2 / 1920.0

    lFactor = 1.0
    lFactor += inTimeMaturity * \
        ((1.0 - inExponent) * inInitVol / lFK) ** 2 / 24.0
    lFactor += 0.25 * inTimeMaturity * inCorr * \
        inInitVol * inVolvol * inExponent / lFK
    lFactor += inTimeMaturity * \
        (2.0 - 3.0 * inCorr * inCorr) * inVolvol * inVolvol / 24.0
    return lResult * lFactor
