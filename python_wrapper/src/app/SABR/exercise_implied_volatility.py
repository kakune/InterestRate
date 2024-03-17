import __init__
import os
import pandas as pd
import matplotlib.pyplot as plt
from lib.cpp.cmake import buildCmakeRelease, runExe
from lib.plot.graph import plotGraph
from lib.finance.exercise import approxImpVol
from lib.utils.parameters import Parameters

gPathCppExe = os.path.join(
    __init__.gPathCppExeDir, "exercise_implied_volatility")
gPathParam = os.path.join(__init__.gPathParamDir, "exercise.ini")
gPathOutput = os.path.join(__init__.gPathOutputDir, "exercise_output.csv")
gPathGraph = os.path.join(__init__.gPathOutputDir, "exercise_graph.png")

gNameSection = "PARAM1"


if __name__ == '__main__':
    buildCmakeRelease(__init__.gPathCppDir)
    runExe(
        gPathCppExe,
        (gPathParam, gNameSection, gPathOutput)
    )

    lDataFrame = pd.read_csv(gPathOutput)
    fig, ax = plotGraph(
        xs=lDataFrame["Strike"],
        ys=lDataFrame["ImpVol"],
        label="Numerical"
    )

    lParam = Parameters()
    lParam.readParameters(gPathParam)
    lParam.setNameCurrentSection(gNameSection)

    lApprox = [
        approxImpVol(
            inStrike=lStrike,
            inInitPrice=lParam("InitPrice"),
            inInitVol=lParam("InitVol"),
            inCorr=lParam("Corr"),
            inExponent=lParam("Exponent"),
            inShift=lParam("Shift"),
            inStandard=lParam("Standard"),
            inVolvol=lParam("Volvol"),
            inTimeMaturity=lParam("TimeMaturity")
        )
        for lStrike in lDataFrame["Strike"]
    ]

    fig, ax = plotGraph(
        xs=lDataFrame["Strike"],
        ys=lApprox,
        label="Approx",
        fig=fig,
        ax=ax
    )

    plt.savefig(gPathGraph)
