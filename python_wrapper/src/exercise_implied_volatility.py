import os
import pandas as pd
import matplotlib.pyplot as plt
from lib.cpp.cmake import buildCmakeRelease, runExe
from lib.plot.graph import plotGraph
from lib.finance.exercise import approxImpVol
from lib.utils.parameters import Parameters

gPathCurrent = os.path.abspath(__file__)
gPathProject = os.path.split(os.path.split(
    os.path.split(gPathCurrent)[0])[0])[0]
gPathCppDir = os.path.join(gPathProject, "cpp_calculator")
gPathCppExe = os.path.join(
    gPathCppDir, "build", "src", "exercise_implied_volatility")

gNameParam = "exercise.ini"
gPathParam = os.path.join(gPathProject, "parameters", gNameParam)

gNameOutput = "exercise_output.csv"
gPathOutput = os.path.join(gPathProject, "output", gNameOutput)
gNameGraph = "exercise_graph.png"
gPathGraph = os.path.join(gPathProject, "output", gNameGraph)

gNameSection = "PARAM1"


if __name__ == '__main__':
    buildCmakeRelease(gPathCppDir)
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
