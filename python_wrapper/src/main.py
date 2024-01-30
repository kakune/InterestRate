import os
import pandas as pd
import matplotlib.pyplot as plt
from lib.cpp.cmake import buildCmake, runExe
from lib.plot.graph import plotGraph
from lib.finance.SABR import approxImpVol
from lib.utils.parameters import Parameters

gPathCurrent = os.path.abspath(__file__)
gPathProject = os.path.split(os.path.split(
    os.path.split(gPathCurrent)[0])[0])[0]
gPathCppDir = os.path.join(gPathProject, "cpp_calculator")
gPathCppExe = os.path.join(gPathCppDir, "build", "src", "main")

gNameParam = "parameters.ini"
gPathParam = os.path.join(gPathProject, "parameters", gNameParam)

gNameOutput = "output.csv"
gPathOutput = os.path.join(gPathProject, "output", gNameOutput)
gNameGraph = "graph.png"
gPathGraph = os.path.join(gPathProject, "output", gNameGraph)

gNameSection = "PARAM1"


if __name__ == '__main__':
    buildCmake(gPathCppDir)
    runExe(gPathCppExe, (gPathParam, gNameSection, gPathOutput))
    lDataFrame = pd.read_csv(gPathOutput)
    fig, ax = plotGraph(
        xs=lDataFrame["Strike"],
        ys=lDataFrame["ImpVol"],
        label="Numerical"
    )
    lParam = Parameters()
    lParam.readParameters(gPathParam)
    lParam.setNameCurrentSection(gNameSection)
    lApprox = [approxImpVol(
        inVolvol=lParam("Volvol"),
        inExponent=lParam("Exponent"),
        inCorr=lParam("Corr"),
        inInitPrice=lParam("InitPrice"),
        inStrike=lStrike,
        inInitVol=lParam("InitVol"),
        inTimeMaturity=lParam("TimeMaturity")
    ) for lStrike in lDataFrame["Strike"]]

    fig, ax = plotGraph(
        xs=lDataFrame["Strike"],
        ys=lApprox,
        label="Approx",
        fig=fig,
        ax=ax
    )

    plt.savefig(gPathGraph)
