import __init__
import os
import pandas as pd
import matplotlib.pyplot as plt
from lib.cpp.cmake import buildCmakeRelease, runExe
from lib.plot.graph import plotGraph
from lib.finance.SABR import approxImpVol
from lib.utils.parameters import Parameters

gPathCppExe = os.path.join(
    __init__.gPathCppExeDir, "SABR_implied_volatility")
gPathParam = os.path.join(__init__.gPathParamDir, "SABR.ini")

gNameSections = ["PARAM1", "PARAM2", "PARAM3"]
gNameOutputs = ["SABR_output_" + gNameSection +
                ".csv" for gNameSection in gNameSections]
gPathOutputs = [os.path.join(__init__.gPathOutputDir, gNameOutput)
                for gNameOutput in gNameOutputs]
gPathGraph = os.path.join(__init__.gPathOutputDir, "SABR_graph_multi.png")

fig = None
ax = None


if __name__ == '__main__':
    buildCmakeRelease(__init__.gPathCppDir)
    for gNameSection, gPathOutput in zip(gNameSections, gPathOutputs):
        runExe(
            gPathCppExe,
            (gPathParam, gNameSection, gPathOutput)
        )

        lParam = Parameters()
        lParam.readParameters(gPathParam)
        lParam.setNameCurrentSection(gNameSection)

        lDataFrame = pd.read_csv(gPathOutput)
        fig, ax = plotGraph(
            xs=lDataFrame["Strike"],
            ys=lDataFrame["ImpVol"] * 100,
            label=r"N: $S_{0} = $"+str(lParam("InitPrice")),
            title=r"$\beta = $"+str(lParam("Exponent")),
            xlabel=r"Strike $K$",
            ylabel=r"ImpVol $\sigma ~ (\mathrm{\%})$",
            fig=fig,
            ax=ax
        )

        lApprox = [
            approxImpVol(
                inStrike=lStrike,
                inInitPrice=lParam("InitPrice"),
                inInitVol=lParam("InitVol"),
                inCorr=lParam("Corr"),
                inExponent=lParam("Exponent"),
                inVolvol=lParam("Volvol"),
                inTimeMaturity=lParam("TimeMaturity")
            ) * 100
            for lStrike in lDataFrame["Strike"]
        ]

        fig, ax = plotGraph(
            xs=lDataFrame["Strike"],
            ys=lApprox,
            label=r"A: $S_{0} = $"+str(lParam("InitPrice")),
            fig=fig,
            ax=ax
        )

    plt.savefig(gPathGraph, bbox_inches='tight')
