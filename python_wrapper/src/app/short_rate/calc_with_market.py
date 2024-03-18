import __init__
import os
import pandas as pd
import matplotlib.pyplot as plt
from lib.cpp.cmake import buildCmakeRelease, runExe
from lib.plot.graph import scatter3D
from lib.utils.parameters import Parameters

gNameModel = "HoLee"
gNameSection = "PARAM1"
gNameSectionMarket = "JGB"
gNameCSVMarket = "JGB"

gPathCppExe = os.path.join(
    __init__.gPathCppExeDir, "calc_with_market")
gPathOutput = os.path.join(__init__.gPathOutputDir,
                           gNameModel + "_with_market_output.csv")
gPathParam = os.path.join(__init__.gPathParamDir, gNameModel + ".ini")
gPathMarketParam = os.path.join(__init__.gPathMarketDir, "parameters.ini")
gPathCSVMarket = os.path.join(__init__.gPathMarketDir, gNameCSVMarket + ".csv")
gPathGraphZCB = os.path.join(__init__.gPathOutputDir,
                             gNameModel + "_with_market_ZCB_graph.png")
gPathGraphFR = os.path.join(__init__.gPathOutputDir,
                            gNameModel + "_with_market_FR_graph.png")

if __name__ == '__main__':
    buildCmakeRelease(__init__.gPathCppDir)
    runExe(
        gPathCppExe,
        (gPathOutput, gNameModel, gNameSection, gPathParam,
         gNameSectionMarket, gPathMarketParam, gPathCSVMarket)
    )

    lDataFrame = pd.read_csv(gPathOutput)
    fig, ax = scatter3D(
        xs=lDataFrame["Maturity"],
        ys=lDataFrame["Start"],
        zs=lDataFrame["PriceZCB"],
        xlabel='Maturity',
        ylabel='Start'
    )
    plt.savefig(gPathGraphZCB, bbox_inches='tight')

    fig, ax = scatter3D(
        xs=lDataFrame["Maturity"],
        ys=lDataFrame["Start"],
        zs=lDataFrame["ForwardRate"],
        xlabel='Maturity',
        ylabel='Start'
    )
    plt.savefig(gPathGraphFR, bbox_inches='tight')
