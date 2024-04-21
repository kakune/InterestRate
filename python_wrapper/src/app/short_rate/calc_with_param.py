import __init__
import os
import pandas as pd
import matplotlib.pyplot as plt
from lib.cpp.cmake import buildCmakeRelease, runExe
from lib.plot.graph import scatter3D
from lib.utils.parameters import Parameters

gNameSection = "ConstantAffine"

gPathCppExe = os.path.join(
    __init__.gPathCppExeDir, "SR_param")
gPathOutput = os.path.join(__init__.gPathOutputDir,
                           gNameSection + "_SR_output.csv")
gPathGraphZCB = os.path.join(__init__.gPathOutputDir,
                             gNameSection + "_SR_ZCB_graph.png")
gPathGraphFR = os.path.join(__init__.gPathOutputDir,
                            gNameSection + "_SR_FR_graph.png")

if __name__ == '__main__':
    buildCmakeRelease(__init__.gPathCppDir)
    runExe(
        gPathCppExe,
        (gPathOutput, gNameSection, __init__.gPathParam)
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
