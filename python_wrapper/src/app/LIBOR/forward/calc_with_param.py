import __init__
import os
import pandas as pd
import matplotlib.pyplot as plt
from lib.cpp.cmake import buildCmakeRelease, runExe
from lib.plot.graph import scatter3D
from lib.utils.parameters import Parameters

gNameSection = "SABR"

gPathCppExe = os.path.join(
    __init__.gPathCppExeDir, "FR_param")
gPathOutput = os.path.join(__init__.gPathOutputDir,
                           gNameSection + "_FR_output.csv")
gPathGraphImpVol = os.path.join(__init__.gPathOutputDir,
                            gNameSection + "_FR_ImpVol.png")

if __name__ == '__main__':
    buildCmakeRelease(__init__.gPathCppDir)
    runExe(
        gPathCppExe,
        (gPathOutput, gNameSection, __init__.gPathParam)
    )

    lDataFrame = pd.read_csv(gPathOutput)

    fig, ax = scatter3D(
        ys=lDataFrame["Start"],
        xs=lDataFrame["Strike"],
        zs=lDataFrame["ImpVol"],
        ylabel='Start',
        xlabel='Strike'
    )
    plt.savefig(gPathGraphImpVol, bbox_inches='tight')
