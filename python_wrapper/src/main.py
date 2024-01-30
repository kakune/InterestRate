import os
from lib.cpp.cpp_module import buildCmake, runExe

gPathCurrent = os.getcwd()
gPathProject = os.path.split(gPathCurrent)[0]
gPathCppDir = os.path.join(gPathProject, "cpp_calculator")
gPathCppExe = os.path.join(gPathCppDir, "build", "src", "main")

gNameParam = "parameters.ini"
gPathParam = os.path.join(gPathProject, "parameters", gNameParam)

gNameOutput = "output.csv"
gPathOutput = os.path.join(gPathProject, "output", gNameOutput)


if __name__ == '__main__':
    buildCmake(gPathCppDir)
    runExe(gPathCppExe, (gPathParam, "PARAM1", gPathOutput))
