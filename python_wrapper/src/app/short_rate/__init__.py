import os
import sys
sys.path.append(os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', '..')))
gPathCurrent = os.path.abspath(__file__)
gPathProject = os.path.split(os.path.split(os.path.split(os.path.split(
    os.path.split(gPathCurrent)[0])[0])[0])[0])[0]
gPathCppDir = os.path.join(gPathProject, "cpp_calculator")
gPathCppExeDir = os.path.join(gPathCppDir, "build", "src", "app", "short_rate")
gPathParam = os.path.join(gPathProject, "parameters", "short_rate.ini")
gPathMarketDir = os.path.join(gPathProject, "parameters", "market")
gPathOutputDir = os.path.join(gPathProject, "output", "short_rate")

os.makedirs(gPathOutputDir, exist_ok=True)
