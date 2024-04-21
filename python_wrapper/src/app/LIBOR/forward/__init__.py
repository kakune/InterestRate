import os
import sys
sys.path.append(os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')))
gPathCurrent = os.path.abspath(__file__)
gPathProject = os.path.split(os.path.split(os.path.split(os.path.split(os.path.split(
    os.path.split(gPathCurrent)[0])[0])[0])[0])[0])[0]
gPathCppDir = os.path.join(gPathProject, "cpp_calculator")
gPathCppExeDir = os.path.join(gPathCppDir, "build", "src", "app", "LIBOR", "forward")
gPathParam = os.path.join(gPathProject, "parameters", "LIBOR", "forward.ini")
gPathMarketDir = os.path.join(gPathProject, "parameters", "market")
gPathOutputDir = os.path.join(gPathProject, "output", "LIBOR", "forward")

os.makedirs(gPathOutputDir, exist_ok=True)
