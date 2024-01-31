import subprocess


def buildCmakeRelease(
    inPathSource: str
):
    try:
        lMakeCommand = ["cmake", "-S",
                        inPathSource, "-B", inPathSource + "/build", "-DCMAKE_BUILD_TYPE=Release"]
        lBuildCommand = ["cmake", "--build", inPathSource + "/build"]
        subprocess.run(lMakeCommand, check=True,
                       capture_output=True, text=True)
        subprocess.run(lBuildCommand, check=True,
                       capture_output=True, text=True)
        # print(process.stdout)
    except subprocess.CalledProcessError as e:
        print(f"returncode:{e.returncode}")
        print(e.stderr)
        exit()


def buildCmakeDebug(
    inPathSource: str
):
    try:
        lMakeCommand = ["cmake", "-S",
                        inPathSource, "-B", inPathSource + "/build", "-DCMAKE_BUILD_TYPE=Debug"]
        lBuildCommand = ["cmake", "--build", inPathSource + "/build"]
        subprocess.run(lMakeCommand, check=True,
                       capture_output=True, text=True)
        subprocess.run(lBuildCommand, check=True,
                       capture_output=True, text=True)
        # print(process.stdout)
    except subprocess.CalledProcessError as e:
        print(f"returncode:{e.returncode}")
        print(e.stderr)
        exit()


def runExe(
    inPathExe: str,
    args: tuple[str] = ()
):
    try:
        process = subprocess.run(
            [inPathExe, *args], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"returncode:{e.returncode}")
        print(e.stderr)
        exit()
    return process
