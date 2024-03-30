import subprocess
import os


def buildCmakeRelease(
    inPathSource: str
):
    try:
        os.chdir(inPathSource)
        subprocess.run(["make", "build"], check=True)

    except subprocess.CalledProcessError as e:
        print(f"returncode:{e.returncode}")
        print(e.stderr)
        exit()


def buildCmakeReleaseCUDA(
    inPathSource: str
):
    try:
        os.chdir(inPathSource)
        subprocess.run(["make", "build", "CUDA_ENABLED=1"], check=True)

    except subprocess.CalledProcessError as e:
        print(f"returncode:{e.returncode}")
        print(e.stderr)
        exit()


def buildCmakeDebug(
    inPathSource: str
):
    try:
        os.chdir(inPathSource)
        subprocess.run(["make", "bebug"], check=True)

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
