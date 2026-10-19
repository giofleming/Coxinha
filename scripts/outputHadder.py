import glob
import argparse
from pathlib import Path
import subprocess


def getArgumentParser():

    parser = argparse.ArgumentParser(
        description="Command line arguments for full chain setup"
    )

    parser.add_argument(
        "--indir",
        "-i",
        dest="indir",
        help="Input dir where the outputs are",
        type=str,
        default="./output/",
    )

    parser.add_argument(
        "--regex",
        "-r",
        dest="regex",
        help="Regex to find outputs",
        type=str,
        default="FullChain",
    )

    parser.add_argument(
        "--outdir",
        "-o",
        dest="outdir",
        help="Output folder",
        type=str,
        default="mergedOutput",
    )

    return parser


def haddOutputs(indir=None,
                regex=None,
                outdir=None,
                iterdirs=["iter_0", "iter_1", "iter_2"]):

    inputFolders = glob.glob(indir+regex)
    print(inputFolders)

    # Create the hadded outputFolder
    Path(outdir).mkdir(parents=True, exist_ok=True)

    # List of files to be hadded
    iterfiles = ["*performance_*seed*.root",
                 "*performance_*ambi_tracks*.root"]

    # Merge the iterations

    for iterfile in iterfiles:

        for iterdir in iterdirs:

            print("Checking ... ", iterdir)
            filesToMerge = []
            # Create the iteration sub-folders
            Path(outdir+"/"+iterdir).mkdir(parents=True, exist_ok=True)

            # Get the inputs
            for inFolder in inputFolders:
                inFileName = glob.glob(inFolder+"/"+iterdir+"/"+iterfile)[0]
                filesToMerge.append(inFileName)

            outFileName = outdir+"/"+iterdir+"/"+inFileName.split("/")[-1]
            print("Merging", filesToMerge, "-->", outFileName)

            subprocess.run(["hadd", "-f", outFileName]+filesToMerge, check=True)

    print("---iterations merging done---")

    # Merge the merged files
    filesToMerge = []
    for iterfile in iterfiles:

        for inFolder in inputFolders:
            inFileName = glob.glob(inFolder+"/"+iterfile)[0]
            filesToMerge.append(inFileName)

        outFileName = outdir+"/"+inFileName.split("/")[-1]

        print("Merging", filesToMerge, "--->", outFileName)

        subprocess.run(["hadd", "-f", outFileName]+filesToMerge, check=True)

    print("---merged outputs done---")


def main():

    parser = getArgumentParser()
    args = parser.parse_args()
    haddOutputs(args.indir,
                args.regex,
                args.outdir)


if __name__ == "__main__":
    main()
