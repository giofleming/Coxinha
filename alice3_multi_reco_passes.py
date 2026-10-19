import os
import sys
import pathlib
import acts

from job_configs import ChainConfig
from alice3_full_chain import runRecoChain, getArgumentParser, addSimReadersForReco, load_config, print_cfg, setMagneticField


def runMultiRecoPasses(cfg=None, args=None):
    """
    Runs the reconstruction chain multiple times with changed configurations.
    """
    # Geometry configuration in case cfg not set
    if (cfg is None):
        cfg = ChainConfig.Config()

    if (args is None):
        parser = getArgumentParser()
        args = parser.parse_args()

    # loading geometry only ONCE for all reco runs
    geo_dir = pathlib.Path(os.path.abspath(os.path.dirname(__file__))) / cfg.general.geo_dir
    print("Loading geometry from ... ", str(geo_dir))
    sys.path.insert(0, str(pathlib.Path(geo_dir).resolve()))

    import buildALICE3Geometry

    detector = buildALICE3Geometry.buildALICE3Geometry(
        geo_dir, cfg.general.enableMaterial, False, acts.logging.INFO
    )
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()
    field = setMagneticField(cfg)
    rnd = acts.examples.RandomNumbers(seed=args.seed)

    # sim directory (in this script, we run reco only)
    simdir = None if args.simdir is None else pathlib.Path(args.simdir)

    #################################
    # establish parameter ranges:
    #################################
    arr_chi2meas = (10, 15, 20)  # example

    # loop over the parameter range and run reco multiple times
    saveConfigs = True
    for chi2meas in arr_chi2meas:
        # apply changes to the default config (before it's used later):
        cfg.tracking.ckfChi2Measurement = chi2meas

        # output directory
        outputDirName = (
            "output/"
            + args.out_dir_prefix
            + "_nEv"
            + str(args.nEvents)
            + "_nMeasMin"
            + str(cfg.tracking.nMeasurementsMin)
            + "_nOutliers"
            + str(cfg.tracking.nOutliers)
            + "_nHoles"
            + str(cfg.tracking.nHoles)
            + "_ckfChi2Meas"
            + str(cfg.tracking.ckfChi2Measurement)
            # + "_ckfMeasPerSurf"
            # + str(cfg.tracking.ckfMeasPerSurf)
        )
        outputDir = pathlib.Path.cwd() / outputDirName
        if not outputDir.exists():
            outputDir.mkdir(mode=0o777, parents=True, exist_ok=True)

        common_cfg = dict(cfg=cfg,
                          geo_dir=geo_dir,
                          rnd=rnd,
                          field=field,
                          trackingGeometry=trackingGeometry,
                          outputDir=outputDir)

        s = acts.examples.Sequencer(events=args.nEvents,
                                    skip=args.skip,
                                    numThreads=args.nThreads)
        s = addSimReadersForReco(s, simdir)
        s = runRecoChain(s, **common_cfg)
        s.run()

        # save this script to the output directory
        if saveConfigs:
            outputDirForConfig = outputDir / "configs"
            if not outputDirForConfig.exists():
                outputDirForConfig.mkdir()
            import shutil
            shutil.copyfile(sys.argv[0], outputDirForConfig / pathlib.Path(sys.argv[0]).name)
            shutil.copyfile("job_configs/ChainConfig.py", outputDirForConfig / "ChainConfig.py")
            shutil.copyfile("alice3_full_chain.py", outputDirForConfig / "alice3_full_chain.py")
            shutil.copyfile("alice3/performance/generator.py", outputDirForConfig / "generator.py")
            shutil.copyfile("alice3/performance/seeding.py", outputDirForConfig / "seeding.py")
            shutil.copyfile("alice3/performance/simulation.py", outputDirForConfig / "simulation.py")
            shutil.copyfile("alice3/performance/writers.py", outputDirForConfig / "writers.py")

            # save command line used
            from datetime import datetime
            cmd = " ".join(sys.argv)
            with open(outputDirForConfig / "script_command_line_args_used.log", "a") as f:
                f.write(f"{datetime.now()}\n")
                f.write(f"{cmd}\n")
            saveConfigs = False  # save only once for the first reco run


if __name__ == "__main__":
    parser = getArgumentParser()
    args = parser.parse_args()
    cfg = load_config(args.config)
    # print_cfg(cfg)

    runMultiRecoPasses(cfg=cfg, args=args)
