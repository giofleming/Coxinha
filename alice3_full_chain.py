#!/usr/bin/env python3

import importlib.util
from job_configs import ChainConfig

import alice3.performance.simulation as alice3_simulation
import alice3.performance.writers as alice3_writers
import alice3.performance.generator as alice3_generator
import alice3.performance.alignment as alice3_alignment

import os
import argparse
import sys
from acts.examples.reconstruction import (
    TrackSelectorConfig,
    addVertexFitting,
    VertexFinder,
)

from acts.examples.simulation import (
    addGenParticleSelection,
    addSimParticleSelection,
    addFatras,
    ParticleSelectorConfig,
    addDigitization,
    addDigiParticleSelection,
)

from acts.examples.root import (
    RootParticleWriter,
    RootParticleReader,
    RootSimHitReader,
    RootVertexReader,
    RootTrackSummaryWriter,
)

import pathlib
import acts
import acts.examples
import acts.examples.geant4


# To automatically unzip stuff
from zipfile import ZipFile

# Move it to an utility tool


def unzipFile(zipfile: pathlib.Path):
    if zipfile.exists():
        # unzip
        with ZipFile(zipfile, 'r') as zip_ref:
            zip_ref.extractall(zipfile.parent)
            print(f"Extracted {zipfile}...")
    else:
        raise FileNotFoundError(f"{zipfile} doesn't exist!")


u = acts.UnitConstants


def load_config(path: str):

    config_dir = str(pathlib.Path(path).parent.resolve())
    if config_dir not in sys.path:
        sys.path.insert(0, config_dir)

    spec = importlib.util.spec_from_file_location("job_config", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    return mod.cfg


def print_cfg(cfg):

    import dataclasses
    from pprint import pprint

    for name, value in vars(type(cfg)).items():
        if name.startswith("_"):
            continue
        if dataclasses.is_dataclass(value):
            print(f"--- {name} ---")
            pprint(dataclasses.asdict(value))
            print()


def setMagneticField(cfg):
    """Set MF using a map. If no map provided - set a constant MF along z"""
    if cfg.general.fieldMap != "":
        print(">>> !", cfg.general.fieldMap, "field map is used in this sim+rec run.")
        with open(cfg.general.fieldMap) as magFile:  # Checking if it's RZ or XYZ coordinates
            for l in magFile:
                l = l.strip()
                if l.startswith("#"):
                    continue
                while "  " in l:
                    l = l.replace("  ", " ")
                l = l.split(" ")
                if len(l) == 4:
                    field = acts.MagneticFieldMapRz(cfg.general.fieldMap)
                else:
                    field = acts.MagneticFieldMapXyz(cfg.general.fieldMap)
                break
    else:
        print(">>> !Using Constant B-Field")
        field = acts.ConstantBField(acts.Vector3(0.0, 0.0, cfg.general.MF * u.T))
    return field


def applyMaterialMapGen3(trackingGeometry, geo_dir):
    """Gen3 only: apply the material map AFTER construction.

    Gen3 has no material-decorator hook (construct() takes none), so instead of
    the Gen1 `buildALICE3Geometry(..., material=True)` argument the map is loaded
    by walking the geometry and handing every surface to the decorator. Reads
    geo_dir/material_map_gen3.json (produced by material_mapping.py
    --geomGeneration gen3).
    """
    matFile = geo_dir / "material_map.json"
    if not matFile.exists():
        unzipFile(geo_dir / "material_map.json.zip")

    decorator = acts.IMaterialDecorator.fromFile(str(matFile))

    class _MaterialLoader(acts.TrackingGeometryMutableVisitor):
        def __init__(self):
            super().__init__()
            self.nWithMaterial = 0

        def visitSurface(self, surface):
            decorator.decorate(surface)
            if surface.surfaceMaterial is not None:
                self.nWithMaterial += 1

    loader = _MaterialLoader()
    trackingGeometry.apply(loader)
    print(f">>> Gen3 material: {loader.nWithMaterial} surfaces carry material "
          f"from {matFile}")
    if loader.nWithMaterial == 0:
        raise RuntimeError(
            f"{matFile} matched no surface in this geometry. Geometry IDs "
            "change whenever the blueprint structure changes - regenerate the "
            "map with material_mapping.py --geomGeneration gen3."
        )


def addLutMakerAlgorithm(
    s: acts.examples.Sequencer,
    *,
    inputTracks: str,
    inputParticles: str,
    inputTrackParticleMatching: str,
    magneticField: float,
    lutTag: str,
    useResiduals: bool = True,
    storeCovarianceSigma: bool = False,
    produceOutputTree: bool = True,
    logLevel=None,
    outputDirectory: str = "luts"
):
    """Add the custom LUT maker to the live reconstruction chain.

    The input track collection is the final ACTS ConstTrackContainer. The
    accompanying LutMaker.hpp/LutMaker.cpp patch in this bundle changes the
    algorithm input handle accordingly and performs the lookup with
    track.index().
    """
    if not inputTrackParticleMatching:
        raise ValueError("LutMaker requires an explicit track-particle matching collection")

    from AliceActsPythonBindings import LutMakerAlgorithm

    lutMakerCfg = LutMakerAlgorithm.Config()
    lutMakerCfg.inputTracks = inputTracks
    lutMakerCfg.inputParticles = inputParticles
    lutMakerCfg.inputTrackParticleMatching = inputTrackParticleMatching
    lutMakerCfg.lutTag = lutTag
    lutMakerCfg.produceOutputTree = produceOutputTree
    lutMakerCfg.useResiduals = useResiduals
    lutMakerCfg.magneticField = float(magneticField)
    lutMakerCfg.storeCovarianceSigma = storeCovarianceSigma
    lutMakerCfg.outputDirectory = os.path.abspath(outputDirectory)
    customLogLevel = acts.examples.defaultLogging(s, logLevel)
    s.addAlgorithm(
        LutMakerAlgorithm(
            config=lutMakerCfg,
            level=customLogLevel(),
        )
    )


def addLutProduction(
    s: acts.examples.Sequencer,
    *,
    cfg,
    outputDir: pathlib.Path,
    finalTracks: str,
    matchingRatio: float = 0.5,
    useResiduals: bool = True,
    writeTruthSummary: bool = True,
    outputDirectory: str = "luts",
):
    if not finalTracks:
        raise ValueError("The final track collection name must not be empty")

    selectedParticles = "particles_digitized_selected"
    trackParticleMatching = "lut_track_particle_matching"
    particleTrackMatching = "lut_particle_track_matching"

    print(f"LUT maker input tracks: {finalTracks}")

    s.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=acts.logging.INFO,
            inputTracks=finalTracks,
            inputParticles=selectedParticles,
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching=trackParticleMatching,
            outputParticleTrackMatching=particleTrackMatching,
            matchingRatio=matchingRatio,
            doubleMatching=False,
        )
    )

    if writeTruthSummary:
        s.addWriter(
            RootTrackSummaryWriter(
                level=acts.logging.INFO,
                inputTracks=finalTracks,
                inputParticles=selectedParticles,
                inputTrackParticleMatching=trackParticleMatching,
                filePath=str(
                    outputDir /
                    "tracksummary_ambi-tracks-merged-with-truth.root"
                ),
                writeCovMat=True,
            )
        )

    addLutMakerAlgorithm(
        s,
        inputTracks=finalTracks,
        inputParticles=selectedParticles,
        inputTrackParticleMatching=trackParticleMatching,
        magneticField=cfg.general.MF,
        lutTag=cfg.general.geo_dir.rstrip("/").split("/")[-2],
        outputDirectory=outputDirectory,
        useResiduals=useResiduals,
        storeCovarianceSigma=False,
        produceOutputTree=True,
        logLevel=acts.logging.INFO,
    )


def getArgumentParser():
    """Get arguments from command line"""
    parser = argparse.ArgumentParser(
        description="Command line arguments for full chain setup"
    )

    # General parameters
    parser.add_argument(
        "--config",
        "-c",
        dest="config",
        help="Job config",
        type=str,
        default="./job_configs/pp.py"
    )
    parser.add_argument(
        "--nThreads",
        "-nthr",
        dest="nThreads",
        help="Number of threads",
        type=int,
        default=-1
    )
    parser.add_argument(
        "--nEv",
        "-n",
        dest="nEvents",
        help="nEvents",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--skip",
        "-s",
        dest="skip",
        help="Skip events",
        type=int,
        default=0,
    )

    parser.add_argument(
        "--out_dir_prefix",
        "-o",
        dest="out_dir_prefix",
        help="Dir for output",
        type=str,
        default="test",
    )

    parser.add_argument(
        "--full_output_path",
        "-f",
        dest="full_output_path",
        help="Full path for output",
        type=str,
        default="",
    )

    parser.add_argument(
        "--seed",
        "--randomSeed",
        dest="seed",
        help="seed to use in the random number generator",
        type=int,
        default=42)

    parser.add_argument(
        "--inGenParticlesFile",
        dest="inGenParticlesFile",
        help="Input particles file for loading generated particles",
        type=str,
        default="")

    parser.add_argument(
        "--inParticlesFile",
        dest="inParticlesFile",
        help="Input particles file for loading simulated particles",
        type=str,
        default="")

    parser.add_argument(
        "--inSimHitsFile",
        dest="inSimHitsFile",
        help="Input sim-hits file for loading simulated hits",
        type=str,
        default="")

    parser.add_argument(
        "--inVtxsFile",
        dest="inVtxsFile",
        help="Input vertices file for loading truth vertices/collisions",
        type=str,
        default="")

    parser.add_argument(
        "--simdir",
        dest="simdir",
        help="Location of the simulation output directory",
        type=str,
        default="")

    parser.add_argument(
        "--runMode",
        dest="runMode",
        help="Running mode: sim_only, reco_only, full_chain",
        type=str,
        choices=["sim_only", "reco_only", "full_chain"],
        default="full_chain",
    )

    parser.add_argument(
        "--enableLut", "--enable-lut", "--enable_lut",
        action="store_true",
        help="Truth-match final tracks and produce LUT files",
    )

    #####
    return parser


def runRecoChain(s,
                 cfg,
                 geo_dir,
                 rnd,
                 field,
                 trackingGeometry,
                 outputDir,
                 enableLut):

    # Run Digi
    s = addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=cfg.general.digi_file,
        outputDirRoot=outputDir if cfg.digi.storeMeasurements else None,
        rnd=rnd,
        logLevel=acts.logging.INFO,
        doMerge=True,   # controls cross-hit cluster merging in the geometric digitization, and its effective default is true
        mergeCommonCorner=True,  # for clustering with segmented digi, 8-cell connectivity
        minEnergyDeposit=1.0 * u.keV,  # 3.65 keV (1000 e/h pairs × 3.65 eV per e/h pair in Si, PDG 2023 Table 35.10)
    )

    if cfg.digi.injectPerPixelNoise:
        alice3_simulation.addPerPixelNoise(
            s,
            trackingGeometry,
            digiConfigFile=cfg.general.digi_file,
            # layersByVolume=None, # {volumeId: [layerIds]}; [] => all layers; None/{} => no-op
            # layersByVolume={14: [],  19: [], 22: [],   18: [], 21: [],   5: [], 17: []}, # gen 1
            layersByVolume={20: [],  30: [], 40: [],   50: [], 60: [],   70: []},  # gen 3
            noiseProbPerPixel=cfg.digi.noiseProbPerPixel,       # fake-hit probability per pixel per event
            # time window for pure-noise time (0 => t=0); only used if digi smears time.
            # Matches the collision spread: U(0, nBC*spacing) with discrete BC times, t=0 otherwise
            noiseTimeWindow=(cfg.pythia.bunchSpacing * cfg.pythia.nBunchCrossings * u.ns
                             if cfg.pythia.useDiscreteBCtime else 0.0),
            # noiseActivation=1.0, # activation value assigned to a noise cell
            # rnd=None,
            # mergeCommonCorner=True,       # clustering parameter (mirror the digitizer)
            # mergeNsigma=1.0,              # clustering parameter (mirror the digitizer)
            # inputSuffix="_raw",           # the digitizer must emit "<name><inputSuffix>"
            # outputSuffix="_with_noise",   # the digitizer must emit "<name><outputSuffix>"
            # logLevel=None,
            outputDirRoot=outputDir if cfg.digi.storeNoiseMeasurements else None,
        )

    # Digitization writes:
    #   "measurements"       -> MeasurementContainer
    #   "measurement_subset" -> MeasurementSubset containing all measurements
    #
    # SpacePointMaker, ProtoTracksToTracks and TrackFinding should use
    # MeasurementSubset keys, not the raw "measurements" container.
    initialMeasurementSubset = "measurement_subset"

    alice3_simulation.addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            eta=(-cfg.detSim.particleSelectionEta,
                 cfg.detSim.particleSelectionEta),
            pt=(20 * u.MeV, None),
            removeNeutral=True,
            removeSecondaries=True,  # May 2026
            hits=(cfg.digi.particleSelectionMinHits, None),
            rho=tuple(x * u.mm for x in cfg.digi.digi_particleSelectionRho),
            absZ=tuple(x * u.mm for x in cfg.digi.digi_particleSelectionZ),
        ),
        measurementLayers=cfg.digi.particleSelectionLayers,
    )

    if (cfg.digi.storeDigiParticles):
        s.addWriter(
            RootParticleWriter(
                level=acts.logging.INFO,
                inputParticles="particles_digitized_selected",
                bField=field,
                writeHelixParameters=False,
                filePath=str(outputDir / "particles_digitized_selected.root"),
            )
        )

    ###### Reconstruction ######

    # ---- optional independent "pre-vertexing" pass -------------------------
    # A single seeding+CKF+ambi pass followed by vertex finding, used ONLY to
    # localize the primary vertices in z. The found vertex z-positions (+ their
    # uncertainties) are then handed to the main iterative reco to constrain the
    # seeding (see addIterativeTracking(inputVertices=...) + the GridTriplet
    # experimentCuts vertex-z window). It is fully namespaced ("prevtx_") so it
    # does not interfere with the main reconstruction.
    preVertexCollection = ""
    if getattr(cfg.general, "doPreVertexing", False) and cfg.general.iterations > 0:
        prevtxDir = pathlib.Path(outputDir) / "prevtx"
        prevtxDir.mkdir(parents=True, exist_ok=True)
        alice3_writers.addIterativeTracking(
            s,
            cfg=cfg,
            trackingGeometry=trackingGeometry,
            geo_dir=geo_dir,
            field=field,
            iterations=1,                        # one pass is enough to find vertices
            inputTracks="tracks",
            inputMeasurementSubset=initialMeasurementSubset,
            outputDir=prevtxDir,
            namePrefix="prevtx_",                # namespace all collections/outputs
            preVertexing=True,                   # use the lighter dedicated seeding config
        )
        # Vertex finding on the pre-pass tracks, built DIRECTLY (not via
        # addVertexFitting) so every collection is "prevtx_"-namespaced and no
        # VertexTruthMatcher is pulled in (it would read the global
        # "track_particle_matching", which the single-pass pre-reco never
        # creates). The main reco then consumes ONLY the resulting prevtx_vertices.
        acts.examples.reconstruction.addTrackSelection(
            s,
            TrackSelectorConfig(eta=(-1.7, 1.7), pt=(0.4, 100.0), nMeasurementsMin=7),
            inputTracks="prevtx_ambi_tracks",
            outputTracks="prevtx_selectedTracksVertexing",
            logLevel=acts.logging.INFO,
        )
        s.addAlgorithm(acts.examples.TracksToParameters(
            level=acts.logging.INFO,
            inputTracks="prevtx_selectedTracksVertexing",
            outputTrackParameters="prevtx_selectedTracksParametersVertexing",
        ))
        s.addAlgorithm(acts.examples.AdaptiveMultiVertexFinderAlgorithm(
            level=acts.logging.INFO,
            inputTrackParameters="prevtx_selectedTracksParametersVertexing",
            outputProtoVertices="prevtx_protovertices",
            outputVertices="prevtx_vertices",
            bField=field,
            seedFinder=acts.examples.VertexSeedFinder.AdaptiveGridSeeder,
            **acts.examples.defaultKWArgs(
                useTime=True,
                spatialBinExtent=0.05 * u.mm,
                temporalBinExtent=5.0 * u.ns,
            ),
        ))
        preVertexCollection = "prevtx_vertices"

    if cfg.general.iterations > 0:
        alice3_writers.addIterativeTracking(s,
                                            cfg=cfg,
                                            trackingGeometry=trackingGeometry,
                                            geo_dir=geo_dir,
                                            field=field,
                                            iterations=cfg.general.iterations,
                                            inputTracks="tracks",
                                            inputMeasurementSubset=initialMeasurementSubset,
                                            outputDir=outputDir,
                                            inputVertices=preVertexCollection)

        # The LUT must use the final ambiguity-resolved/merged collection.
        # Truth matching is performed while tracks, measurements and truth maps
        # are all still present on the whiteboard.

        if enableLut:
            addLutProduction(
                s,
                cfg=cfg,
                outputDir=outputDir,
                finalTracks="ambi-tracks-merged",
                matchingRatio=0.5,
                useResiduals=True,
                writeTruthSummary=False,
            )

        # Do track re-fitting with GSF if requested (relevant for electrons)
        if cfg.general.doGSFtrackRefit:
            alice3_writers.addGsfRefit(
                s,
                trackingGeometry,
                field,
                inputTracks="tracks",   # the final ambiguity-resolved/merged collection
                maxComponents=12,
                weightCutoff=1.0e-4,
                componentMergeMethod=acts.examples.ComponentMergeMethod.maxWeight,
                mixtureReductionAlgorithm=acts.examples.MixtureReductionAlgorithm.KLDistance,
                reverseFilteringCovarianceScaling=10,  # 100.0,
                outputDirRoot=outputDir if cfg.general.writeGSFTrackSummary else None,  # tracksummary_gsf_refit.root
                writeStates=False,                      # if we also want trackstates
                logLevel=acts.logging.INFO,
            )

        # Do vertex finding & fitting if requested
        if cfg.general.doVertexReco:
            s = addVertexFitting(
                s,
                field,
                vertexFinder=VertexFinder.AMVF,
                outputDirRoot=outputDir if cfg.general.writeVertexPerformance else None,
                seeder=acts.examples.VertexSeedFinder.AdaptiveGridSeeder,
                writeTrackInfo=False,
                useTime=True,
                trackSelectorConfig=TrackSelectorConfig(
                    eta=(-1.7, 1.7),                   # |eta| window for tracks
                    pt=(0.4, 100.0),
                    # time=(-0.2 * u.ns, 0.2 * u.ns),  # time window
                    nMeasurementsMin=7,                # minimum number of measurements
                ),
                tracksMaxZinterval=1.0 * u.mm,              # z window for adding tracks to a vertex seed
                # spatialWindow=(250 * u.mm, 250 * u.mm),   # (min,max) z range of the seeding density map; should be fine enough to separate neighbouring vertices
                spatialBinExtent=0.05 * u.mm,               # z granularity (default 15 µm is also fine)
                # [this python param needs porting in ACTS!!!] wide time range to cover many BCs; default ±10 ns would clip almost everything
                temporalWindow=(-50 * u.ns, (50 + cfg.pythia.bunchSpacing * (cfg.pythia.nBunchCrossings-1)) * u.ns),
                temporalBinExtent=(5. * u.ns),     # t bin width of the seeding grid
                # initialVariances=acts.Vector4(10, 10, 500, (30 * u.ns)**2), # priors (data-driven; just be generous, especially in t)  # x,y,z ~ (10 mm)², t ~ (100 ns)²
                maxIterations=20,
                # logLevel = acts.logging.VERBOSE,
            )

            # AliceVertexTruthMatcher for vertex truth matching with vertices from many bunch crossings
            alice3_writers.addVertexPostProcessing(
                s,
                field=field,
                outputDir=outputDir,
                writePerformance=cfg.general.writeVertexPerformance,
            )

    return s


def addSimReadersForReco(s, simdir):
    """Add root readers needed when running reco from existing sim output."""

    # s.addReader(
    #    RootParticleReader(
    #        level=acts.logging.WARNING,
    #        outputParticles="particles_generated_selected",
    #        filePath=simdir / "particles.root",
    #    )
    #
    # )

    s.addReader(
        RootVertexReader(
            level=acts.logging.WARNING,
            outputVertices="vertices_truth",
            filePath=simdir / "vertices.root",
        )
    )

    s.addReader(
        RootParticleReader(
            level=acts.logging.WARNING,
            outputParticles="particles_simulated_selected",
            filePath=simdir / "particles_simulation.root",
        )
    )

    s.addWhiteboardAlias("particles", "particles_simulated_selected")
    # needed for truth matching with Geant4 secondaries
    s.addWhiteboardAlias("particles_simulated", "particles_simulated_selected")
    # will be rewritten by addDigiParticleSelection
    s.addWhiteboardAlias("particles_selected", "particles_simulated_selected")

    s.addReader(
        RootSimHitReader(
            level=acts.logging.WARNING,
            outputSimHits="simhits",
            treeName="hits",
            filePath=simdir / "hits.root",
        )
    )
    return s


def runSimChain(s: acts.examples.Sequencer = None,
                cfg=None,
                geo_dir=None,
                rnd=None,
                field=None,
                trackingGeometry=None,
                outputDir=None
                ):

    ###### Simulation ######

    # generator and Fatras output dir, None = no output
    simOutputDir = outputDir if cfg.detSim.writeSimOutput else None

    if not cfg.pythia.usePythia:
        alice3_generator.addAlice3ParticleGun(s, rnd=rnd, cfg=cfg,
                                              outputDir=simOutputDir)
    else:
        s = alice3_generator.addAlice3Pythia8(s, rnd=rnd, cfg=cfg,
                                              outputDir=simOutputDir)

    addGenParticleSelection(
        s,
        ParticleSelectorConfig(
            eta=(-cfg.detSim.particleSelectionEta,
                 cfg.detSim.particleSelectionEta),
            pt=(0.02 * u.MeV, None),
            removeNeutral=False,
            rho=(0.0, 24 * u.mm),
            # absZ=(0.0, 1.0 * u.m),
        ),
    )

    if cfg.detSim.simulation == "Fatras":

        addFatras(
            s,
            trackingGeometry,
            field,
            enableInteractions=True,
            rnd=rnd,
            pMin=cfg.detSim.minFatrasPt,  # GeV, May 2025
            outputDirRoot=simOutputDir,
            logLevel=acts.logging.INFO,
        )

    elif cfg.detSim.simulation == "Geant4":

        alice3_gdml = geo_dir / "o2sim_geometry.gdml"

        # Check if we need to unzip
        if not alice3_gdml.exists():
            unzipFile(geo_dir / "o2sim_geometry.gdml.zip")

        gdml_detector = acts.examples.geant4.GdmlDetector(path=str(alice3_gdml))

        alice3_simulation.addAlice3Geant4(
            s,
            gdml_detector,
            trackingGeometry,
            field,
            # materialMappings=["TRK_SILICON", "FT3_SILICON", "FT3_Silicon", "TF3_Silicon"],
            volumeMappings=["TRKSensor", "FT3Sensor", "FT3Sensor_Active",
                            "TOFTRKSensor", "TOFSensor", "radTile", "RICHPhotoTile"],
            outputDirRoot=outputDir,
            rnd=rnd,
            logLevel=acts.logging.INFO,
            killVolume=trackingGeometry.highestTrackingVolume,
            # allows simulation of many bunch crossings with G4 (e.g. 500 ns for 20 BCs with 25 ns spacing),
            killAfterTime=40. * u.ns if cfg.pythia.useDiscreteBCtime == False else (
                40 + cfg.pythia.bunchSpacing * (cfg.pythia.nBunchCrossings - 1)) * u.ns,
            recordHitsOfSecondaries=True,
            # IA: no need, to save space; set to True for secondary vertexing from G4 and V0/cascade studies
            keepParticlesWithoutHits=cfg.detSim.keepParticlesWithoutHits,
            killSecondaries=False,
            propagatorLargestAcceptableStep=150,
            writeParticles=cfg.detSim.writeSimOutput,
            writeHits=cfg.detSim.writeSimOutput,
        )

        # if secondary vertices from Geant (with id-s of mothers and outcoming particles) are needed (e.g. for V0 analysis):
        if cfg.detSim.storeGenAndGeantVertices:
            alice3_writers.addGenAndGeantVertices(
                s,
                inputParticles="particles_simulated",   # FULL simulated set (incl. G4 secondaries)
                outputVertices="vertices_gen_and_geant",
                outputDirRoot=outputDir,                 # writes vertices_gen_and_geant.root
                logLevel=acts.logging.INFO,
            )

    addSimParticleSelection(
        s,
        ParticleSelectorConfig(
            eta=(-cfg.detSim.particleSelectionEta, cfg.detSim.particleSelectionEta),
            pt=(10 * u.MeV, None),
            removeNeutral=True,  # May 2025
            hits=(5, None),  # Nov 2025
        ),
    )

    # IA: better to store all particles from Geant ("particles_simulation.root"), to have full information about origins of hits
    # s.addWriter(
    #         RootParticleWriter(
    #             level=acts.logging.INFO,
    #             inputParticles="particles_simulated_selected",
    #             bField=field,
    #             writeHelixParameters=False,
    #             filePath=str(outputDir / "particles_simulated_selected.root"),
    #         )
    #     )

    return s


def runFullChain(cfg=None, args=None):
    """
    Run the full ALICE3 simulation and reconstruction chain.

    Parameters
    ----------
    cfg : ChainConfig.Config
        Configuration object for the chain
    args : argparse.Namespace
        Command line arguments (same format as getArgumentParser().parse_args())
    """
    # Geometry configuration in case cfg not set
    if (cfg is None):
        cfg = ChainConfig.Config()

    if (args is None):
        parser = getArgumentParser()
        args = parser.parse_args()

    geo_dir = pathlib.Path(os.path.abspath(os.path.dirname(__file__))) / cfg.general.geo_dir
    print("Loading geometry from ... ", str(geo_dir))
    if cfg.general.geomGeneration == "gen1":
        sys.path.insert(0, str(pathlib.Path(geo_dir).resolve()))
        import buildALICE3Geometry
    # gen3: the builder is compiled into AliceActsPythonBindings (imported in the
    # gen3 branch below); no sys.path juggling / buildALICE3Geometry needed.

    ##########################
    # SOME OTHER PARAMS
    ##########################

    # output directory
    if args.full_output_path != "":
        IA_outputDirName = pathlib.Path("output") / args.full_output_path
    else:
        IA_outputDirName = (
            "output/"
            + args.out_dir_prefix
            + "_nEv"
            + str(args.nEvents)
            + "_PID"
            + str(cfg.particleGun.gunPID)
        )

    outputDir = pathlib.Path.cwd() / IA_outputDirName

    if not outputDir.exists():
        outputDir.mkdir(mode=0o777, parents=True, exist_ok=True)

    if cfg.general.geomGeneration == "gen3":
        # Gen3 (Blueprint): build straight from the TGeo file. withMaterial=True
        # keeps the passive/proto-material surfaces the map has entries for; the
        # real material is loaded afterwards by applyMaterialMap().
        from AliceActsPythonBindings import (
            buildALICE3GeometryGen3,
            Gen3NavigationPolicy,
        )
        navPolicy = {
            "tryAll": Gen3NavigationPolicy.TryAll,
            "surfaceArray": Gen3NavigationPolicy.SurfaceArray,
        }[getattr(cfg.general, "navigationPolicy", "tryAll")]
        trackingGeometry = buildALICE3GeometryGen3(
            str(geo_dir / "o2sim_geometry.root"), acts.logging.INFO, "", True,
            "", navPolicy
        )
        if cfg.general.enableMaterial:
            applyMaterialMapGen3(trackingGeometry, geo_dir)
    else:
        detector = buildALICE3Geometry.buildALICE3Geometry(
            geo_dir, cfg.general.enableMaterial, False, acts.logging.INFO
        )
        trackingGeometry = detector.trackingGeometry()
        # decorators = detector.contextDecorators() #Unused

    field = setMagneticField(cfg)

    rnd = acts.examples.RandomNumbers(seed=args.seed)

    s = acts.examples.Sequencer(events=args.nEvents,
                                skip=args.skip,
                                numThreads=args.nThreads,
                                outputDir=outputDir,
                                outputTimingFile="timing.csv")

    common_cfg = dict(cfg=cfg,
                      geo_dir=geo_dir,
                      rnd=rnd,
                      field=field,
                      trackingGeometry=trackingGeometry,
                      outputDir=outputDir
                      )

    run_sim = args.runMode in ("sim_only", "full_chain")
    run_reco = args.runMode in ("reco_only", "full_chain")

    if cfg.alignment.testMisalignment:
        if not run_reco:
            print("== WARNING == testMisalignment only affects reconstruction: nothing to misalign in sim_only")
        else:
            alice3_alignment.addGlobalShifts(s,
                                             trackingGeometry,
                                             cfg.alignment.globalShifts,
                                             outputDir=outputDir)

    print("INFO:: run_sim = ", run_sim, "run_reco=", run_reco)
    if run_sim:
        s = runSimChain(s, **common_cfg)

    if run_reco:

        if not run_sim:  # reco_only: need to read sim output from disk

            simdir = None if args.simdir is None else pathlib.Path(args.simdir)
            s = addSimReadersForReco(s, simdir)

        s = runRecoChain(s, **common_cfg, enableLut=args.enableLut)

    s.run()

    if run_reco:
        alice3_writers.print_track_state_creator_statistics()

    # save this script to the output directory
    outputDirForConfig = outputDir / "configs"
    if not outputDirForConfig.exists():
        outputDirForConfig.mkdir()
    import shutil
    shutil.copyfile(sys.argv[0], outputDirForConfig / pathlib.Path(sys.argv[0]).name)
    shutil.copytree("job_configs", outputDirForConfig / "job_configs", dirs_exist_ok=True)
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


if __name__ == "__main__":
    parser = getArgumentParser()
    args = parser.parse_args()
    cfg = load_config(args.config)

    print_cfg(cfg)

    runFullChain(cfg=cfg, args=args)
