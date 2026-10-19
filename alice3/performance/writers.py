from pathlib import Path
from typing import Optional, Union, List, Literal, Tuple
from enum import Enum
from collections import namedtuple
import pathlib

import acts
import acts.examples
import acts.examples.reconstruction as acts_reco
from acts import UnitConstants as u


from acts.examples.reconstruction import (
    TrackSelectorConfig,
    CkfConfig,
    #    addAmbiguityResolution,
    AmbiguityResolutionConfig,
)

# import config for default params
from job_configs import ChainConfig

# Alice3 specific algorithms
from AliceActsPythonBindings import TrackTruthMatcher
from AliceActsPythonBindings import TrackMergerAlgorithm
from AliceActsPythonBindings import AliceVertexTruthMatcher
from AliceActsPythonBindings import SimVertexBuilder
from AliceActsPythonBindings import SeedAnalysisAlg
from AliceActsPythonBindings import RootUnbiasedResidualWriter
from AliceActsPythonBindings import AliceTrackFindingAlgorithm

# Alice3 seeding
import alice3.performance.seeding as alice3_seeding

# Alice3 plotting
import alice3.performance.plotting as alice3_plotting


def get_track_selector_config(cfg: ChainConfig.Config,
                              iteration: int = 0,
                              **overrides) -> TrackSelectorConfig:
    """
    Get track selector configuration for a given iteration.
    Args:
        cfg: full Config object
        iteration: Iteration number
        overrides: kwargs to override iteration-specific defaults
    """

    tc = cfg.tracking

    if iteration not in tc.iterationParams:
        raise ValueError(f"No trackSelector iterationParams defined for iteration {iteration}")

    it = tc.iterationParams[iteration]

    # Convert to GeV (Shouldn't be needed as u.GeV = 1 in ACTS)

    pt_lo, pt_hi = it['pt']
    pt = (pt_lo * u.GeV if pt_lo is not None else None,
          pt_hi * u.GeV if pt_hi is not None else None)

    params = {
        "nMeasurementsMin": it['nMeasurementsMin'],
        "maxHoles": it['nHoles'],
        "maxOutliers": it['nOutliers'],
        "maxSharedHits": tc.maxSharedHits,
        "pt": pt,
    }

    params.update(overrides)

    return TrackSelectorConfig(**params)


# factories made by make_track_state_creator_factory, kept to report their
# statistics at the end of the run
_trackStateCreatorFactories = []


def make_track_state_creator_factory(cfg: ChainConfig.Config, name: str = "", forceSeedMeasurements: bool = False):
    """
    CKF track state creator factory selected by cfg.tracking.trackStateCreator.
    Returns None for the ACTS default creator.
    forceSeedMeasurements: force the seed measurements in this CKF instance
    (from cfg.tracking.iterationParams[i]['forceSeedMeasurements']).
    """
    tc = cfg.tracking
    if hasattr(tc, "forceSeedMeasurements"):
        raise ValueError("cfg.tracking.forceSeedMeasurements was replaced by "
                         "cfg.tracking.iterationParams[i]['forceSeedMeasurements']")
    if hasattr(tc, "trackFindingAlgorithm"):
        raise ValueError("cfg.tracking.trackFindingAlgorithm was removed: the CKF always runs in "
                         "AliceTrackFindingAlgorithm (identical tracks to the ACTS examples one)")
    creator = getattr(tc, "trackStateCreator", "default")
    if creator == "default":
        if forceSeedMeasurements:
            raise ValueError(f"{name}: forceSeedMeasurements needs cfg.tracking.trackStateCreator = 'gridded'")
        return None
    if creator != "gridded":
        raise ValueError(f"Unknown cfg.tracking.trackStateCreator '{creator}'")
    if tc.stayOnSeed:
        # with stayOnSeed the window can change which candidates are selected
        raise ValueError("The gridded track state creator is not exact with stayOnSeed = True; "
                         "use iterationParams[i]['forceSeedMeasurements'] instead")

    # imported here so that builds without it still run the default creator
    from AliceActsPythonBindings import GriddedTrackStateCreatorFactory

    factoryCfg = GriddedTrackStateCreatorFactory.Config()
    factoryCfg.minMeasurements = tc.trackStateCreatorGridMinMeasurements
    factoryCfg.measurementsPerCell = tc.trackStateCreatorGridMeasurementsPerCell
    factoryCfg.forceSeedMeasurements = forceSeedMeasurements
    level = getattr(acts.logging, getattr(tc, "trackStateCreatorLogLevel", "INFO"))
    factory = GriddedTrackStateCreatorFactory(factoryCfg, level)
    _trackStateCreatorFactories.append((name, factory))
    return factory


def print_track_state_creator_statistics():
    """Print how many measurements the gridded creators calibrated."""
    for name, factory in _trackStateCreatorFactories:
        st = factory.statistics()
        frac = 100. * st.nCandidates / st.nMeasurementsOnSurface if st.nMeasurementsOnSurface else 0.
        print(f"GriddedTrackStateCreator {name}: {st.nGrids} surface grids, {st.nCalls} calls on gridded "
              f"surfaces, calibrated {st.nCandidates} of {st.nMeasurementsOnSurface} measurements ({frac:.3f}%)")
        if st.nForced:
            print(f"    forced seed measurements: {st.nSeeds} seeds, {st.nForced} forced surfaces, "
                  f"{st.nForcedAboveChi2Measurement} above the measurement chi2 cut, "
                  f"{st.nForcedAboveChi2Outlier} above the outlier cut")


def addCKFTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    trackSelectorConfig: Optional[
        Union[acts_reco.TrackSelectorConfig, List[acts_reco.TrackSelectorConfig]]
    ] = None,
    ckfConfig: acts_reco.CkfConfig = acts_reco.CkfConfig(),
    twoWay: bool = True,
    reverseSearch: bool = False,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrackSummary: bool = True,
    writeTrackStates: bool = False,
    writePerformance: bool = True,
    writeCovMat: bool = False,
    inputMeasurements: str = "measurement_subset",
    inputMeasurementParticlesMap: str = "measurement_particles_map",
    inputParticleMeasurementsMap: str = "particle_measurements_map",
    inputInitialTrackParameters: str = "estimatedparameters",
    inputSeeds: str = "estimatedseeds",
    outputTracks: str = "ckf_tracks",
    truthMatchingRatio: float = 0.5,
    writeMatchingDetails: bool = False,
    trackStateCreatorFactory=None,
    secondPassTrimmedCounts: bool = False,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the seeding

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for ROOT output, None triggers no output
    trackSelectorConfig : TrackSelectorConfig(loc0, loc1, time, eta, absEta, pt, phi, minMeasurements)
        TrackSelector configuration. Each range is specified as a tuple of (min,max).
        Specify as a list(TrackSelectorConfig) for eta-dependent cuts, with binning specified by absEta[1].
        Defaults of no cuts specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackSelector.hpp
    writeTrackSummary : bool, True
        write tracksummary_ckf.root ntuple?
    writeTrackStates : bool, False
        write trackstates_ckf.root ntuple? This can be quite large.
    writePerformance : bool, True
        write performance_fitting_ckf.root and performance_finding_ckf.root ntuples?
    writeCovMat : bool, False
        write covaraiance matrices to tracksummary_ckf.root ntuple?
    truthMatchingRatio : float, 0.5
        the ratio used for truth matching
    trackStateCreatorFactory : TrackFindingAlgorithm.TrackStateCreatorFactory, None
        CKF track state creator, None for the ACTS default
    secondPassTrimmedCounts : bool, False
        see cfg.tracking.secondPassTrimmedCounts
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    tslist = (
        []
        if trackSelectorConfig is None
        else (
            [trackSelectorConfig]
            if isinstance(trackSelectorConfig, acts_reco.TrackSelectorConfig)
            else trackSelectorConfig
        )
    )

    if len(tslist) > 1:
        cutSets = []
        for c in tslist:
            defKW = acts_reco.trackSelectorDefaultKWArgs(c)
            defKW.pop("absEtaMax", None)
            cutSets += [acts.TrackSelector.Config(**defKW)]
    else:
        cutSets = [
            acts.TrackSelector.Config(**acts_reco.trackSelectorDefaultKWArgs(c)) for c in tslist
        ]

    if len(tslist) == 0:
        trkSelCfg = None
    elif len(tslist) == 1:
        trkSelCfg = cutSets[0]
    else:
        trkSelCfg = acts.TrackSelector.EtaBinnedConfig(
            cutSets=cutSets,
            absEtaEdges=[cutSets[0].absEtaMin] + [c.absEta[1] for c in tslist],
        )

    finderKWArgs = dict(
        level=customLogLevel(),
        measurementSelectorCfg=acts.MeasurementSelector.Config(
            [
                (
                    acts.GeometryIdentifier(),
                    (
                        [],
                        [ckfConfig.chi2CutOffMeasurement],
                        [ckfConfig.chi2CutOffOutlier],
                        [ckfConfig.numMeasurementsCutOff],
                    ),
                )
            ]
        ),
        inputMeasurements=inputMeasurements,
        inputInitialTrackParameters=inputInitialTrackParameters,
        inputSeeds=(
            inputSeeds
            if ckfConfig.seedDeduplication or ckfConfig.stayOnSeed or (
                # forcing the seed measurements needs the seeds
                trackStateCreatorFactory is not None and trackStateCreatorFactory.config.forceSeedMeasurements)
            else ""
        ),
        outputTracks=outputTracks,
        findTracks=acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
            trackingGeometry, field, customLogLevel()
        ),
        **acts.examples.defaultKWArgs(
            trackingGeometry=trackingGeometry,
            magneticField=field,
            trackSelectorCfg=trkSelCfg,
            maxSteps=ckfConfig.maxSteps,
            twoWay=twoWay,
            reverseSearch=reverseSearch,
            seedDeduplication=ckfConfig.seedDeduplication,
            stayOnSeed=ckfConfig.stayOnSeed,
            pixelVolumeIds=ckfConfig.pixelVolumes,
            stripVolumeIds=ckfConfig.stripVolumes,
            maxPixelHoles=ckfConfig.maxPixelHoles,
            maxStripHoles=ckfConfig.maxStripHoles,
            trimTracks=ckfConfig.trimTracks,
            constrainToVolumeIds=ckfConfig.constrainToVolumes,
            endOfWorldVolumeIds=ckfConfig.endOfWorldVolumes,
            trackStateCreatorFactory=trackStateCreatorFactory,
        ),
    )
    # the CKF runs in actso2's AliceTrackFindingAlgorithm (a copy of the ACTS
    # examples one with the ALICE 3 options); only the finder function above
    # comes from the ACTS examples
    level = finderKWArgs.pop("level")
    finderCfg = AliceTrackFindingAlgorithm.Config()
    for key, value in finderKWArgs.items():
        setattr(finderCfg, key, value)
    finderCfg.secondPassTrimmedCounts = secondPassTrimmedCounts
    trackFinder = AliceTrackFindingAlgorithm(config=finderCfg, level=level)
    s.addAlgorithm(trackFinder)

    if outputTracks == "ckf_tracks":
        s.addWhiteboardAlias("tracks", trackFinder.config.outputTracks)

    truthMatchCfg = TrackTruthMatcher.Config()
    truthMatchCfg.inputTracks = outputTracks
    truthMatchCfg.inputParticles = "particles_selected"
    truthMatchCfg.inputMeasurementParticlesMap = inputMeasurementParticlesMap
    truthMatchCfg.outputTrackParticleMatching = f"{outputTracks}_particle_matching"
    truthMatchCfg.outputParticleTrackMatching = f"particle_{outputTracks}_matching"
    truthMatchCfg.matchingRatio = truthMatchingRatio
    truthMatchCfg.doubleMatching = False
    truthMatchCfg.looperProtection = True
    truthMatchCfg.loop_absEta = 1.5
    truthMatchCfg.loop_maxPt = 0.8 * u.GeV
    truthMatchCfg.loop_maxParticleHits = 9

    matchAlg = TrackTruthMatcher(
        config=truthMatchCfg,
        level=customLogLevel(),
    )
    s.addAlgorithm(matchAlg)

    if outputTracks == "ckf_tracks":
        s.addWhiteboardAlias(
            "track_particle_matching", truthMatchCfg.outputTrackParticleMatching
        )
        s.addWhiteboardAlias(
            "particle_track_matching", truthMatchCfg.outputParticleTrackMatching
        )

    addTrackWriters(
        s,
        name=trackFinder.config.outputTracks,
        tracks=trackFinder.config.outputTracks,
        inputTrackParticleMatching=truthMatchCfg.outputTrackParticleMatching,
        inputParticleTrackMatching=truthMatchCfg.outputParticleTrackMatching,
        inputParticleMeasurementsMap=inputParticleMeasurementsMap,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeSummary=writeTrackSummary,
        writeStates=writeTrackStates,
        writeFitterPerformance=False,
        writeFinderPerformance=writePerformance,
        writeCovMat=writeCovMat,
        writeMatchingDetails=writeMatchingDetails,
        logLevel=logLevel,
    )

    return s


def addTrackTruthMatcher(
        s: acts.examples.Sequencer,
        inputTracks: str,
        inputParticles: str,
        inputMeasurementParticlesMap: str,
        outputTrackParticleMatching: str,
        outputParticleTrackMatching: str,
        looperProtection: bool = True,
        loop_absEta: float = 1.5,
        loop_maxPt: float = 0.8,
        loop_maxParticleHits: int = 11,
        truthMatchingRatio: float = 0.5,
        logLevel: Optional[acts.logging.Level] = None,
):
    truthMatchCfg = TrackTruthMatcher.Config()
    truthMatchCfg.inputTracks = inputTracks
    truthMatchCfg.inputParticles = inputParticles
    truthMatchCfg.inputMeasurementParticlesMap = inputMeasurementParticlesMap
    truthMatchCfg.outputTrackParticleMatching = outputTrackParticleMatching
    truthMatchCfg.outputParticleTrackMatching = outputParticleTrackMatching
    truthMatchCfg.matchingRatio = truthMatchingRatio
    truthMatchCfg.doubleMatching = False
    truthMatchCfg.looperProtection = looperProtection
    # IA: if doubleMatching is False, no effect from this looper protection (?),
    # https://github.com/acts-project/acts/blob/main/Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackTruthMatcher.cpp#L104
    truthMatchCfg.loop_absEta = loop_absEta
    truthMatchCfg.loop_maxPt = loop_maxPt
    truthMatchCfg.loop_maxParticleHits = loop_maxParticleHits

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    matchAlg = TrackTruthMatcher(
        config=truthMatchCfg,
        level=customLogLevel(),
    )

    s.addAlgorithm(matchAlg)


def addResidualWriters(
    s: acts.examples.Sequencer,
    name: str,
    tracks: str = "tracks",
    outputDirRoot: Optional[Union[Path, str]] = None,
    prefix=None
):

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)

    if not outputDirRoot.exists():
        outputDirRoot.mkdir(parents=True, exists_ok=True)

    # Unbiased residual writer
    rootUnbiasedResidualWriterCfg = RootUnbiasedResidualWriter.Config()
    rootUnbiasedResidualWriterCfg.inputTracks = tracks
    rootUnbiasedResidualWriterCfg.filePath = str(Path(outputDirRoot) / f"{prefix}unbiasedResiduals.root")
    rootUnbiasedResidualWriterCfg.plotToolConfig = alice3_plotting.unbiasedResidualPlotToolConfig
    rootUnbiasedResidualWriter = RootUnbiasedResidualWriter(rootUnbiasedResidualWriterCfg, level=acts.logging.INFO)

    s.addWriter(rootUnbiasedResidualWriter)


def addTrackWriters(
    s: acts.examples.Sequencer,
    name: str,
    tracks: str = "tracks",
    inputTrackParticleMatching: str = "track_particle_matching",
    inputParticleTrackMatching: str = "particle_track_matching",
    inputParticleMeasurementsMap: str = "particle_measurements_map",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeSummary: bool = True,
    writeStates: bool = False,
    writeFitterPerformance: bool = False,
    writeFinderPerformance: bool = False,
    logLevel: Optional[acts.logging.Level] = None,
    writeCovMat: bool = False,
    writeMatchingDetails: bool = False,
):
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir(parents=True, exist_ok=True)

        if writeSummary:
            trackSummaryWriter = acts.examples.root.RootTrackSummaryWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching=inputTrackParticleMatching,
                filePath=str(outputDirRoot / f"tracksummary_{name}.root"),
                treeName="tracksummary",
                writeCovMat=writeCovMat,
            )
            s.addWriter(trackSummaryWriter)

        if writeStates:
            trackStatesWriter = acts.examples.root.RootTrackStatesWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching=inputTrackParticleMatching,
                inputSimHits="simhits",
                inputMeasurementSimHitsMap="measurement_simhits_map",
                filePath=str(outputDirRoot / f"trackstates_{name}.root"),
                treeName="trackstates",
            )
            s.addWriter(trackStatesWriter)

        if writeFitterPerformance:
            trackFitterPerformanceWriter = acts.examples.root.RootTrackParameterPerformanceWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching=inputTrackParticleMatching,
                resPlotToolConfig=alice3_plotting.resPlotToolConfig,
                effPlotToolConfig=alice3_plotting.effPlotToolConfig,
                trackSummaryPlotToolConfig=alice3_plotting.trackSummaryPlotToolConfig,
                filePath=str(outputDirRoot / f"performance_fitting_{name}.root"),
            )
            s.addWriter(trackFitterPerformanceWriter)

        if writeFinderPerformance:
            trackFinderPerfWriter = acts.examples.root.RootPatternRecognitionPerformanceWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching=inputTrackParticleMatching,
                inputParticleTrackMatching=inputParticleTrackMatching,
                inputParticleMeasurementsMap=inputParticleMeasurementsMap,
                label="track",
                effPlotToolConfig=alice3_plotting.effPlotToolConfig,
                fakePlotToolConfig=alice3_plotting.fakePlotToolConfig,
                duplicationPlotToolConfig=alice3_plotting.duplicationPlotToolConfig,
                trackQualityPlotToolConfig=alice3_plotting.trackQualityPlotToolConfig,
                trackSummaryPlotToolConfig=alice3_plotting.trackSummaryPlotToolConfig,
                filePath=str(outputDirRoot / f"performance_finding_{name}.root"),
                writeMatchingDetails=writeMatchingDetails,
            )
            s.addWriter(trackFinderPerfWriter)


def addTrackMerger(
        sequence: acts.examples.Sequencer,
        inputTrackCollections,
        outputTrackCollection: str,
        logLevel: acts.logging.Level = None,):

    customLogLevel = acts.examples.defaultLogging(sequence, logLevel)

    cfgMerger = TrackMergerAlgorithm.Config()
    cfgMerger.inputTrackCollections = inputTrackCollections
    cfgMerger.outputTrackCollection = outputTrackCollection

    trackMergerAlg = TrackMergerAlgorithm(
        config=cfgMerger,
        level=customLogLevel(),
    )

    sequence.addAlgorithm(trackMergerAlg)


def addIterativeTracking(
        s: acts.examples.Sequencer = None,
        cfg: ChainConfig.Config = ChainConfig.Config(),
        geo_dir: pathlib.Path = None,
        trackingGeometry: acts.TrackingGeometry = None,
        field: Literal[acts.MagneticFieldMapRz, acts.MagneticFieldMapXyz] = None,
        iterations: int = 0,
        inputTracks: str = "ckf_tracks",
        inputMeasurementSubset: str = "measurement_subset",
        outputDir: Optional[Union[Path, str]] = None,
        inputVertices: str = "",
        namePrefix: str = "",
        preVertexing: bool = False,):

    if iterations == 0:
        print("WARNING, addIterativeTracking(): No iterative tracking configured, skipping addIterativeTracking.")
        return

    seedTrackCollectionForMerging = []
    ckfTrackCollectionForMerging = []

    mergedSeedTrackCollection = f"{namePrefix}seed-tracks-merged"
    mergedTrackCollection = f"{namePrefix}ambi-tracks-merged"

    previousTracks = inputTracks
    previousMeasurementSubset = inputMeasurementSubset

    for iteration in range(0, iterations):
        # namePrefix namespaces the whole pass (used by the independent
        # pre-vertexing run, e.g. namePrefix="prevtx_"); empty for the main reco.
        prefix = namePrefix + (f"iter_{iteration}_" if iteration > 0 else "")

        iterOutputDir = Path(outputDir) / f"iter_{iteration}" if iterations > 1 else Path(outputDir)
        iterOutputDir.mkdir(parents=True, exist_ok=True)

        # Iteration 0 seeds off the initial (passed-in) measurement subset; later
        # iterations seed off the filtered leftover subset created below. (Using
        # inputMeasurementSubset explicitly rather than a name that only happens
        # to match when namePrefix is empty.)
        currentMeasurementSubset = (
            inputMeasurementSubset if iteration == 0
            else f"{prefix}measurement_subset"
        )
        currentSeedTracks = f"{prefix}seed-tracks"
        currentCkfTracks = f"{prefix}ckf_tracks"
        currentAmbiTracks = f"{prefix}ambi_tracks"

        # Each iteration of tracking uses left over measurements
        if iteration > 0:
            s.addAlgorithm(
                acts.examples.MeasurementFilterAlgorithm(
                    level=acts.logging.INFO,
                    inputTracks=previousTracks,
                    inputMeasurementSubset=previousMeasurementSubset,
                    outputMeasurementSubset=currentMeasurementSubset,
                    includeOutliers=False,
                )
            )

        seeding_cfg, geoSelectionConfigFile = alice3_seeding.get_seed_finder_config(
            cfg, iteration, preVertexing=preVertexing)

        # Per-iteration vertex-z seeding constraint (from the pre-vertexing pass).
        # Enabled only when a vertex collection is provided AND this iteration
        # opts in via iterationParams['useVertexZ']; k*sigma_z + margin are also
        # read per iteration (with sensible defaults).
        iterParams = cfg.seeding.iterationParams.get(iteration, {})
        useVtxZ = bool(inputVertices) and iterParams.get("useVertexZ", False)
        seedInputVertices = inputVertices if useVtxZ else ""
        seedVertexZNSigma = iterParams.get("vertexZNSigma", 3.0)
        seedVertexZMargin = iterParams.get("vertexZMargin", 0.0)

        # Per-iteration seeding algorithm. Defaults to GridTriplet; an entry
        # carrying 'seedingAlgorithm': 'Gbts' runs the ALICE 3 graph-based
        # seeder instead (see docs/GBTS.md). preVertexing always uses the
        # default algorithm, since its own params block has no override.
        seedingAlgorithmName = (
            "GridTriplet" if preVertexing
            else iterParams.get("seedingAlgorithm", "GridTriplet")
        )
        seedingAlgorithm = getattr(
            alice3_seeding.SeedingAlgorithm, seedingAlgorithmName
        )

        alice3_seeding.addSeeding(
            s,
            trackingGeometry,
            field,
            seedingAlgorithm=seedingAlgorithm,
            gbtsConfig=getattr(cfg, "gbts", None),
            gbtsOverrides=iterParams.get("gbts"),
            geoSelectionConfigFile=geo_dir / "../seedingConfigurations" / geoSelectionConfigFile,
            seedFinderConfigArg=seeding_cfg,
            seedFinderOptionsArg=alice3_seeding.DefaultSeedFinderOptionsArg,
            seedFilterConfigArg=alice3_seeding.testSeedFilterConfigArg,
            spacePointGridConfigArg=alice3_seeding.testSeedingGrid,
            seedingAlgorithmConfigArg=alice3_seeding.testSeedingAlgorithmConfigArg,

            # seedFinderOptionsArg=alice3_seeding.DefaultSeedFinderOptionsArg,
            # seedFilterConfigArg=alice3_seeding.DefaultSeedFilterConfigArg,
            # spacePointGridConfigArg=alice3_seeding.DefaultSpacePointGridConfigArg,
            # seedingAlgorithmConfigArg=alice3_seeding.DefaultSeedingAlgorithmConfigArg,
            outputDirRoot=iterOutputDir,
            # CKF start covariance. Its loc0/loc1 sigmas set the chi2 window on
            # the seed's own surface, i.e. how many measurements the CKF tests there.
            initialSigmas=getattr(cfg.seeding, "initialSigmas", None) or [
                1 * u.mm,
                1 * u.mm,
                1 * u.degree,
                1 * u.degree,
                0.1 * u.e / u.GeV,
                100 * u.ns,
            ],
            initialSigmaPtRel=0.1,
            initialVarInflation=alice3_seeding.PavelInitialVarInflation,
            particleHypothesis=acts.ParticleHypothesis.pion,
            inputMeasurements=currentMeasurementSubset,
            outputSpacePoints="spacepoints",
            prefix=prefix,
            writeMatchingDetails=cfg.seeding.writeMatchingDetails,
            writeSeedPerformance=cfg.seeding.writeSeedPerformance,
            inputVertices=seedInputVertices,
            vertexZNSigma=seedVertexZNSigma,
            vertexZMargin=seedVertexZMargin,
        )

        forceSeed = (not preVertexing) and cfg.tracking.iterationParams[iteration].get('forceSeedMeasurements', False)
        if forceSeed and "_VD" not in str(geoSelectionConfigFile):
            # forcing pays off for VD seeds; on ML seeds it keeps inconsistent (fake) seeds alive
            print(f"WARNING: iteration {iteration} forces the seed measurements but seeds on "
                  f"{geoSelectionConfigFile}, not on the VD (see perf_study/README.md section 9)")

        # run CKF
        addCKFTracks(
            s,
            trackingGeometry,
            field,
            get_track_selector_config(cfg=cfg, iteration=iteration),
            CkfConfig(
                seedDeduplication=cfg.tracking.seedDeduplication,
                stayOnSeed=cfg.tracking.stayOnSeed,
                chi2CutOffMeasurement=cfg.tracking.iterationParams[iteration]['ckfChi2Measurement'],
                chi2CutOffOutlier=cfg.tracking.iterationParams[iteration]['ckfChi2Outlier'],
                numMeasurementsCutOff=cfg.tracking.ckfMeasPerSurf,
            ),
            twoWay=cfg.tracking.twoWayCKF,
            outputDirRoot=iterOutputDir,
            writeTrackSummary=cfg.tracking.writeTrackSummary,
            writeTrackStates=False,
            writePerformance=cfg.tracking.writeCKFPerformance,
            writeMatchingDetails=cfg.tracking.writeMatchingDetailsCKF,
            inputMeasurements=currentMeasurementSubset,
            inputMeasurementParticlesMap="measurement_particles_map",
            inputParticleMeasurementsMap="particle_measurements_map",
            inputInitialTrackParameters=f"{prefix}estimatedparameters",
            inputSeeds=f"{prefix}estimatedseeds",
            outputTracks=currentCkfTracks,
            # IA: default in ACTS TrackTruthMatcher is 0.5, currently not configurable via reconstruction.py
            truthMatchingRatio=cfg.tracking.truthMatchingRatio,
            # the pre-vertexing pass seeds on the ML layers: never forced
            trackStateCreatorFactory=make_track_state_creator_factory(
                cfg, currentCkfTracks, forceSeedMeasurements=forceSeed),
            secondPassTrimmedCounts=getattr(cfg.tracking, "secondPassTrimmedCounts", False),
            logLevel=acts.logging.INFO,
        )

        s.addWhiteboardAlias(f"{prefix}tracks", currentCkfTracks)

        # Per-seed SeedAnalysisAlg tree (kinematics, quality, fake flags).
        # Switch on via cfg.seeding.runSeedAnalysis = True.
        if cfg.seeding.runSeedAnalysis:
            seedAnalysisCfg = SeedAnalysisAlg.Config()
            seedAnalysisCfg.inputTracks = currentSeedTracks
            seedAnalysisCfg.inputSeeds = f"{prefix}estimatedseeds"
            # Truth inputs for the two per-seed fake flags:
            #  - seed-track matching  -> isFake_seedTrack
            #  - measurement->particles map -> isFake_sharedParticle (3 SPs share a particle)
            seedAnalysisCfg.inputTrackParticleMatching = f"{prefix}seed_particle_matching"
            seedAnalysisCfg.inputMeasurementParticlesMap = "measurement_particles_map"
            seedAnalysisCfg.filePath = str(Path(outputDir) / f"{prefix}estimatedseeds_analysis.root")
            seedAnalysisAlg = SeedAnalysisAlg(seedAnalysisCfg, level=acts.logging.INFO)

            s.addWriter(seedAnalysisAlg)

        s = addAmbiguityResolution(
            s,
            cfg=cfg,
            ambiConfig=AmbiguityResolutionConfig(
                maximumSharedHits=cfg.tracking.maxSharedHits,
                nMeasurementsMin=cfg.tracking.iterationParams[iteration]['nMeasurementsMin'],
            ),
            prefix=prefix,
            outputDirRoot=iterOutputDir,
            writeTrackSummary=cfg.tracking.writeAmbiTrackSummary,
            writePerformance=cfg.tracking.writeAmbiPerformance,
            writeMatchingDetails=cfg.tracking.writeMatchingDetailsAmbi,
            logLevel=acts.logging.INFO,
            inputParticleMeasurementsMap="particle_measurements_map",
        )

        # Add the seed tracks for merging and the measurement mapping for this iteration
        seedTrackCollectionForMerging.append(currentSeedTracks)
        ckfTrackCollectionForMerging.append(currentCkfTracks)

        previousTracks = currentAmbiTracks
        previousMeasurementSubset = currentMeasurementSubset

        if cfg.monitoring.addUnbiasedResidualWriters:
            addResidualWriters(s,
                               name="AmbiTracksResiduals",
                               outputDirRoot=Path(outputDir),
                               tracks=currentAmbiTracks,
                               prefix=prefix)

    if iterations > 1:  # just one iteration - nothing to merge.

        # Merge seeds
        addTrackMerger(
            s,
            seedTrackCollectionForMerging,
            mergedSeedTrackCollection,
            acts.logging.INFO,
        )

        # Match seeds
        addTrackTruthMatcher(
            s,
            inputTracks=mergedSeedTrackCollection,
            inputParticles="particles_selected",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching="seed_merged_particle_matching",
            outputParticleTrackMatching="particle_seed_merged_matching",
            truthMatchingRatio=1.0,  # IA: as in reconstruction.py for all seeding algorithms
            logLevel=acts.logging.INFO,
        )

        # Merge CKF tracks from all iterations.
        mergedCkfTrackCollection = "ckf-tracks-merged"

        addTrackMerger(
            s,
            ckfTrackCollectionForMerging,
            mergedCkfTrackCollection,
            acts.logging.INFO,
        )

        # Run one final/global ambiguity resolution on the merged CKF tracks.
        from acts.examples import GreedyAmbiguityResolutionAlgorithm

        finalAmbiAlg = GreedyAmbiguityResolutionAlgorithm(
            level=acts.logging.INFO,
            inputTracks=mergedCkfTrackCollection,
            outputTracks=mergedTrackCollection,
            **acts.examples.defaultKWArgs(
                maximumSharedHits=cfg.tracking.maxSharedHits,
                nMeasurementsMin=cfg.tracking.nMeasurementsMinForFinalAmbiRes,
                maximumIterations=AmbiguityResolutionConfig().maximumIterations,
            ),
        )
        s.addAlgorithm(finalAmbiAlg)
        s.addWhiteboardAlias("tracks", mergedTrackCollection)

        # Match final merged ambiguity-resolved tracks.
        mergedOutputTrackParticleMatching = "ambi_tracks_merged_particle_matching"
        mergedOutputParticleTrackMatching = "particle_ambi_tracks_merged_matching"
        addTrackTruthMatcher(
            s,
            inputTracks=mergedTrackCollection,
            # matching only with particles_digitized_selected or all particles (also from Geant)
            inputParticles="particles_selected" if not cfg.tracking.doTruthMatchingAlsoWithParticlesFromGeant else "particles_simulated",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching=mergedOutputTrackParticleMatching,
            outputParticleTrackMatching=mergedOutputParticleTrackMatching,
            # IA: default in ACTS TrackTruthMatcher is 0.5, currently not configurable via reconstruction.py
            truthMatchingRatio=cfg.tracking.truthMatchingRatio,
            logLevel=acts.logging.INFO,
        )

        s.addWhiteboardAlias(
            "track_particle_matching", mergedOutputTrackParticleMatching  # IA: crucial for vertex truth matching!
        )
        s.addWhiteboardAlias(
            "particle_track_matching", mergedOutputParticleTrackMatching
        )

        # tracksummary_ambi-tracks-merged.root
        addTrackWriters(
            s,
            name=mergedTrackCollection,
            tracks=mergedTrackCollection,
            inputTrackParticleMatching=mergedOutputTrackParticleMatching,
            inputParticleTrackMatching=mergedOutputParticleTrackMatching,
            inputParticleMeasurementsMap="particle_measurements_map",
            outputDirCsv=None,
            outputDirRoot=outputDir,
            writeSummary=cfg.tracking.writeMergedTrackSummary,
            writeStates=cfg.tracking.writeTrackStates,
            writeFitterPerformance=False,
            writeFinderPerformance=False,
            writeCovMat=False,
            logLevel=acts.logging.INFO,
        )

        if cfg.tracking.writePatternRecoPerformance:
            s.addWriter(
                acts.examples.root.RootPatternRecognitionPerformanceWriter(
                    level=acts.logging.INFO,
                    inputTracks=mergedSeedTrackCollection,
                    inputParticles="particles_selected",
                    inputTrackParticleMatching="seed_merged_particle_matching",
                    inputParticleTrackMatching="particle_seed_merged_matching",
                    inputParticleMeasurementsMap="particle_measurements_map",
                    label="track",
                    effPlotToolConfig=alice3_plotting.effPlotToolConfig,
                    fakePlotToolConfig=alice3_plotting.fakePlotToolConfig,
                    duplicationPlotToolConfig=alice3_plotting.duplicationPlotToolConfig,
                    trackQualityPlotToolConfig=alice3_plotting.trackQualityPlotToolConfig,
                    trackSummaryPlotToolConfig=alice3_plotting.trackSummaryPlotToolConfig,
                    filePath=str(Path(outputDir) / "performance_merged_seed.root"),
                )
            )

            s.addWriter(
                acts.examples.root.RootPatternRecognitionPerformanceWriter(
                    level=acts.logging.INFO,
                    inputTracks=mergedTrackCollection,
                    inputParticles="particles_selected",
                    inputTrackParticleMatching=mergedOutputTrackParticleMatching,
                    inputParticleTrackMatching=mergedOutputParticleTrackMatching,
                    inputParticleMeasurementsMap="particle_measurements_map",
                    label="track",
                    effPlotToolConfig=alice3_plotting.effPlotToolConfig,
                    fakePlotToolConfig=alice3_plotting.fakePlotToolConfig,
                    duplicationPlotToolConfig=alice3_plotting.duplicationPlotToolConfig,
                    trackQualityPlotToolConfig=alice3_plotting.trackQualityPlotToolConfig,
                    trackSummaryPlotToolConfig=alice3_plotting.trackSummaryPlotToolConfig,
                    filePath=str(Path(outputDir) / "performance_merged_ambi_tracks.root"),
                )
            )

    return s


def addAmbiguityResolution(
        s,
        cfg: ChainConfig.Config = ChainConfig.Config(),
        ambiConfig: AmbiguityResolutionConfig = AmbiguityResolutionConfig(),
        tracks: Optional[str] = None,
        outputDirCsv: Optional[Union[Path, str]] = None,
        outputDirRoot: Optional[Union[Path, str]] = None,
        writeTrackSummary: bool = True,
        writeTrackStates: bool = False,
        writePerformance: bool = True,
        writeMatchingDetails: bool = False,
        writeCovMat=False,
        logLevel: Optional[acts.logging.Level] = None,
        inputParticleMeasurementsMap="particle_measurements_map",
        inputMeasurementParticlesMap="measurement_particles_map",
        prefix: str = "",
) -> None:

    from acts.examples import GreedyAmbiguityResolutionAlgorithm

    if tracks is None:
        tracks = f"{prefix}tracks"

    customLogLevel = acts.examples.defaultLogging(s, logLevel)
    outputAmbiTracks = f"{prefix}ambi_tracks"

    alg = GreedyAmbiguityResolutionAlgorithm(
        level=customLogLevel(),
        inputTracks=tracks,
        outputTracks=outputAmbiTracks,
        **acts.examples.defaultKWArgs(
            maximumSharedHits=ambiConfig.maximumSharedHits,
            nMeasurementsMin=ambiConfig.nMeasurementsMin,
            maximumIterations=ambiConfig.maximumIterations,
        ),
    )
    s.addAlgorithm(alg)
    s.addWhiteboardAlias(f"{prefix}tracks", alg.config.outputTracks)

    truthMatchCfg = TrackTruthMatcher.Config()
    truthMatchCfg.inputTracks = outputAmbiTracks
    truthMatchCfg.inputParticles = "particles_selected"
    truthMatchCfg.inputMeasurementParticlesMap = inputMeasurementParticlesMap
    truthMatchCfg.outputTrackParticleMatching = f"{outputAmbiTracks}_particle_matching"
    truthMatchCfg.outputParticleTrackMatching = f"particle_{outputAmbiTracks}_matching"
    truthMatchCfg.matchingRatio = cfg.tracking.truthMatchingRatio
    truthMatchCfg.doubleMatching = False
    truthMatchCfg.looperProtection = True
    truthMatchCfg.loop_absEta = 1.5
    truthMatchCfg.loop_maxPt = 0.8 * u.GeV
    truthMatchCfg.loop_maxParticleHits = 11

    matchAlg = TrackTruthMatcher(
        config=truthMatchCfg,
        level=customLogLevel(),
    )

    s.addAlgorithm(matchAlg)

    s.addWhiteboardAlias(
        "track_particle_matching", truthMatchCfg.outputTrackParticleMatching  # IA: crucial for vertex truth matching!
    )
    s.addWhiteboardAlias(
        "particle_track_matching", truthMatchCfg.outputParticleTrackMatching
    )
    addTrackWriters(
        s,
        name=outputAmbiTracks+"_writer",
        tracks=outputAmbiTracks,
        inputTrackParticleMatching=truthMatchCfg.outputTrackParticleMatching,
        inputParticleTrackMatching=truthMatchCfg.outputParticleTrackMatching,
        inputParticleMeasurementsMap=inputParticleMeasurementsMap,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeSummary=writeTrackSummary,
        writeStates=writeTrackStates,
        writeFitterPerformance=False,
        writeFinderPerformance=writePerformance,
        writeCovMat=writeCovMat,
        writeMatchingDetails=writeMatchingDetails,
        logLevel=logLevel,
    )

    return s


def addVertexPostProcessing(
    s,
    field: acts.MagneticFieldProvider,
    outputDir: Optional[Union[Path, str]] = None,
    writePerformance: bool = True,
):
    cfg = AliceVertexTruthMatcher.Config()
    cfg.inputVertices = "vertices"
    cfg.inputTracks = "tracks"
    cfg.inputParticles = "particles"
    cfg.inputTrackParticleMatching = "track_particle_matching"
    cfg.inputTruthVertices = "vertices_truth"
    cfg.outputVertexTruthMatching = "vertex_truth_matching_alice3"
    cfg.outputTruthVertexMatching = "truth_vertex_matching_alice3"

    cfg.vertexMatchThreshold = 0.999   # default is 0.7
    cfg.collapseSecondaryVertices = True   # fold HF/decay daughters into their primary
    # use time window (in fact, around nominal bunch crossings) for vertex truth matching
    cfg.useTimeWindow = True
    cfg.timeWindow = 120.0 * u.ns             # half-width; ~3–5x track time resolution
    cfg.printReco = False
    cfg.printTruth = False
    cfg.maxListPerVertex = -1

    s.addAlgorithm(AliceVertexTruthMatcher(cfg, acts.logging.INFO))

    if not writePerformance:
        return s

    # write  for vertex truth matching with vertices from many bunch crossings
    s.addWriter(
        acts.examples.root.RootVertexNTupleWriter(
            level=acts.logging.INFO,
            inputVertices="vertices",
            inputTracks="tracks",
            inputTruthVertices="vertices_truth",
            inputParticles="particles",
            inputSelectedParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            inputVertexTruthMatching="vertex_truth_matching_alice3",
            bField=field,
            writeTrackInfo=True,
            treeName="vertexing_alice3",
            filePath=str(outputDir / "performance_vertexing_alice3.root"),
        )
    )

    return s


def addGsfRefit(
        s: acts.examples.Sequencer,
        trackingGeometry,
        field,
        inputTracks: str = "tracks",
        outputTracks: str = "gsf_refit_tracks",
        inputParticles: str = "particles_selected",
        inputMeasurementParticlesMap: str = "measurement_particles_map",
        # ---- GSF options (defaults tuned for low-pT electrons) ----
        maxComponents: int = 12,
        weightCutoff: float = 1.0e-4,
        componentMergeMethod=None,        # default: ComponentMergeMethod.maxWeight
        mixtureReductionAlgorithm=None,   # default: MixtureReductionAlgorithm.KLDistance
        reverseFilteringCovarianceScaling: float = 100.0,
        initialVarInflation=None,         # default: 6 * [100.0]
        # ---- optional ROOT output of the refit tracks ----
        outputDirRoot: Optional[Union[Path, str]] = None,
        writeStates: bool = False,
        inputSimHits: str = "simhits",
        inputMeasurementSimHitsMap: str = "measurement_simhits_map",
        logLevel: Optional[acts.logging.Level] = None,
):
    """
    Refit an existing (ideally merged / ambiguity-resolved) track collection
    with the ACTS Gaussian-Sum Filter (GSF), which models Bremsstrahlung energy
    loss - useful for electron momentum reconstruction.

    Note: the GSF only re-fits the hits already on each input track (it does not
    re-find hits), so run it on the FINAL merged 'tracks' collection.

    Produces:
      <outputTracks>                            (default 'gsf_refit_tracks')
      <outputTracks>_particle_matching          truth matching of the refit tracks
      particle_<outputTracks>_matching
      and, if outputDirRoot is set, tracksummary_gsf_refit.root (+ trackstates).
    """
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if componentMergeMethod is None:
        componentMergeMethod = acts.examples.ComponentMergeMethod.maxWeight
    if mixtureReductionAlgorithm is None:
        mixtureReductionAlgorithm = acts.examples.MixtureReductionAlgorithm.KLDistance
    if initialVarInflation is None:
        initialVarInflation = 6 * [100.0]

    # Bethe-Heitler energy-loss approximation. clampToRange silences warnings for
    # detectors without a dedicated parametrisation (e.g. ALICE 3...)
    bha = acts.examples.AtlasBetheHeitlerApprox.makeDefault(clampToRange=True)

    gsf = acts.examples.makeGsfFitterFunction(
        trackingGeometry,
        field,
        betheHeitlerApprox=bha,
        maxComponents=maxComponents,
        weightCutoff=weightCutoff,
        componentMergeMethod=componentMergeMethod,
        mixtureReductionAlgorithm=mixtureReductionAlgorithm,
        reverseFilteringCovarianceScaling=reverseFilteringCovarianceScaling,
        level=customLogLevel(),
    )

    s.addAlgorithm(
        acts.examples.RefittingAlgorithm(
            customLogLevel(),
            inputTracks=inputTracks,
            outputTracks=outputTracks,
            initialVarInflation=initialVarInflation,
            fit=gsf,
        )
    )

    # Truth-match the refit tracks (reuse the ALICE 3 matcher wrapper above).
    addTrackTruthMatcher(
        s,
        inputTracks=outputTracks,
        inputParticles=inputParticles,
        inputMeasurementParticlesMap=inputMeasurementParticlesMap,
        outputTrackParticleMatching=f"{outputTracks}_particle_matching",
        outputParticleTrackMatching=f"particle_{outputTracks}_matching",
        logLevel=logLevel,
    )

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        outputDirRoot.mkdir(parents=True, exist_ok=True)

        s.addWriter(
            acts.examples.root.RootTrackSummaryWriter(
                level=customLogLevel(),
                inputTracks=outputTracks,
                inputParticles=inputParticles,
                inputTrackParticleMatching=f"{outputTracks}_particle_matching",
                filePath=str(outputDirRoot / "tracksummary_gsf_refit.root"),
            )
        )
        if writeStates:
            s.addWriter(
                acts.examples.root.RootTrackStatesWriter(
                    level=customLogLevel(),
                    inputTracks=outputTracks,
                    inputParticles=inputParticles,
                    inputTrackParticleMatching=f"{outputTracks}_particle_matching",
                    inputSimHits=inputSimHits,
                    inputMeasurementSimHitsMap=inputMeasurementSimHitsMap,
                    filePath=str(outputDirRoot / "trackstates_gsf_refit.root"),
                )
            )

    return s


def addGenAndGeantVertices(
        s: acts.examples.Sequencer,
        inputParticles: str = "particles_simulated",
        outputVertices: str = "vertices_gen_and_geant",
        outputDirRoot: Optional[Union[Path, str]] = None,
        fileName: str = "vertices_gen_and_geant.root",
        logLevel: Optional[acts.logging.Level] = None,
):
    """
    Build a vertex collection that includes BOTH generator-level and Geant4-level
    vertices (e.g. weak-decay vertices) from the simulated-particle collection,
    and (optionally) write it with the standard RootVertexWriter.

    inputParticles must be the FULL simulated set including the Geant4
    secondaries (e.g. "particles_simulated"), NOT a selected / primary-only
    collection, otherwise the secondary (decay) vertices are lost.
    """
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    cfg = SimVertexBuilder.Config()
    cfg.inputParticles = inputParticles
    cfg.outputVertices = outputVertices
    s.addAlgorithm(SimVertexBuilder(cfg, customLogLevel()))

    if outputDirRoot is not None:
        try:
            from acts.examples.root import RootVertexWriter
        except Exception:
            from acts.examples import RootVertexWriter
        outputDirRoot = Path(outputDirRoot)
        outputDirRoot.mkdir(parents=True, exist_ok=True)
        s.addWriter(
            RootVertexWriter(
                level=customLogLevel(),
                inputVertices=outputVertices,
                filePath=str(outputDirRoot / fileName),
            )
        )

    return s
