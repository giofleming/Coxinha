from pathlib import Path
from typing import Optional, Union, List, Tuple
from enum import Enum
from collections import namedtuple
import os
import array as arr

import acts
from acts import UnitConstants as u
import acts.examples
import acts.examples.json
import acts.examples.reconstruction as acts_reco

import alice3.performance.plotting as alice3_plotting

from acts.examples.reconstruction import (
    addSeeding,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedFilterConfigArg,
    SpacePointGridConfigArg,
    SeedingAlgorithmConfigArg,
    SeedingAlgorithm,
    # addSeedFilterML,
    # SeedFilterMLDBScanConfig,
)

from enum import Enum, auto


class SeedingMode(Enum):
    HIGH_PT = auto()
    LOW_PT = auto()
    LOOPERS = auto()
    STRANGENESS = auto()


#############################################
#### IGOR's SEEDING: SOME OTHER PARAMS   ####
#############################################

collisionRegion_forSeeds = 250  # 1000 # mm; large values - for V0 daughter reconstruction
seedingAlg = SeedingAlgorithm.GridTriplet
maxSeedsPerSpMConf = 2
maxQualitySeedsPerSpMConf = 2
enableMaterial = True
seedingGeoSelection_dir = Path(os.environ["MAIN_DIR"])
print(seedingGeoSelection_dir)

# Parameters to be passed
sigmaScattering = 5.0  # Tuned with particle gun, combine radLengthPerSeed = 0.05
radLengthPerSeed = 0.05
minSeedPt = 0.2  # to be changed per iteration
impParForSeeds = 1.0  # mm
bFieldInZ = 2.0  # T


##############################################
#### PAVEL's SEEDING: SOME OTHER PARAMS   ####
##############################################
seedParamOption = 0
deltaRmin = arr.array("f", [1., 1.0592])
deltaRmax = arr.array("f", [60., 61.587])
maxSeedsPerMiddleSp = arr.array("I", [1, 8])
maxSigmaScattering = arr.array("f", [5., 3.163])
radiationLengthPerSeed = arr.array("f", [0.1, 0.0785])
maxImpact = arr.array("f", [3., 22.29])
maxCotTheta = arr.array("f", [27.2899, 32.57108])
enableSeedconfirmation = False
initVarInflFactor = 1.0
PavelInitialVarInflation = [initVarInflFactor] * 6
enableMat = True


def get_seed_finder_config(cfg,
                           iteration: int = 0, preVertexing: bool = False, **overrides) -> Tuple[SeedFinderConfigArg, str]:
    """
    Get seeding configuration for a given iteration.
    Args:
    cfg: full Config object
    iteration: Iteration number (0,1,2)
    preVertexing: use the dedicated (lighter) pre-vertexing config instead of
                  the per-iteration one
    overriders: kwargs to override iteration-specific defaults
    """
    sc = cfg.seeding

    if preVertexing:
        it = dict(sc.preVertexingParams)  # dedicated lighter config
    else:
        if iteration not in sc.iterationParams:
            raise ValueError(f"No seeding iterationParams defined for iteration {iteration}")
        it = dict(sc.iterationParams[iteration])  # copy, don't mutate cfg
    seedingConfigFile = it.pop('seedingConfigFile')

    params = {
        "deltaR": tuple(v * u.mm for v in sc.deltaR),
        "collisionRegion": (-sc.collisionRegion * u.mm, sc.collisionRegion * u.mm),
        "cotThetaMax": sc.cotThetaMax,
        "seedConfirmation": sc.seedConfirmation,
        "centralSeedConfirmationRange": acts.SeedConfirmationRangeConfig(
            zMinSeedConf=sc.centralSeedConfirmation['zMinSeedConf'] * u.mm,
            zMaxSeedConf=sc.centralSeedConfirmation['zMaxSeedConf'] * u.mm,
            rMaxSeedConf=sc.centralSeedConfirmation['rMaxSeedConf'] * u.mm,
            nTopForLargeR=sc.centralSeedConfirmation['nTopForLargeR'],
            nTopForSmallR=sc.centralSeedConfirmation['nTopForSmallR'],
        ),
        "forwardSeedConfirmationRange": acts.SeedConfirmationRangeConfig(
            zMinSeedConf=sc.forwardSeedConfirmation['zMinSeedConf'] * u.mm,
            zMaxSeedConf=sc.forwardSeedConfirmation['zMaxSeedConf'] * u.mm,
            rMaxSeedConf=sc.forwardSeedConfirmation['rMaxSeedConf'] * u.mm,
            nTopForLargeR=sc.forwardSeedConfirmation['nTopForLargeR'],
            nTopForSmallR=sc.forwardSeedConfirmation['nTopForSmallR'],
        ),
        "useVariableMiddleSPRange": sc.useVariableMiddleSPRange,
        "deltaRMiddleSPRange": tuple(v * u.mm for v in sc.deltaRMiddleSPRange),
        "minPt": it['minPt'] * u.GeV,
        "impactMax": it['impactMax'] * u.mm,
        "z": tuple(v * u.mm if v is not None else None for v in it['z']),
        "r": tuple(v * u.mm if v is not None else None for v in it['r']),
        "maxSeedsPerSpM": it['maxSeedsPerSpM'],
        "sigmaScattering": it['sigmaScattering'],
        "radLengthPerSeed": it['radLengthPerSeed'],
        "deltaRTopSP": it['deltaRTopSP'],
        "deltaRBottomSP": it['deltaRBottomSP'],
        "rRangeMiddleSP": it['rRangeMiddleSP'],
        "useTimeDoubletCut": sc.useTimeInSeeding,
        "useTimeTripletCut": sc.useTimeInSeeding,
    }

    # Layer on any explicit call-site overrides (highest priority)
    params.update(overrides)

    return SeedFinderConfigArg(**params), seedingConfigFile


DefaultSeedFinderOptionsArg = SeedFinderOptionsArg(bFieldInZ=bFieldInZ * u.T, beamPos=(0 * u.mm, 0 * u.mm))

DefaultSeedFilterConfigArg = SeedFilterConfigArg(
    seedConfirmation=True if impParForSeeds < 4.0 else False,  # mm
    # If seedConfirmation is true we classify seeds as "high-quality" seeds.
    # Seeds that are not confirmed as "high-quality" are only selected if no
    # other "high-quality" seed has been found for that inner-middle doublet
    # Maximum number of normal seeds (not classified as "high-quality" seeds)
    # in seed confirmation
    maxSeedsPerSpMConf=maxSeedsPerSpMConf,  # 1 - USED_FOR_AUG_2025,#3, # CRUCIAL!!!!!!
    # 1 - USED_FOR_AUG_2025,   # Core/include/Acts/Seeding/SeedFilterConfig.hpp
    maxQualitySeedsPerSpMConf=maxQualitySeedsPerSpMConf,
    # Maximum number of "high-quality" seeds for each inner-middle SP-dublet in
    # seed confirmation. If the limit is reached we check if there is a normal
    # quality seed to be replaced
)


testSeedFilterConfigArg = SeedFilterConfigArg(
    # If seedConfirmation is true we classify seeds as "high-quality" seeds.
    # Seeds that are not confirmed as "high-quality" are only selected if no
    # other "high-quality" seed has been found for that inner-middle doublet
    seedConfirmation=True if impParForSeeds < 4.0 else False,  # mm
    # Maximum number of normal seeds (not classified as "high-quality" seeds)
    # in seed confirmation
    maxSeedsPerSpMConf=maxSeedsPerSpMConf,  # 1 - USED_FOR_AUG_2025,#3, # CRUCIAL!!!!!!
    # Maximum number of "high-quality" seeds for each inner-middle SP-dublet in
    # seed confirmation. If the limit is reached we check if there is a normal
    # quality seed to be replaced
    # 1 - USED_FOR_AUG_2025,   # Core/include/Acts/Seeding/SeedFilterConfig.hpp
    maxQualitySeedsPerSpMConf=maxQualitySeedsPerSpMConf,
)

testSeedingGrid = SpacePointGridConfigArg(
    zBinEdges=[
        -3000.,
        -1400.,
        -900.,
        -700.,
        -250.,
        250.,
        700.,
        900.,
        1400.,
        3000.,
    ],
    impactMax=1.0 * u.mm,  # not sure what it does
    phiBinDeflectionCoverage=3,  # 3
)

testSeedingAlgorithmConfigArg = SeedingAlgorithmConfigArg(

    zBinNeighborsTop=[
        [0, 0],
        [0, 0],  # The outerdisks are currently excluded
        [-1, 0],
        [-1, 0],  # Maybe -2, 0 could make sense if we want to include a middle space point in the second ML disk. We could also make the z bin such that each bin is a disk... We could also try [-1, 1] to be more inclusive on the barrel or [-2, 1]
        [-1, 1],
        [0, 1],
        [0, 1],
        [0, 0],
        [0, 0]
    ],
    zBinNeighborsBottom=[
        [0, 0],
        [0, 0],  # The outerdisks are currently excluded
        [0, 1],
        [0, 1],  # Maybe -2, 0 could make sense if we want to include a middle space point in the second ML disk. We could also make the z bin such that each bin is a disk... We could also try [-1, 1] to be more inclusive on the barrel or [-2, 1]
        [0, 0],
        [-1, 0],
        [-1, 0],
        [0, 0],
        [0, 0]
    ],
    # numPhiNeighbors=1,
)

DefaultSpacePointGridConfigArg = SpacePointGridConfigArg(
    # zBinEdges=[
    # -4000.0,
    # -2500.0,
    # -2000.0,
    # -1320.0,
    # -625.0,
    # -350.0,
    # -250.0,
    # 250.0,
    # 350.0,
    # 625.0,
    # 1320.0,
    # 2000.0,
    # 2500.0,
    # 4000.0,
    # ],
    impactMax=1.0 * u.mm,
    phiBinDeflectionCoverage=3,  # 3
)
DefaultSeedingAlgorithmConfigArg = SeedingAlgorithmConfigArg(
    # zBinNeighborsTop=[
    # [0, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 1],
    # [0, 1],
    # [0, 1],
    # [0, 1],
    # [0, 1],
    # [0, 1],
    # [0, 0],
    # ],
    # zBinNeighborsBottom=[
    # [0, 1],
    # [0, 1],
    # [0, 1],
    # [0, 1],
    # [0, 1],
    # [0, 1],
    # [0, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 0],
    # ],
    # numPhiNeighbors=1,
)


def addSpacePointsMaking(
        sequence: acts.examples.Sequencer,
        trackingGeometry: acts.TrackingGeometry,
        geoSelectionConfigFile: Union[Path, str],
        stripGeoSelectionConfigFile: Union[Path, str],
        inputMeasurements: str,
        outputSpacePoints: str,
        logLevel: acts.logging.Level = None,
):
    """adds space points making
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    # logLevel=acts.logging.DEBUG
    spAlg = acts.examples.SpacePointMaker(
        level=logLevel,
        inputMeasurements=inputMeasurements,
        outputSpacePoints=outputSpacePoints,
        trackingGeometry=trackingGeometry,
        geometrySelection=acts.examples.json.readJsonGeometryList(
            str(geoSelectionConfigFile)
        ),
        stripGeometrySelection=(
            acts.examples.json.readJsonGeometryList(str(stripGeoSelectionConfigFile))
            if stripGeoSelectionConfigFile
            else []
        ),
    )
    sequence.addAlgorithm(spAlg)
    return spAlg.config.outputSpacePoints


def _make_prefix(prefix: str = "", iterationIndex: int = 0) -> str:
    """
    ACTS v46.4.0 multi-pass convention.

    For prefix="iter_1_", addSeeding/addCKFTracks should use:
      iter_1_measurement_subset
      iter_1_spacepoints
      iter_1_seeds
      iter_1_estimatedparameters
      iter_1_estimatedseeds
      iter_1_seed-prototracks
      iter_1_seed-tracks

    iterationIndex is kept only for backward compatibility.
    """
    if prefix:
        return prefix

    if iterationIndex > 0:
        return f"iter_{iterationIndex}_"

    return ""


def addAlice3GbtsSeeding(
        sequence: acts.examples.Sequencer,
        spacePoints: str,
        trackingGeometry: acts.TrackingGeometry,
        gbtsCfg,
        overrides: Optional[dict] = None,
        outputSeeds: str = "seeds",
        bFieldInZ: float = 2.0 * u.T,
        logLevel: acts.logging.Level = None,
):
    """Schedule the ALICE 3 GBTS (graph-based) seeder.

    Parameters
    ----------
    spacePoints : name of the space point collection to seed from. Its geo
        selection must cover the same layers as the GBTS layer set, otherwise
        the surplus space points are simply dropped (logged as "not on a GBTS
        layer").
    gbtsCfg : the cfg.gbts block of ChainConfig
    overrides : per-iteration overrides, i.e. iterationParams[i]['gbts']
    See docs/GBTS.md for what each knob does.
    """
    from AliceActsPythonBindings import (
        Alice3GbtsSeedingAlgorithm,
        GbtsSeedFinderConfig,
        GbtsTrackingFilterConfig,
        GbtsLayerSpec,
    )

    overrides = dict(overrides or {})

    def get(name):
        return overrides.get(name, getattr(gbtsCfg, name))

    # The layer set is either given explicitly as (volume, layer) pairs, or
    # named among the presets on the gbts config block.
    layers = overrides.get("layers")
    if layers is None:
        layers = getattr(gbtsCfg, get("layerSet"))

    seedFinderConfig = GbtsSeedFinderConfig()
    seedFinderConfig.minPt = get("minPt") * u.GeV
    seedFinderConfig.minSeedLevel = get("minSeedLevel")
    seedFinderConfig.minDeltaRadius = get("minDeltaRadius") * u.mm
    seedFinderConfig.maxAbsTau = get("maxAbsTau")
    seedFinderConfig.maxOuterRadius = get("maxOuterRadius") * u.mm
    seedFinderConfig.nMaxPhiSlice = get("nMaxPhiSlice")
    seedFinderConfig.nMaxEdges = get("nMaxEdges")
    seedFinderConfig.hitShareThreshold = get("hitShareThreshold")
    seedFinderConfig.validateTriplets = get("validateTriplets")
    seedFinderConfig.matchBeforeCreate = get("matchBeforeCreate")
    seedFinderConfig.tauRatioPrecut = get("tauRatioPrecut") * ptScale
    seedFinderConfig.useAdaptiveCuts = get("useAdaptiveCuts")
    seedFinderConfig.addTriplets = get("addTriplets")
    seedFinderConfig.maxAbsEtaAddTriplets = get("maxAbsEtaAddTriplets")
    seedFinderConfig.edgeMaskMinEta = get("edgeMaskMinEta")
    seedFinderConfig.d0Max = get("d0Max") * u.mm
    seedFinderConfig.z0HistogramMaxBarrelOrder = get("z0HistogramMaxBarrelOrder")

    # ACTS scales its curvature and phi windows by tuningPt / minPt but leaves
    # the edge-to-edge linking tolerances absolute, so at low pt they, not the
    # doublet cuts, are what stops a graph from forming. Apply the same factor
    # here. See docs/GBTS.md.
    minPt = get("minPt")
    tuningPt = get("tuningPt")
    ptScale = (tuningPt / minPt) if get("scaleLinkingCutsWithPt") and minPt > 0 else 1.0
    seedFinderConfig.tuningPt = tuningPt * u.GeV
    seedFinderConfig.tauRatioCut = get("tauRatioCut") * ptScale
    seedFinderConfig.tauRatioCorr = get("tauRatioCorr") * ptScale
    seedFinderConfig.cutDPhiMax = get("cutDPhiMax") * ptScale
    seedFinderConfig.cutDCurvMax = get("cutDCurvMax") * ptScale / u.mm
    # No cluster widths or local positions are available from the ALICE 3
    # digitisation, so the machine-learning tau narrowing cannot run.
    seedFinderConfig.useClusterWidthCuts = False

    trackingFilterConfig = GbtsTrackingFilterConfig()
    trackingFilterConfig.sigmaX = get("filterSigmaX") * u.mm
    trackingFilterConfig.sigmaY = get("filterSigmaY") * u.mm
    trackingFilterConfig.maxDChi2X = get("filterMaxDChi2X")
    trackingFilterConfig.maxDChi2Y = get("filterMaxDChi2Y")
    trackingFilterConfig.addHit = get("filterAddHit")
    trackingFilterConfig.radLen = get("filterRadLen")
    # The filter's own curvature bound is absolute too, and its 1e-3/mm default
    # corresponds to ~300 MeV: every lower-pt chain is thrown away after the
    # first update unless it is opened up. |curvature| ~ 1 / (2 R), with
    # R[mm] = 1000 * pt[GeV] / (0.3 * B[T]).
    filterMaxCurvature = get("filterMaxCurvature")
    if filterMaxCurvature is None:
        bFieldTesla = bFieldInZ / u.T
        radiusMm = 1000.0 * minPt / (0.3 * bFieldTesla)
        filterMaxCurvature = get("filterMaxCurvatureSafety") / (2.0 * radiusMm)
    trackingFilterConfig.maxCurvature = filterMaxCurvature / u.mm
    # Keep the scattering model fixed while maxCurvature moves (see the
    # filterMsPerCurvature comment in ChainConfig).
    trackingFilterConfig.sigmaMS = get("filterMsPerCurvature") * filterMaxCurvature

    roiEta = get("roiEta")
    roiZ = get("roiZ")
    binTableZ0 = get("binTableZ0")

    # The binding macro does not install a kwargs constructor, so build the
    # Config object explicitly - the convention this package already uses for
    # its own algorithms.
    algConfig = Alice3GbtsSeedingAlgorithm.Config()
    algConfig.inputSpacePoints = spacePoints
    algConfig.outputSeeds = outputSeeds
    algConfig.trackingGeometry = trackingGeometry
    algConfig.layers = [GbtsLayerSpec(volume, layer) for volume, layer in layers]
    algConfig.layerConnections = list(overrides.get("layerConnections", []))
    algConfig.maxLayerSkip = get("maxLayerSkip")
    algConfig.connectionScanZ0Steps = get("connectionScanZ0Steps")
    algConfig.connectionScanEtaSteps = get("connectionScanEtaSteps")
    algConfig.etaBinWidth = get("etaBinWidth")
    algConfig.binTableZ0Min = binTableZ0[0] * u.mm
    algConfig.binTableZ0Max = binTableZ0[1] * u.mm
    algConfig.roiEtaMin = roiEta[0]
    algConfig.roiEtaMax = roiEta[1]
    algConfig.roiZMin = roiZ[0] * u.mm
    algConfig.roiZMax = roiZ[1] * u.mm
    algConfig.bFieldInZ = bFieldInZ
    algConfig.seedFinderConfig = seedFinderConfig
    algConfig.trackingFilterConfig = trackingFilterConfig

    alg = Alice3GbtsSeedingAlgorithm(
        algConfig, acts.examples.defaultLogging(sequence, logLevel)()
    )
    sequence.addAlgorithm(alg)
    return alg.config.outputSeeds


def addSeeding(
        s: acts.examples.Sequencer,
        trackingGeometry: acts.TrackingGeometry,
        field: acts.MagneticFieldProvider,
        geoSelectionConfigFile: Optional[Union[Path, str]] = None,
        stripGeoSelectionConfigFile: Optional[Union[Path, str]] = None,
        layerMappingConfigFile: Optional[Union[Path, str]] = None,
        connectorInputConfigFile: Optional[Union[Path, str]] = None,
        lutInputConfigFile: Optional[Union[Path, str]] = None,
        seedingAlgorithm: SeedingAlgorithm = SeedingAlgorithm.GridTriplet,
        trackSmearingSigmas: acts.examples.reconstruction.TrackSmearingSigmas = acts.examples.reconstruction.TrackSmearingSigmas(),
        initialSigmas: Optional[list] = None,
        initialSigmaQoverPt: Optional[float] = None,
        initialSigmaPtRel: Optional[float] = None,
        initialVarInflation: Optional[list] = None,
        seedFinderConfigArg: SeedFinderConfigArg = SeedFinderConfigArg(),
        seedFinderOptionsArg: SeedFinderOptionsArg = SeedFinderOptionsArg(),
        seedFilterConfigArg: SeedFilterConfigArg = SeedFilterConfigArg(),
        spacePointGridConfigArg: SpacePointGridConfigArg = SpacePointGridConfigArg(),
        seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg = SeedingAlgorithmConfigArg(),
        houghTransformConfig: acts.examples.HoughTransformSeeder.Config = acts.examples.HoughTransformSeeder.Config(),
        adaptiveHoughTransformConfig: Optional[
            acts.examples.AdaptiveHoughTransformSeeder.Config
        ] = None,
        hashingTrainingConfigArg: Optional[
            acts.examples.reconstruction.HashingTrainingConfigArg
        ] = acts.examples.reconstruction.HashingTrainingConfigArg(),
        hashingAlgorithmConfigArg: Optional[
            acts.examples.reconstruction.HashingAlgorithmConfigArg
        ] = acts.examples.reconstruction.HashingAlgorithmConfigArg(),
        particleHypothesis: Optional[
            acts.ParticleHypothesis
        ] = acts.ParticleHypothesis.pion,
        inputParticles: str = "particles",
        selectedParticles: str = "particles_selected",
        outputDirRoot: Optional[Union[Path, str]] = None,
        outputDirCsv: Optional[Union[Path, str]] = None,
        writeMatchingDetails: bool = False,
        logLevel: Optional[acts.logging.Level] = None,
        rnd: Optional[acts.examples.RandomNumbers] = None,
        inputMeasurements: str = "measurement_subset",
        outputSpacePoints: str = "spacepoints",
        outputSeeds: str = "seeds",
        iterationIndex: int = 0,
        prefix: str = "",
        # Optional vertex-z constraint on the seeding (GridTriplet only). When
        # inputVertices is a non-empty collection name, doublets whose z-origin
        # is not within k*sigma_z (+ margin) of a reconstructed vertex are cut.
        inputVertices: str = "",
        vertexZNSigma: float = 3.0,
        vertexZMargin: float = 0.0,
        writeSeedPerformance: bool = True,
        # GBTS only: the cfg.gbts block and this iteration's overrides
        # (iterationParams[i]['gbts']).
        gbtsConfig=None,
        gbtsOverrides: Optional[dict] = None,
) -> None:
    """This function steers the seeding
    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    geoSelectionConfigFile : Path|str, path, None
        Json file for space point geometry selection. Not required for SeedingAlgorithm.TruthSmeared.
    stripGeoSelectionConfigFile : Path|str, path, None
        Json file for space point geometry selection in strips. Needed for SpacePoint making.
    seedingAlgorithm : SeedingAlgorithm, Default
        seeding algorithm to use: one of Default (no truth information used), TruthSmeared, TruthEstimated
    trackSmearingSigmas : TrackSmearingSigmas(loc0, loc0PtA, loc0PtB, loc1, loc1PtA, loc1PtB, time, phi, theta, ptRel)
        TrackSmearing configuration.
        Defaults specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackParameterSmearing.hpp
    initialSigmas : list
        Sets the initial covariance matrix diagonal. This is ignored in case of TruthSmearing.
        Defaults specified in Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp
    initialVarInflation : list
        List of 6 scale factors to inflate the initial covariance matrix
        Defaults (all 1) specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackParameterSmearing.hpp
    seedFinderConfigArg : SeedFinderConfigArg(maxSeedsPerSpM, cotThetaMax, sigmaScattering, radLengthPerSeed, minPt, impactMax, deltaPhiMax, interactionPointCut, deltaZMax, maxPtScattering, zBinEdges, zBinsCustomLooping, rRangeMiddleSP, useVariableMiddleSPRange, binSizeR, seedConfirmation, centralSeedConfirmationRange, forwardSeedConfirmationRange, deltaR, deltaRBottomSP, deltaRTopSP, deltaRMiddleSPRange, collisionRegion, r, z)
        SeedFinderConfig settings. deltaR, deltaRBottomSP, deltaRTopSP, deltaRMiddleSPRange, collisionRegion, r, z.
        Defaults specified in Core/include/Acts/Seeding/SeedFinderConfig.hpp
    seedFinderOptionsArg :  SeedFinderOptionsArg(bFieldInZ, beamPos)
        Defaults specified in Core/include/Acts/Seeding/SeedFinderConfig.hpp
    seedFilterConfigArg : SeedFilterConfigArg(compatSeedWeight, compatSeedLimit, numSeedIncrement, seedWeightIncrement, seedConfirmation, maxSeedsPerSpMConf, maxQualitySeedsPerSpMConf, useDeltaRorTopRadius)
                                Defaults specified in Core/include/Acts/Seeding/SeedFilterConfig.hpp
    spacePointGridConfigArg : SpacePointGridConfigArg(rMax, zBinEdges, phiBinDeflectionCoverage, phi, maxPhiBins, impactMax)
                                SpacePointGridConfigArg settings. phi is specified as a tuple of (min,max).
        Defaults specified in Core/include/Acts/Seeding/SpacePointGrid.hpp
    seedingAlgorithmConfigArg : SeedingAlgorithmConfigArg(allowSeparateRMax, zBinNeighborsTop, zBinNeighborsBottom, numPhiNeighbors, useExtraCuts)
                                Defaults specified in Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/SeedingAlgorithm.hpp
    hashingTrainingConfigArg : HashingTrainingConfigArg(annoySeed, f)
                                Defaults specified in Plugins/Hashing/include/ActsPlugins/Hashing/HashingTrainingConfig.hpp
    hashingAlgorithmConfigArg : HashingAlgorithmConfigArg(bucketSize, zBins, phiBins)
                                Defaults specified in Plugins/Hashing/include/ActsPlugins/Hashing/HashingAlgorithmConfig.hpp
    particleHypothesis : Optional[acts.ParticleHypothesis]
        The hypothesis used for track finding. Defaults to pion.
    inputParticles : str, "particles"
        input particles name in the WhiteBoard
    selectedParticles : str, "particles_selected"
        selected particles name in the WhiteBoard
    outputDirRoot : Path|str, path, None
        the output folder for ROOT output, None triggers no output
    logLevel : acts.logging.Level, None
        logging level to override setting given in `s`
    rnd : RandomNumbers, None
        random number generator. Only used by SeedingAlgorithm.TruthSmeared.
    """

    collection_prefix = _make_prefix(prefix=prefix, iterationIndex=iterationIndex)

    outputSpacePointsWithPrefix = collection_prefix + outputSpacePoints
    outputSeedsWithPrefix = collection_prefix + outputSeeds

    outputEstimatedParameters = collection_prefix + "estimatedparameters"
    outputEstimatedSeeds = collection_prefix + "estimatedseeds"
    outputProtoTracks = collection_prefix + "seed-prototracks"
    outputSeedTracks = collection_prefix + "seed-tracks"

    outputSeedParticleMatching = collection_prefix + "seed_particle_matching"
    outputParticleSeedMatching = collection_prefix + "particle_seed_matching"

    trackFinderWriterOutName = collection_prefix + "performance_seeding.root"
    trackParamsWriterOutName = collection_prefix + "estimatedparameters.root"

    logLevel = acts.examples.defaultLogging(s, logLevel)()
    # logLevel = acts.logging.VERBOSE

    # Create starting parameters from either particle smearing or combined seed
    # finding and track parameters estimation
    if seedingAlgorithm == SeedingAlgorithm.TruthSmeared:
        acts_reco.addTruthSmearedSeeding(
            s=s,
            rnd=rnd,
            selectedParticles=selectedParticles,
            trackSmearingSigmas=trackSmearingSigmas,
            initialSigmas=initialSigmas,
            initialSigmaQoverPt=initialSigmaQoverPt,
            initialSigmaPtRel=initialSigmaPtRel,
            initialVarInflation=initialVarInflation,
            particleHypothesis=particleHypothesis,
            logLevel=logLevel,
        )
    else:

        spacePoints = addSpacePointsMaking(
            sequence=s,
            trackingGeometry=trackingGeometry,
            geoSelectionConfigFile=geoSelectionConfigFile,
            stripGeoSelectionConfigFile=stripGeoSelectionConfigFile,
            inputMeasurements=inputMeasurements,
            outputSpacePoints=outputSpacePointsWithPrefix,
            logLevel=logLevel,

        )
        seeds = None
        perSeedParticleHypothesis = None
        # proto tracks the estimation has to filter alongside the seeds
        unestimatedProtoTracks = None
        # Run either: truth track finding or seeding
        if seedingAlgorithm == SeedingAlgorithm.TruthEstimated:
            logger.info("Using truth track finding from space points for seeding")
            (
                seeds,
                perSeedParticleHypothesis,
                unestimatedProtoTracks,
            ) = addTruthEstimatedSeeding(
                s,
                spacePoints,
                selectedParticles,
                particleHypothesis=particleHypothesis,
                logLevel=logLevel,
            )
        elif seedingAlgorithm == SeedingAlgorithm.HoughTransform:
            houghTransformConfig.inputSpacePoints = [spacePoints]
            houghTransformConfig.inputMeasurements = inputMeasurements
            houghTransformConfig.outputProtoTracks = collection_prefix + "prototracks"
            houghTransformConfig.outputSeeds = outputSeedsWithPrefix
            houghTransformConfig.trackingGeometry = trackingGeometry
            seeds = acts_reco.addHoughTransformSeeding(s, houghTransformConfig, logLevel)

        elif seedingAlgorithm == SeedingAlgorithm.AdaptiveHoughTransform:
            if adaptiveHoughTransformConfig is None:
                adaptiveHoughTransformConfig = (
                    acts.examples.AdaptiveHoughTransformSeeder.Config()
                )

            adaptiveHoughTransformConfig.inputSpacePoints = [spacePoints]
            adaptiveHoughTransformConfig.outputProtoTracks = (
                collection_prefix + "prototracks"
            )
            adaptiveHoughTransformConfig.outputSeeds = outputSeedsWithPrefix
            adaptiveHoughTransformConfig.trackingGeometry = trackingGeometry
            adaptiveHoughTransformConfig.threshold = 4
            adaptiveHoughTransformConfig.noiseThreshold = 12
            adaptiveHoughTransformConfig.phiMinBinSize = 3.14 / (2.0 * 257.0)
            adaptiveHoughTransformConfig.qOverPtMinBinSize = 1.1 / (2.0 * 257.0)
            adaptiveHoughTransformConfig.qOverPtMin = 1.1
            adaptiveHoughTransformConfig.doSecondPhase = True
            adaptiveHoughTransformConfig.zMinBinSize = 1 * u.mm
            adaptiveHoughTransformConfig.cotThetaMinBinSize = 0.1
            adaptiveHoughTransformConfig.deduplicate = True
            seeds = acts_reco.addAdaptiveHoughTransformSeeding(
                s, adaptiveHoughTransformConfig, logLevel=logLevel
            )
        elif seedingAlgorithm == SeedingAlgorithm.Gbts:
            # ALICE 3 GBTS driver (AliceActsTrk::Alice3GbtsSeedingAlgorithm),
            # not the upstream ATLAS/ITk one: it builds the GBTS geometry from
            # the ALICE 3 tracking geometry instead of reading the ATLAS layer
            # mapping CSV and connector file. See docs/GBTS.md.
            if gbtsConfig is None:
                raise ValueError(
                    "SeedingAlgorithm.Gbts needs gbtsConfig (cfg.gbts)"
                )
            seeds = addAlice3GbtsSeeding(
                s,
                spacePoints,
                trackingGeometry,
                gbtsConfig,
                overrides=gbtsOverrides,
                outputSeeds=outputSeedsWithPrefix,
                bFieldInZ=seedFinderOptionsArg.bFieldInZ
                if seedFinderOptionsArg.bFieldInZ is not None
                else bFieldInZ * u.T,
                logLevel=logLevel,
            )
        elif seedingAlgorithm == SeedingAlgorithm.Hashing:
            seeds, buckets = acts_reco.addHashingSeeding(
                s,
                spacePoints,
                seedingAlgorithmConfigArg,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                spacePointGridConfigArg,
                hashingTrainingConfigArg,
                hashingAlgorithmConfigArg,
                logLevel,
            )
        elif seedingAlgorithm == SeedingAlgorithm.GridTriplet:
            seeds = acts.examples.reconstruction.addGridTripletSeeding(
                s,
                spacePoints,
                seedingAlgorithmConfigArg,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                spacePointGridConfigArg,
                logLevel,
                outputSeeds=outputSeedsWithPrefix,
                inputVertices=inputVertices,
                vertexZNSigma=vertexZNSigma,
                vertexZMargin=vertexZMargin,
            )

        elif seedingAlgorithm == SeedingAlgorithm.OrthogonalTriplet:
            seeds = acts_reco.addOrthogonalTripletSeeding(
                s,
                spacePoints,
                seedingAlgorithmConfigArg,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                spacePointGridConfigArg,
                logLevel,
            )
        else:
            pass

        parEstimateAlg = acts.examples.TrackParamsEstimationAlgorithm(
            level=logLevel,
            inputSeeds=seeds,
            inputParticleHypotheses=perSeedParticleHypothesis,
            outputTrackParameters=outputEstimatedParameters,
            outputSeeds=outputEstimatedSeeds,
            trackingGeometry=trackingGeometry,
            magneticField=field,
            **acts.examples.defaultKWArgs(
                initialSigmas=initialSigmas,
                initialSigmaQoverPt=initialSigmaQoverPt,
                initialSigmaPtRel=initialSigmaPtRel,
                initialVarInflation=initialVarInflation,
                particleHypothesis=particleHypothesis,
            ),
        )
        s.addAlgorithm(parEstimateAlg)

        s.addAlgorithm(
            acts.examples.SeedsToProtoTracks(
                level=logLevel,
                inputSeeds=outputEstimatedSeeds,
                outputProtoTracks=outputProtoTracks,
            )
        )

        s.addAlgorithm(
            acts.examples.ProtoTracksToTracks(
                level=logLevel,
                inputProtoTracks=outputProtoTracks,
                inputTrackParameters=outputEstimatedParameters,
                inputMeasurements=inputMeasurements,
                outputTracks=outputSeedTracks,
            )
        )

        print("PF::DEBUG collection_prefix", collection_prefix)
        print("PF::DEBUG inputMeasurements", inputMeasurements)
        print("PF::DEBUG measurement_particles_map")
        print("PF::DEBUG", outputSeedParticleMatching)
        print("PF::DEBUG", outputParticleSeedMatching)

        inputMeasurementParticlesMap = "measurement_particles_map"
        inputParticleMeasurementMap = "particle_measurements_map"

        s.addAlgorithm(
            acts.examples.TrackTruthMatcher(
                level=logLevel,
                inputTracks=outputSeedTracks,
                inputParticles=selectedParticles,
                inputMeasurementParticlesMap=inputMeasurementParticlesMap,
                outputTrackParticleMatching=outputSeedParticleMatching,
                outputParticleTrackMatching=outputParticleSeedMatching,
                matchingRatio=1.0,
                doubleMatching=False,
            )
        )

        if outputDirRoot is not None and writeSeedPerformance:
            addSeedPerformanceWriters(
                s,
                outputDirRoot,
                outputSeedTracks,
                outputProtoTracks,
                selectedParticles,
                inputParticles,
                inputParticleMeasurementMap,
                outputSeedParticleMatching,
                outputParticleSeedMatching,
                parEstimateAlg.config.outputTrackParameters,
                trackFinderWriterOutName,
                trackParamsWriterOutName,
                writeMatchingDetails,
                logLevel,
            )

    return s


# This adds only the performance algorithms to check the performance on a set of seed-tracks
# - Schedule the TruthMatcher
# - Schedule the TrackFinderPerformanceWriter for plots

def addSeedPerformanceAlgorithms(
        sequence: acts.examples.Sequencer,
        outputDirRoot: Union[Path, str],
        tracks: str,
        inputParticles: str,
        inputMeasurementParticlesMap: str,
        outputTrackParticleMatching: str,
        outputParticleTrackMatching: str,
        logLevel: acts.logging.Level = None,
):

    customLogLevel = acts.examples.defaultLogging(sequence, logLevel)

    sequence.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=customLogLevel(),
            inputTracks=tracks,
            inputParticles=inputParticles,
            inputMeasurementParticlesMap=inputMeasurementParticlesMap,
            outputTrackParticleMatching=outputTrackParticleMatching,
            outputParticleTrackMatching=outputParticleTrackMatching,
            matchingRatio=1.0,
            doubleMatching=False,
        )
    )


def addSeedPerformanceWriters(
        sequence: acts.examples.Sequencer,
        outputDirRoot: Union[Path, str],
        tracks: str,
        prototracks: str,
        selectedParticles: str,
        inputParticles: str,
        inputParticleMeasurementsMap: str,
        inputTrackParticleMatching: str,
        inputParticleTrackMatching: str,
        outputTrackParameters: str,
        trackFinderWriterOutName: str,
        trackParamsWriterOutName: str,
        writeMatchingDetails: bool = False,
        logLevel: acts.logging.Level = None,
):
    """Writes seeding related performance output"""
    customLogLevel = acts.examples.defaultLogging(sequence, logLevel)
    outputDirRoot = Path(outputDirRoot)
    if not outputDirRoot.exists():
        outputDirRoot.mkdir(parents=True, exist_ok=True)

    sequence.addWriter(
        acts.examples.root.RootPatternRecognitionPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks=tracks,
            inputParticles=selectedParticles,
            inputTrackParticleMatching=inputTrackParticleMatching,
            inputParticleTrackMatching=inputParticleTrackMatching,
            inputParticleMeasurementsMap=inputParticleMeasurementsMap,
            # label="seed",
            effPlotToolConfig=alice3_plotting.effPlotToolConfig,
            fakePlotToolConfig=alice3_plotting.fakePlotToolConfig,
            writeMatchingDetails=writeMatchingDetails,
            filePath=str(outputDirRoot / trackFinderWriterOutName),
        )
    )

    sequence.addWriter(
        acts.examples.root.RootTrackParameterWriter(
            level=customLogLevel(),
            inputTrackParameters=outputTrackParameters,
            inputProtoTracks=prototracks,
            inputParticles=inputParticles,
            inputSimHits="simhits",
            inputMeasurementParticlesMap="measurement_particles_map",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDirRoot / trackParamsWriterOutName),
            treeName="estimatedparams",
        )
    )
