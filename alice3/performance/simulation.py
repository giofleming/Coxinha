from typing import Optional, Union, List, Literal, Tuple, Any
from pathlib import Path

import acts
import acts.examples
import acts.examples.reconstruction as acts_reco
from acts import UnitConstants as u

from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addPythia8,
    addGeant4,
    addGenParticleSelection,
    addSimParticleSelection,
    ParticleSelectorConfig,
)

from acts.examples.root import (
    RootParticleWriter,
    RootSimHitWriter,
)

from AliceActsPythonBindings import PerPixelNoiseAlgorithm

# Load config
from job_configs import ChainConfig
cfg = ChainConfig.Config()


def addDigiParticleSelection(
        s: acts.examples.Sequencer,
        config: ParticleSelectorConfig,
        measurementLayers: List[Tuple[int, int]] = [],
        logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """
    This function steers the particle selection after digitization.

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the ParticleSelector
    config: ParticleSelectorConfig
        the particle selection configuration
    measurementLayers: list of (VolID,lyID) to where adding the measurement counter. 
    at least 1 measurement will be requested
    """
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # Get filtered kwargs
    kwargs = acts.examples.defaultKWArgs(**acts.examples.simulation._getParticleSelectionKWargs(config))

    # Create config object
    cfg = acts.examples.ParticleSelector.Config()

    # Apply all kwargs to config
    for key, value in kwargs.items():
        setattr(cfg, key, value)

    # Override specific fields
    cfg.inputParticles = "particles_simulated_selected"
    cfg.inputParticleMeasurementsMap = "particle_measurements_map"
    cfg.inputMeasurements = "measurements"
    cfg.outputParticles = "tmp_particles_digitized_selected"

    geoIDs = []
    for tup in measurementLayers:
        geoIDs.append(acts.GeometryIdentifier(volume=tup[0], layer=tup[1]))

    if len(measurementLayers) > 0:
        cfg.measurementCounter.addCounter(geoIDs, 1, 1)

    selector = acts.examples.ParticleSelector(
        config=cfg,
        level=customLogLevel(),
    )

    s.addAlgorithm(selector)

    s.addWhiteboardAlias("particles_selected", selector.config.outputParticles)
    s.addWhiteboardAlias(
        "particles_digitized_selected", selector.config.outputParticles
    )


# algorithm to add per-pixel noise to the digitized measurements (June 2026)
def addPerPixelNoise(
    s,
    trackingGeometry,
    digiConfigFile,
    layersByVolume=None,          # {volumeId: [layerIds]}; [] => all layers; None/{} => no-op
    noiseProbPerPixel=1e-7,       # fake-hit probability per pixel per event
    noiseActivation=1.0,
    rnd=None,
    mergeCommonCorner=True,
    mergeNsigma=1.0,
    # pure-noise time window in ACTS native units (pass e.g. X * u.ns); 0 => t=0; only used if digi smears time:
    noiseTimeWindow=0.0,
    inputSuffix="",          # the digitizer must emit "<name><inputSuffix>"
    outputSuffix="_with_noise",          # the digitizer must emit "<name><outputSuffix>"
    logLevel=None,
    outputDirRoot: Optional[Union[Path, str]] = None,
):
    customLogLevel = acts.examples.defaultLogging(s, logLevel)
    rnd = rnd or acts.examples.RandomNumbers()
    layersByVolume = layersByVolume or {}

    cfg = PerPixelNoiseAlgorithm.Config()
    cfg.trackingGeometry = trackingGeometry
    cfg.digitizationConfigs = acts.examples.json.readDigiConfigFromJson(str(digiConfigFile))
    cfg.randomNumbers = rnd
    cfg.noiseProbPerPixel = noiseProbPerPixel
    cfg.noiseActivation = noiseActivation
    cfg.layersByVolume = {int(k): [int(x) for x in v] for k, v in layersByVolume.items()}
    cfg.mergeCommonCorner = mergeCommonCorner
    cfg.mergeNsigma = mergeNsigma
    # noiseTimeWindow is already in ACTS native time units (c*t); callers pass a
    # ns value multiplied by u.ns (see alice3_full_chain.py). Only affects
    # surfaces whose digi config smears time.
    cfg.noiseTimeWindow = noiseTimeWindow

    # Inputs: the renamed digitizer outputs.
    cfg.inputMeasurements = "measurements" + inputSuffix
    cfg.inputClusters = "clusters" + inputSuffix
    cfg.inputMeasurementParticlesMap = "measurement_particles_map" + inputSuffix
    cfg.inputMeasurementSimHitsMap = "measurement_simhits_map" + inputSuffix

    # Outputs: the canonical names consumed by seeding / CKF / writers.
    cfg.outputMeasurements = "measurements" + outputSuffix
    cfg.outputMeasurementSubset = "measurement_subset" + outputSuffix
    cfg.outputClusters = "clusters" + outputSuffix
    cfg.outputMeasurementParticlesMap = "measurement_particles_map" + outputSuffix
    cfg.outputMeasurementSimHitsMap = "measurement_simhits_map" + outputSuffix
    cfg.outputParticleMeasurementsMap = "particle_measurements_map" + outputSuffix
    cfg.outputSimHitMeasurementsMap = "simhit_measurements_map" + outputSuffix

    s.addAlgorithm(PerPixelNoiseAlgorithm(config=cfg, level=customLogLevel()))

    # create aliases to redirect to objects with noise
    s.addWhiteboardAlias("measurements", "measurements" + outputSuffix)
    s.addWhiteboardAlias("measurement_subset", "measurement_subset" + outputSuffix)
    s.addWhiteboardAlias("clusters", "clusters" + outputSuffix)
    s.addWhiteboardAlias("measurement_particles_map", "measurement_particles_map" + outputSuffix)
    s.addWhiteboardAlias("measurement_simhits_map", "measurement_simhits_map" + outputSuffix)
    s.addWhiteboardAlias("particle_measurements_map", "particle_measurements_map" + outputSuffix)
    s.addWhiteboardAlias("simhit_measurements_map", "simhit_measurements_map" + outputSuffix)

    if outputDirRoot is None:
        return s

    rmwCfg = acts.examples.root.RootMeasurementWriter.Config(
        inputMeasurements="measurements",               # canonical = noise-augmented
        inputSimHits="simhits",
        inputMeasurementSimHitsMap="measurement_simhits_map",
        filePath=str(Path(outputDirRoot) / f"measurements{outputSuffix}.root"),   # don't clobber measurements_raw.root
        surfaceByIdentifier=trackingGeometry.geoIdSurfaceMap(),
    )
    s.addWriter(acts.examples.root.RootMeasurementWriter(rmwCfg, acts.logging.INFO))

    return s


# Wrapper around acts addGeant4 to provide better control on what output to store
def addAlice3Geant4(
        s: acts.examples.Sequencer,
        detector: Optional[Any],
        trackingGeometry: acts.TrackingGeometry,
        field: acts.MagneticFieldProvider,
        rnd: acts.examples.RandomNumbers,
        volumeMappings: List[str] = [],
        materialMappings: List[str] = ["Silicon"],
        inputParticles: str = "particles_generated_selected",
        outputParticles: str = "particles_simulated",
        outputSimHits: str = "simhits",
        recordHitsOfSecondaries=True,
        recordPropagationSummaries=False,
        keepParticlesWithoutHits=True,
        writeHelixParameters: bool = False,
        outputDirCsv: Optional[Union[Path, str]] = None,
        outputDirRoot: Optional[Union[Path, str]] = None,
        outputDirObj: Optional[Union[Path, str]] = None,
        logLevel: Optional[acts.logging.Level] = None,
        killVolume: Optional[acts.Volume] = None,
        killAfterTime: float = float("inf"),
        killSecondaries: bool = False,
        physicsList: str = "FTFP_BERT",
        detectorConstructionOptions=None,
        propagatorLargestAcceptableStep: float = float("inf"),
        writeParticles=False,
        writeHits=True,

) -> None:

    # Wrap around acts geant4
    addGeant4(
        s,
        detector=detector,
        trackingGeometry=trackingGeometry,
        field=field,
        volumeMappings=volumeMappings,
        materialMappings=materialMappings,
        outputDirRoot=None,  # avoid calling the dumpers
        rnd=rnd,
        logLevel=logLevel,
        killVolume=killVolume,
        killAfterTime=killAfterTime,  # 25 * u.ns,
        recordHitsOfSecondaries=recordHitsOfSecondaries,
        keepParticlesWithoutHits=keepParticlesWithoutHits,
        killSecondaries=killSecondaries,
        propagatorLargestAcceptableStep=propagatorLargestAcceptableStep,
    )

    outputDirRoot = Path(outputDirRoot)
    if not outputDirRoot.exists():
        outputDirRoot.mkdir()

    if writeParticles:

        s.addWriter(
            RootParticleWriter(
                level=logLevel,
                inputParticles="particles_simulated",
                bField=field,
                writeHelixParameters=writeHelixParameters,
                filePath=str(outputDirRoot / "particles_simulation.root"),
            )
        )

    if writeHits:

        s.addWriter(
            RootSimHitWriter(
                level=logLevel,
                inputSimHits="simhits",
                filePath=str(outputDirRoot / "hits.root"),
            )
        )
