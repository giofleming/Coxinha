from typing import Optional, Union, List, Literal
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
    addGenParticleSelection,
    addSimParticleSelection,
    ParticleSelectorConfig,
)

from typing import Iterable

from acts.examples import (
    FixedMultiplicityGenerator,
    EventGenerator,
    RandomNumbers,
)

from acts.examples.root import (
    RootParticleWriter,
    RootVertexWriter
)


from AliceActsPythonBindings import BunchCrossingVertexGenerator

# Load config
from job_configs import ChainConfig


def addAlice3ParticleGun(s: acts.examples.Sequencer = None,
                         rnd: acts.examples.RandomNumbers = None,
                         cfg: ChainConfig.Config = ChainConfig.Config(),
                         outputDir: Optional[Union[Path, str]] = None):

    addParticleGun(
        s=s,
        momentumConfig=MomentumConfig(
            cfg.particleGun.gunPtRange[0] * u.GeV, cfg.particleGun.gunPtRange[1] * u.GeV, transverse=True
        ),
        etaConfig=EtaConfig(cfg.particleGun.gunEtaRange[0], cfg.particleGun.gunEtaRange[1], uniform=True),
        phiConfig=PhiConfig(cfg.particleGun.gunPhiRange[0], cfg.particleGun.gunPhiRange[1]),
        particleConfig=ParticleConfig(
            cfg.particleGun.gunMult,
            acts.PdgParticle(cfg.particleGun.gunPID),
            randomizeCharge=cfg.particleGun.gunRandCharge,
        ),
        vtxGen=acts.examples.GaussianVertexGenerator(
            stddev=acts.Vector4(
                0.0001 * u.mm, 0.0001 * u.mm, cfg.particleGun.sigmaVz * u.mm, 0.0001 * u.ns
            ),
            mean=acts.Vector4(0, 0, cfg.particleGun.meanVz * u.mm, 0),
        ),
        rnd=rnd,
        logLevel=acts.logging.INFO,
        outputDirRoot=outputDir,
    )


def addAlice3Pythia8(s: acts.examples.Sequencer = None,
                     rnd: acts.examples.RandomNumbers = None,
                     cfg: ChainConfig.Config = ChainConfig.Config(),
                     outputDir: Optional[Union[Path, str]] = None):

    beam = acts.PdgParticle.eProton
    hardProcess = []
    cmsEnergy = 0.
    if cfg.pythia.system == "pp":
        beam = acts.PdgParticle.eProton,
        cmsEnergy = 13.6 * acts.UnitConstants.TeV
        hardProcess = [
            "SoftQCD:inelastic = on",
            "Tune:pp = 14",
            "ParticleDecays:limitTau0 = on",
            "ParticleDecays:tau0Max = 10",
        ]  # tune 14 is Monash
        pileupProcess = [
            "SoftQCD:inelastic = on",
            "Tune:pp = 14",
            "ParticleDecays:limitTau0 = on",
            "ParticleDecays:tau0Max = 10",
        ]  # tune 14 is Monash
    elif cfg.pythia.system == "PbPbMB":
        beam = acts.PdgParticle.eLead,
        cmsEnergy = 5.36 * acts.UnitConstants.TeV
        hardProcess = [
            "SoftQCD:inelastic = on",
            "HeavyIon:SigFitErr =  0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0",
            "HeavyIon:SigFitDefPar = 17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0",
            "HeavyIon:SigFitNGen = 20",
            "ParticleDecays:limitTau0 = on",
            "ParticleDecays:tau0Max = 10",
        ]
        pileupProcess = [  # MB as pileup process
            "SoftQCD:inelastic = on",
            "HeavyIon:SigFitErr =  0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0",
            "HeavyIon:SigFitDefPar = 17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0",
            "HeavyIon:SigFitNGen = 20",
            "ParticleDecays:limitTau0 = on",
            "ParticleDecays:tau0Max = 10",
        ]
    elif cfg.pythia.system == "PbPbCentral":
        beam = acts.PdgParticle.eLead,
        cmsEnergy = 5.36 * acts.UnitConstants.TeV
        hardProcess = [
            "SoftQCD:inelastic = on",
            "HeavyIon:bWidth=0.1",
            "HeavyIon:SigFitErr =  0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0",
            "HeavyIon:SigFitDefPar = 17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0",
            "HeavyIon:SigFitNGen = 20",
            "ParticleDecays:limitTau0 = on",
            "ParticleDecays:tau0Max = 10",
        ]
        pileupProcess = [  # MB as pileup process
            "SoftQCD:inelastic = on",
            "HeavyIon:SigFitErr =  0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0",
            "HeavyIon:SigFitDefPar = 17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0",
            "HeavyIon:SigFitNGen = 20",
            "ParticleDecays:limitTau0 = on",
            "ParticleDecays:tau0Max = 10",
        ]
    else:
        print("Pythia process ", cfg.pythia.system, " not available!")

    # primary vertex generator:
    if cfg.pythia.useDiscreteBCtime:  # special vertex generator for discrete bunch-crossing times
        from AliceActsPythonBindings import BunchCrossingVertexGenerator
        vtxGen = BunchCrossingVertexGenerator()
        vtxGen.stddev = acts.Vector4(
            0.000125 * u.mm,
            0.000125 * u.mm,
            cfg.pythia.sigmaVz * u.mm,
            0,  # time component ignored by this generator
        )
        vtxGen.mean = acts.Vector4(0, 0, cfg.pythia.meanVz * u.mm, 0)
        vtxGen.nBunchCrossings = cfg.pythia.nBunchCrossings
        vtxGen.bunchSpacing = cfg.pythia.bunchSpacing * u.ns
    else:  # default Gaussian vertex generator
        vtxGen = acts.examples.GaussianVertexGenerator(
            stddev=acts.Vector4(
                0.000125 * u.mm,
                0.000125 * u.mm,
                cfg.pythia.sigmaVz * u.mm,
                0.0001 * u.ns,  # !! makes 0 sense
            ),
            mean=acts.Vector4(0, 0, cfg.pythia.meanVz * u.mm, 0),
        )

    s = addPythia8(
        s,
        npileup=cfg.pythia.pileup,
        beam=beam[0],
        cmsEnergy=cmsEnergy,
        hardProcess=hardProcess,
        pileupProcess=pileupProcess,
        vtxGen=vtxGen,
        rnd=rnd,
        logLevel=acts.logging.INFO,
        outputDirRoot=outputDir)

    return s


def addAlice3ParticleGunWithPileup(s: acts.examples.Sequencer = None,
                                   rnd: acts.examples.RandomNumbers = None,
                                   cfg: ChainConfig.Config = ChainConfig.Config(),
                                   outputDir: Optional[Union[Path, str]] = None) -> None:
    """This function steers the particle generation using the particle gun

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the particle gun steps (returned from addParticleGun)
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    momentumConfig : MomentumConfig(min, max, transverse, logUniform)
        momentum configuration: minimum momentum, maximum momentum, transverse, log-uniform
    etaConfig : EtaConfig(min, max, uniform)
        pseudorapidity configuration: eta min, eta max, uniform
    phiConfig : PhiConfig(min, max)
        azimuthal angle configuration: phi min, phi max
    particleConfig : ParticleConfig(num, pdg, randomizeCharge, charge, mass)
        particle configuration: number of particles, particle type, charge flip
    multiplicity : int, 1
        number of generated vertices
    vtxGen : VertexGenerator, None
        vertex generator module
    printParticles : bool, False
        print generated particles
    rnd : RandomNumbers, None
        random number generator
    """

    momentumConfig = MomentumConfig(cfg.particleGun.gunPtRange[0] * u.GeV,
                                    cfg.particleGun.gunPtRange[1] * u.GeV,
                                    transverse=True)
    etaConfig = EtaConfig(cfg.particleGun.gunEtaRange[0], cfg.particleGun.gunEtaRange[1], uniform=True)
    phiConfig = PhiConfig(cfg.particleGun.gunPhiRange[0], cfg.particleGun.gunPhiRange[1])
    particleConfig = ParticleConfig(cfg.particleGun.gunMult,
                                    acts.PdgParticle(cfg.particleGun.gunPID),
                                    randomizeCharge=cfg.particleGun.gunRandCharge)
    vtxGen = acts.examples.GaussianVertexGenerator(stddev=acts.Vector4(0.0001 * u.mm, 0.0001 * u.mm,
                                                                       cfg.particleGun.sigmaVz * u.mm, 0.0001 * u.ns),
                                                   mean=acts.Vector4(0, 0, cfg.particleGun.meanVz * u.mm, 0))
    logLevel = acts.logging.INFO
    outputDirRoot = outputDir
    multiplicity = 1
    printParticles = False,
    npileup = 1
    beam = [2212, 2212]  # default: acts.PdgParticle.eProton
    cmsEnergy = 14 * acts.UnitConstants.TeV
    pileupProcess = ["SoftQCD:all = on"]

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    rnd = rnd or RandomNumbers(seed=228)

    generators = []

    generators.append(EventGenerator.Generator(
        multiplicity=FixedMultiplicityGenerator(n=multiplicity),
        vertex=vtxGen
        or acts.examples.GaussianVertexGenerator(mean=acts.Vector4(0, 0, 0, 0),
                                                 stddev=acts.Vector4(0, 0, 0, 0)),
        particles=acts.examples.ParametricParticleGenerator(
            **acts.examples.defaultKWArgs(
                p=(momentumConfig.min, momentumConfig.max),
                pTransverse=momentumConfig.transverse,
                pLogUniform=momentumConfig.logUniform,
                eta=(etaConfig.min, etaConfig.max),
                phi=(phiConfig.min, phiConfig.max),
                etaUniform=etaConfig.uniform,
                numParticles=particleConfig.num,
                pdg=particleConfig.pdg,
                randomizeCharge=particleConfig.randomizeCharge,
                charge=particleConfig.charge,
                mass=particleConfig.mass,
                # Merging particle gun vertices does not make sense
            )
        ),
    ))

    if npileup > 0:
        import acts.examples.pythia8
        generators.append(acts.examples.EventGenerator.Generator(
            multiplicity=acts.examples.FixedMultiplicityGenerator(n=npileup),
            vertex=vtxGen,
            particles=acts.examples.pythia8.Pythia8Generator(
                level=customLogLevel(),
                **acts.examples.defaultKWArgs(
                    pdgBeam0=beam[0],
                    pdgBeam1=beam[1],
                    cmsEnergy=cmsEnergy,
                    settings=pileupProcess,
                ),
            ),
        )
        )

    evGen = EventGenerator(
        level=customLogLevel(),
        generators=generators,
        randomNumbers=rnd,
        outputEvent="particle_gun_event",
    )
    s.addReader(evGen)

    hepmc3Converter = acts.examples.hepmc3.HepMC3InputConverter(
        level=customLogLevel(),
        inputEvent=evGen.config.outputEvent,
        outputParticles="particles_generated",
        outputVertices="vertices_generated",
        mergePrimaries=False
    )
    s.addAlgorithm(hepmc3Converter)

    s.addWhiteboardAlias("particles", hepmc3Converter.config.outputParticles)
    s.addWhiteboardAlias("vertices_truth", hepmc3Converter.config.outputVertices)

    s.addWhiteboardAlias(
        "particles_generated_selected", hepmc3Converter.config.outputParticles
    )

    if printParticles:
        s.addAlgorithm(
            ParticlesPrinter(
                level=customLogLevel(),
                inputParticles=hepmc3Converter.config.outputParticles,
            )
        )

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

        s.addWriter(
            CsvParticleWriter(
                level=customLogLevel(),
                inputParticles=hepmc3Converter.config.outputParticles,
                outputDir=str(outputDirCsv),
                outputStem="particles",
            )
        )

    if outputDirRoot is not None:
        assert (
            ACTS_EXAMPLES_ROOT_AVAILABLE
        ), "ROOT output requested but ROOT is not available"
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

        s.addWriter(
            RootParticleWriter(
                level=customLogLevel(),
                inputParticles=hepmc3Converter.config.outputParticles,
                filePath=str(outputDirRoot / "particles.root"),
            )
        )

        s.addWriter(
            RootVertexWriter(
                level=customLogLevel(),
                inputVertices=hepmc3Converter.config.outputVertices,
                filePath=str(outputDirRoot / "vertices.root"),
            )
        )

    return s
