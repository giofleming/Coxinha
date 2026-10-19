from ChainConfig import Config
cfg = Config()

# general
cfg.general.iterations = 4  # additional iteration for V0, with loose cuts

# pythia
cfg.pythia.usePythia = False

# particle gun
cfg.particleGun.gunPtRange = (1.0, 1.0)
cfg.particleGun.gunEtaRange = (0.4, 0.4)
cfg.particleGun.gunRandCharge = False

# cfg.particleGun.gunPID = 3122 # Lambda0
cfg.particleGun.gunPID = 310  # K0S

# detector response
cfg.detSim.storeGenAndGeantVertices = True
cfg.detSim.keepParticlesWithoutHits = True

# seeding
cfg.seeding.deltaR = (1, 250)  # mm
cfg.seeding.deltaRBottomSP = (3, 310)  # mm

# example how to change one parameter in an existing iteration:
# cfg.seeding.iterationParams[2]['impactMax'] = 24.0

# add special iteration 3 - relaxed parameters to find seeds for V0 daughters (!)
cfg.seeding.iterationParams[3] = {
    'minPt': 0.12,
    'impactMax': 1500.,
    'z': (-2400, 2400),      # wider z to cover seeds in disks
    'r': (None, 350.),
    'maxSeedsPerSpM': 4,
    'sigmaScattering': 8.0,
    'radLengthPerSeed': 0.05,
    'seedingConfigFile': "geoSelectionForSeeding_5ML.json",
}

# tracking
# cfg.tracking.doTruthMatchingAlsoWithParticlesFromGeant = True # this is already the default in ChainConfig.py
