from ChainConfig import Config
cfg = Config()

# general
cfg.general.iterations = 3
# cfg.general.fieldMap = "" # if no map is provided, a constant MF along z is used
# cfg.general.MF = 0.5 # T  # constant MF along z

# pythia
cfg.pythia.usePythia = False

# particle gun
cfg.particleGun.gunPtRange = (0.01, 2.0)
cfg.particleGun.gunEtaRange = (-0.8, 0.8)
cfg.particleGun.gunPID = 11  # electrons

# example how to change one parameter in an existing iteration:
# cfg.seeding.iterationParams[2]['impactMax'] = 24.0

# last iteration: relaxed parameters to recover more electrons
cfg.seeding.iterationParams[2] = {
    'minPt': 0.01,
    'impactMax': 0.8,
    'z': (-300, 300),
    'r': (None, 27),
    'maxSeedsPerSpM': 2,
    'sigmaScattering': 3,
    'radLengthPerSeed': 0.02,
    'seedingConfigFile': "geoSelectionForSeeding_VD.json",
}
