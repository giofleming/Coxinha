from ChainConfig import Config
cfg = Config()

# pythia
cfg.pythia.usePythia = False

# particle gun
cfg.particleGun.gunPtRange = (0.01, 10.0)
cfg.particleGun.gunEtaRange = (-0.8, 0.8)
cfg.particleGun.gunPID = 211  # pions
cfg.particleGun.gunMult = 100
