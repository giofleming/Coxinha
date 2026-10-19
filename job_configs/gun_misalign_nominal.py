# Misalignment study, nominal reference: barrel muons for the momentum resolution.
# Compare with the misaligned configs (e.g. gun_misalign_OTBarrelL2_Tx100um.py),
# which import this one and only add cfg.alignment.
from ChainConfig import Config
cfg = Config()

# particle gun
cfg.pythia.usePythia = False
cfg.particleGun.gunPID = 13  # mu-
cfg.particleGun.gunRandCharge = True
cfg.particleGun.gunMult = 10
cfg.particleGun.gunPtRange = (0.5, 10.0)  # GeV
cfg.particleGun.gunEtaRange = (-1.0, 1.0)  # barrel: crosses all OT barrel layers

cfg.detSim.simulation = "Geant4"
