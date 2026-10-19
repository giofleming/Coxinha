# Very-low-pt pion gun with the GBTS graph-based seeding iteration enabled.
#
# Iterations 0-3 are the usual GridTriplet passes; iteration 4 runs
# AliceActsTrk::Alice3GbtsSeedingAlgorithm on the leftover measurements, over
# the vertex detector plus the first three inner-tracker layers - the reach of
# a track below ~150 MeV/c at B = 2 T. See docs/GBTS.md.
from ChainConfig import Config

cfg = Config()

# pythia
cfg.pythia.usePythia = False

# particle gun: push the spectrum down to where the extra iteration should pay
cfg.particleGun.gunPtRange = (0.02, 0.5)
cfg.particleGun.gunEtaRange = (-1.5, 1.5)
cfg.particleGun.gunPID = 211  # pions
cfg.particleGun.gunMult = 100

# run the GBTS iteration (iterationParams[4]) on top of the four default passes
cfg.general.iterations = 5
