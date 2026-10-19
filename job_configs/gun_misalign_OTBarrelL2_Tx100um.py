# Misalignment study: OuterTrackerBarrel_L2 (gen3 vol 50, layer 5) shifted by
# +100 um in global x in the reconstruction geometry. Sim/digi stay nominal, so
# the real detector is effectively displaced by -100 um (see addGlobalShifts).
from gun_misalign_nominal import cfg

cfg.alignment.testMisalignment = True
cfg.alignment.globalShifts = [
    {'volume': 50, 'layer': 5, 'shift': (0.1, 0., 0.)},  # mm
]
