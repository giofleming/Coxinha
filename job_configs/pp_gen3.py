from ChainConfig import Config

cfg = Config()

cfg.general.digi_file = "geometries/geometry_2026_08_06_realistic_OT_with_eosCards_4percBarrel_4percDiscs_Gen3/digiConfigurations/digi-segm_config_pitch20um_with_barrel_TOFs-gen3_with_time.json"
cfg.general.geomGeneration = "gen3"

cfg.general.doVertexReco = True
# cfg.general.iterations = 2

# pythia
cfg.pythia.usePythia = True
cfg.pythia.system = "pp"
cfg.pythia.pileup = 29
cfg.pythia.nBunchCrossings = 20

# noise
cfg.digi.injectPerPixelNoise = True
cfg.digi.noiseProbPerPixel = 1e-7


# seeding with time information
cfg.seeding.useTimeInSeeding = True
