from ChainConfig import Config
cfg = Config()

#general
cfg.general.doPreVertexing = True  # turn On for Pb-Pb
cfg.general.digi_file = "geometries/geometry_2026_08_06_realistic_OT_with_eosCards_4percBarrel_4percDiscs_Gen3/digiConfigurations/digi-segm_config_pitch20um_with_barrel_TOFs-gen3_with_time.json"
cfg.general.geomGeneration = "gen3"

#pythia
cfg.pythia.usePythia = True
cfg.pythia.system = "PbPbCentral"  # "PbPbCentral"
cfg.pythia.pileup = 0

# noise
cfg.digi.injectPerPixelNoise = True
cfg.digi.noiseProbPerPixel = 1e-7

# seeding with time information
cfg.seeding.runSeedAnalysis = False
cfg.seeding.useTimeInSeeding = True

# tracking parameters (experimental)
# cfg.tracking.iterationParams[0]['pt'] = (0.2, None)
# cfg.tracking.iterationParams[1]['pt'] = (0.1, None)

# cfg.tracking.iterationParams[0]['nMeasurementsMin'] = 9
# cfg.tracking.iterationParams[1]['nMeasurementsMin'] = 7
# cfg.tracking.iterationParams[2]['nMeasurementsMin'] = 7
# cfg.tracking.iterationParams[3]['nMeasurementsMin'] = 7

# cfg.tracking.iterationParams[0]['nHoles'] = 1

# cfg.tracking.iterationParams[0]['ckfChi2Measurement'] = 25
# cfg.tracking.iterationParams[1]['ckfChi2Measurement'] = 25
# cfg.tracking.iterationParams[2]['ckfChi2Measurement'] = 25
# cfg.tracking.iterationParams[3]['ckfChi2Measurement'] = 25

# cfg.tracking.iterationParams[0]['ckfChi2Outlier'] = 100
# cfg.tracking.iterationParams[1]['ckfChi2Outlier'] = 120
# cfg.tracking.iterationParams[2]['ckfChi2Outlier'] = 120

# cfg.tracking.iterationParams[1]['nOutliers'] = 3
# cfg.tracking.iterationParams[2]['nOutliers'] = 3
# cfg.tracking.iterationParams[2]['ckfChi2Measurement'] = 25
