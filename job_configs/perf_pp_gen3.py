from ChainConfig import Config

cfg = Config()

cfg.general.digi_file = "geometries/geometry_2026_08_06_realistic_OT_with_eosCards_4percBarrel_4percDiscs_Gen3/digiConfigurations/digi-segm_config_pitch20um_with_barrel_TOFs-gen3_with_time.json"
cfg.general.geomGeneration = "gen3"

cfg.general.doVertexReco = True
# cfg.general.iterations = 2

cfg.pythia.usePythia = True
cfg.pythia.system = "pp"
cfg.pythia.pileup = 29
cfg.pythia.nBunchCrossings = 20

# noise
cfg.digi.injectPerPixelNoise = True
cfg.digi.noiseProbPerPixel = 1e-7

# seeding with time information
cfg.seeding.useTimeInSeeding = True


# Monitoring: all ROOT writers off, to profile the reconstruction
# sim / generator output (inputs of reco_only: turn back on for sim_only)
cfg.detSim.writeSimOutput = False
# digitization
cfg.digi.storeMeasurements = False
cfg.digi.storeDigiParticles = False
cfg.digi.storeNoiseMeasurements = False
# seeding, per iteration
cfg.seeding.writeSeedPerformance = False
cfg.seeding.runSeedAnalysis = False
# CKF and ambiguity resolution, per iteration
cfg.tracking.writeTrackSummary = False
cfg.tracking.writeCKFPerformance = False
cfg.tracking.writeAmbiTrackSummary = False
cfg.tracking.writeAmbiPerformance = False
# merged tracks
cfg.tracking.writeMergedTrackSummary = False
cfg.tracking.writeTrackStates = False
cfg.tracking.writePatternRecoPerformance = False
# vertexing, GSF, residuals
cfg.general.writeVertexPerformance = False
cfg.general.writeGSFTrackSummary = False
cfg.monitoring.addUnbiasedResidualWriters = False
