import acts
import acts.examples
import acts.examples.root

# from AliceActsPythonBindings import ResPlotToolConfig
from AliceActsPythonBindings import UnbiasedResidualPlotToolConfig

# HitOnTrack Residuals
unbiasedResidualPlotToolConfig = UnbiasedResidualPlotToolConfig()
binning = unbiasedResidualPlotToolConfig.varBinning
binning["ResidualLoc0"] = acts.Axis.regular(200, -2000, 2000, "res loc0 [um]")
binning["ResidualLoc1"] = acts.Axis.regular(200, -2000, 2000, "res loc1 [um]")
unbiasedResidualPlotToolConfig.varBinning = binning

# Residual Plots
resPlotToolConfig = acts.examples.ResPlotToolConfig()
# Update the defaults in place: ResPlotTool also needs binnings for the q/pT
# parameters ("Residual_qopt", "Residual_qopt_rel") and throws if any key is
# missing, so replacing the whole map wholesale is not safe.
binning = resPlotToolConfig.varBinning
binning["Eta"] = acts.Axis.regular(80,   -4,    4,     "#eta")
binning["Phi"] = acts.Axis.regular(70,   -3.5,  3.5,   "#phi")
binning["Pt"] = acts.Axis.regular(1000,  0,    50,    "pT [GeV/c]")
binning["Residual_phi"] = acts.Axis.regular(200,  -0.02, 0.02,  "r_{#phi} [rad]")
binning["Residual_theta"] = acts.Axis.regular(200,  -0.02, 0.02,  "r_{#theta} [rad]")
binning["Residual_qop"] = acts.Axis.regular(200,  -0.2,  0.2,   "r_{q/p} [c/GeV]")
binning["Residual_qopt"] = acts.Axis.regular(100,  -0.1,  0.1,   "r_{q/pT} [c/GeV]")
binning["Residual_qopt_rel"] = acts.Axis.regular(100,  -0.1,  0.1,   "r_{rel q/pT} [%]")
binning["Pull"] = acts.Axis.regular(100,  -5,    5,     "pull")
binning["Residual_d0"] = acts.Axis.regular(100,  -0.5,  0.5,   "r_{d0} [mm]")
binning["Residual_z0"] = acts.Axis.regular(100,  -0.5,  0.5,   "r_{z0} [mm]")
binning["Residual_t"] = acts.Axis.regular(100,  -1000, 1000,  "r_{t} [mm/c]")
resPlotToolConfig.varBinning = binning


# Duplication Plots
duplicationPlotToolConfig = acts.examples.DuplicationPlotToolConfig()
binning = duplicationPlotToolConfig.varBinning
binning["Eta"] = acts.Axis.regular(80, -4,    4,     "#eta")
binning["Phi"] = acts.Axis.regular(100, -3.15, 3.15, "#phi")
binning["Pt"] = acts.Axis.regular(1000,  0,    50,   "pT [GeV/c]")
binning["Num"] = acts.Axis.regular(30, -0.5,  29.5,  "N")
duplicationPlotToolConfig.varBinning = binning

# Efficiency Plots
effPlotToolConfig = acts.examples.EffPlotToolConfig()
binning = effPlotToolConfig.varBinning
binning["Eta"] = acts.Axis.regular(80, -4,    4,     "#eta")
binning["Pt"] = acts.Axis.regular(1000,  0,    10,   "pT [GeV/c]")
binning["D0"] = acts.Axis.regular(500, -0.0005, 0.0005, "d_0 [mm]")
effPlotToolConfig.varBinning = binning
effPlotToolConfig.minTruthPt = 0.02

# Fake rate Plots
fakePlotToolConfig = acts.examples.FakePlotToolConfig()
binning = fakePlotToolConfig.varBinning
binning["Eta"] = acts.Axis.regular(80, -4,    4,     "#eta")
binning["Pt"] = acts.Axis.regular(1000,  0,    50,   "pT [GeV/c]")
fakePlotToolConfig.varBinning = binning

# TrackQuality Plots
trackQualityPlotToolConfig = acts.examples.TrackQualityPlotToolConfig()
binning = trackQualityPlotToolConfig.varBinning
binning["Eta"] = acts.Axis.regular(80, -4,    4,     "#eta")
binning["Pt"] = acts.Axis.regular(1000,  0,    50,   "pT [GeV/c]")
trackQualityPlotToolConfig.varBinning = binning


# Track Summary Plots
trackSummaryPlotToolConfig = acts.examples.TrackSummaryPlotToolConfig()
binning = trackSummaryPlotToolConfig.varBinning
binning["Eta"] = acts.Axis.regular(80, -4,    4,     "#eta")
binning["Pt"] = acts.Axis.regular(1000,  0,    50,   "pT [GeV/c]")
trackSummaryPlotToolConfig.varBinning = binning
