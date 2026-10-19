from typing import Literal, Optional, Tuple, List
from dataclasses import dataclass, field


@dataclass
class General:

    fieldMap: str = "geometries/magneticFieldMaps/solenoid_R1625_L6500_B2T_scaled.txt"
    MF: float = 2.0  # T

    # Geo 1
    # geo_dir: str = "geometries/geometry_2026_01_30_mockupTiledDisks_for_Geant_tests/geom"
    # digi_file: str = "geometries/geometry_2026_01_30_mockupTiledDisks_for_Geant_tests/digiConfigurations/digi-smearing-config_no_TOFs_iTOF_removed.json"

    # Geo 2
    # geo_dir: str = "geometries/geometry_2026_02_12_default_v3_cylindrical_barrel_tiled_disks/geom"
    # digi_file: str = "geometries/geometry_2026_02_12_default_v3_cylindrical_barrel_tiled_disks/digiConfigurations/digi-smearing-config_with_TOFs_noEndcapTOFs_noTimeInTOFs.json"

    # Marco Peacock
    # geo_dir : str = "geometries/geometry_2026_02_23_MarcoVL_peacock_iriscylinder_or_segmented/geom"
    # digi_file : str = "geometries/geometry_2026_02_23_MarcoVL_peacock_iriscylinder_or_segmented/digiConfigurations/digi-smearing-config_with_TOFs_noEndcapTOFs_noTimeInTOFs.json"

    # Segmented iTOF
    # geo_dir : str = "geometries/geometry_2026_03_13_peacock_segmented_iTof/geom"
    # digi_file : str = "geometries/geometry_2026_03_13_peacock_segmented_iTof/digiConfigurations/digi-smearing-config_with_TOFs_noEndcapTOFs_noTimeInTOFs.json"

    # Segmented June 11 with RICH
    # geo_dir : str = "geometries/geometry_2026_06_10_segmented_withRICH_noEndcapTOFs/geom"
    # digi_file : str = "geometries/geometry_2026_06_10_segmented_withRICH_noEndcapTOFs/digiConfigurations/digi-segm_config_pitch20um_with_barrel_TOFs_RICH.json"

    # Segmented June 11 without RICH
    # geo_dir : str = "geometries/geometry_2026_06_10_segmented_withoutRICH_noEndcapTOFs/geom"
    # digi_file : str = "geometries/geometry_2026_06_10_segmented_withoutRICH_noEndcapTOFs/digiConfigurations/digi-segm_config_pitch20um_with_barrel_TOFs_RICH.json"

    # Segmented June 19 without RICH, without endcap TOFs, and with FT3 stave tolerance
    # geo_dir : str = "geometries/geometry_2026_06_19_segmented_FT3stavetolerance_withoutRICH_noEndcapsTOFs/geom"
    # digi_file : str = "geometries/geometry_2026_06_19_segmented_FT3stavetolerance_withoutRICH_noEndcapsTOFs/digiConfigurations/digi-segm_config_pitch20um_with_barrel_TOFs_RICH.json"

    # Segmented August 6 - GEN3 geometry, w/o RICH, w/o endcap TOFs, end-of-stave cards for barrel and disks, OT walls
    geo_dir: str = "geometries/geometry_2026_08_06_realistic_OT_with_eosCards_4percBarrel_4percDiscs_Gen3/geom"
    digi_file: str = "geometries/geometry_2026_08_06_realistic_OT_with_eosCards_4percBarrel_4percDiscs_Gen3/digiConfigurations/digi-segm_config_pitch20um_with_barrel_TOFs-gen3.json"

    # Geometry generation:
    #   "gen1" -> python buildALICE3Geometry imported from geo_dir
    #   "gen3" -> compiled Blueprint builder buildALICE3GeometryGen3, reading
    #             geo_dir/o2sim_geometry.root; material map (if enableMaterial)
    #             from geo_dir/material_map_gen3.json
    geomGeneration: Literal["gen1", "gen3"] = "gen3"

    # Navigation policy for the Gen3 tracking geometry.
    #   "tryAll"       - ACTS blueprint default: every navigation initialisation
    #                    enumerates EVERY surface confined in the volume. With
    #                    ~19k sensors on an OT layer and ~54k on the oTOF this is
    #                    O(10^6) intersection candidates per event and dominates
    #                    CKF time (measured: 1.7M candidates in one pion-gun
    #                    event, ~0.73 us each).
    #   "surfaceArray" - binned lookup for the sensitive surfaces, plus the
    #                    portals/passives (collected once) so the material (42 of
    #                    the 45 material surfaces are passive) is still found.
    # See docs/GBTS.md section 10.
    navigationPolicy: Literal["tryAll", "surfaceArray"] = "surfaceArray"

    enableMaterial: bool = True
    iterations: int = 4  # this steers the iteration numbers. This number needs to be >1 to run extra iterations
    doVertexReco: bool = False
    # performance_vertexing.root (ACTS) and performance_vertexing_alice3.root
    writeVertexPerformance: bool = True

    doGSFtrackRefit: bool = False
    # tracksummary_gsf_refit.root
    writeGSFTrackSummary: bool = True
    # Run an independent pre-vertexing pass (seeding+CKF+ambi+vertexing) before
    # the main iterations, and feed the found vertex z-positions back to
    # constrain the seeding. Per-iteration opt-in via
    # seeding.iterationParams[i]['useVertexZ'] (+ optional 'vertexZNSigma',
    # 'vertexZMargin').
    doPreVertexing: bool = False


@dataclass
class ParticleGunConfig:
    gunMult: int = 1
    gunPtRange: Tuple[float, float] = (1.0, 1.00001)
    gunEtaRange: Tuple[float, float] = (-3, 3)
    gunPhiRange: Tuple[float, float] = (-3.1415926, 3.1415926)
    gunPID: int = 211
    gunRandCharge: bool = True
    meanVz: float = 0.
    sigmaVz: float = 50.
    # sigmaT: 0.180


@dataclass
class PythiaConfig:
    usePythia: bool = False
    system: Literal["pp", "PbPbMB", "PbPbCentral"] = "pp"
    pileup: int = 29
    meanVz: float = 0.
    sigmaVz: float = 50.

    # Discrete bunch-crossing vertex time (instead of Gaussian time).
    # When True, vertex times are drawn uniformly among nBunchCrossings slots
    # spaced by bunchSpacing (ns), starting at 0 (0, 25, ..., 475 ns). x,y,z stay Gaussian.
    useDiscreteBCtime: bool = True
    nBunchCrossings: int = 20
    bunchSpacing: float = 25.   # ns


@dataclass
class detSimConfig:
    simulation: Literal["Fatras", "Geant4"] = "Geant4"
    minFatrasPt: float = 0.001
    particleSelectionEta: float = 4.2
    storeGenAndGeantVertices = False  # set to True if we want also G4 secondary vertices, not only generator-level vertices
    # usually, no need, to save space; set to True for secondary vertexing from G4, e.g. for V0/cascade studies
    keepParticlesWithoutHits = False
    # Generator (particles.root, vertices.root) and simulation (particles_simulation.root,
    # hits.root) output. These are the inputs of a reco_only run: keep on for sim_only.
    writeSimOutput: bool = True


@dataclass
class digiConfig:
    digi_particleSelectionRho: Tuple[float, float] = (0.0, 4)  # mm keep only primaries
    digi_particleSelectionZ:   Tuple[float, float] = (0.0, 1000.)  # mm keep only primaries
    # minimum number of measurements a particle must leave to enter
    # particles_digitized_selected, i.e. the denominator of every efficiency
    particleSelectionMinHits: int = 7
    # [volume,layer], layer 7 is [(21,10)] => field(default_factory=lambda: [(21, 10)])
    particleSelectionLayers: List[Tuple[int, int]] = field(default_factory=list)
    # particleSelectionLayers: List[Tuple[int,int]] = field(default_factory=lambda: [(21, 10)])  # [volume,layer], layer 7 is [(21,10)] => field(default_factory=lambda: [(21, 10)])
    storeMeasurements = False
    storeDigiParticles = True
    # measurements_with_noise.root, written when injectPerPixelNoise is on
    storeNoiseMeasurements: bool = True

    injectPerPixelNoise: bool = False  # for studies of the effect of noise on tracking performance
    noiseProbPerPixel: float = 1e-7


@dataclass
class seedingConfig:
    # bField : float = 2.0 #T
    # seedingLayers: str = "geoSelectionForSeeding_VD_ML.json"  #helps at ~ 80 MeV
    # seedingAlgo: Literal["GridTriplet","TruthSmeared"] = "GridTriplet"
    # minSeedPt: float = 0.12
    # impParForSeeds: float = 1.0
    # sigmaScattering: float = 4
    # radLengthPerSeed: float = 0.04
    # maxSeedsPerSpM: int = 2
    # filter_maxSeedsPerSpMConf: int = 1
    # filter_maxQualitySeedsPerSpMConf: int = 1
    writeMatchingDetails: bool = False

    # Write the per-seed SeedAnalysisAlg tree (kinematics, quality, fake flags).
    runSeedAnalysis: bool = False
    # Per-iteration seeding performance (<prefix>performance_seeding.root, <prefix>estimatedparameters.root)
    writeSeedPerformance: bool = True

    # baseline seed-finder params (same for all iterations unless overridden)
    # The seed finder first forms two doublets: bottom-middle and middle-top. The deltaR* cuts limit the allowed radial gap between the two points of each doublet.
    # deltaR — the *general* min/max radial distance between any two doublet space points.
    # deltaRBottomSP: (3., 50.) — the min/max radial distance for the *inner-middle* doublet. If set (not NaN), it overrides the general deltaR for that doublet.
    deltaR: Tuple[float, float] = (1., 200.)          # mm
    # deltaRBottomSP: Tuple[float, float] = (3., 50.)   # mm
    collisionRegion: float = 250.  # mm
    cotThetaMax: float = 13.538    # 10.018 ==> eta=3.0 #13.538 ==> eta=3.3      #27.2899 ==> eta=4
    seedConfirmation: bool = True
    centralSeedConfirmation: dict = field(default_factory=lambda: {
        'zMinSeedConf': -620.,
        'zMaxSeedConf': 620.,
        'rMaxSeedConf': 4.9,
        'nTopForLargeR': 1,
        'nTopForSmallR': 2,
    })
    forwardSeedConfirmation: dict = field(default_factory=lambda: {
        'zMinSeedConf': -1220.,
        'zMaxSeedConf': 1220.,
        'rMaxSeedConf': 26.,
        'nTopForLargeR': 1,
        'nTopForSmallR': 2,
    })
    useVariableMiddleSPRange: bool = False
    deltaRMiddleSPRange: Tuple[float, float] = (0.2, 1.0)  # mm

    # use time information in seeding
    useTimeInSeeding: bool = True

    # Sigmas (loc0, loc1, phi, theta, q/p, t) of the covariance the CKF starts from, in
    # ACTS units (use acts.UnitConstants: time is in mm/c, 1 ns = 299.79). None = the default
    # in writers.addIterativeTracking: 1 mm, 1 mm, 1 deg, 1 deg, 0.1 e/GeV, 100 ns.
    initialSigmas: Optional[List[float]] = None

    # Dedicated, lighter seeding config for the independent pre-vertexing pass
    # (general.doPreVertexing). High minPt + few seeds -> fast and high-purity,
    # which is all that is needed to localize the primary vertices in z. Must
    # carry the same keys as an iterationParams entry.
    preVertexingParams: dict = field(default_factory=lambda: {
        'minPt': 0.7, 'impactMax': 1.5, 'z': (-1400, 1400), 'r': (None, 350),
        'maxSeedsPerSpM': 1, 'sigmaScattering': 5, 'radLengthPerSeed': 0.05,
        'deltaRTopSP': (20., 150.),
        'deltaRBottomSP': (20., 50.),
        'seedingConfigFile': "geoSelectionForSeeding_3firstML.json",
        'rRangeMiddleSP': [[0.0, 0.0],
                           [200., 350.],
                           [80., 350.],
                           [80., 350.],
                           [80., 350.],
                           [80., 350.],
                           [80., 350.],
                           [200., 350.],
                           [0.0, 0.0],
                           ]
    })

    iterationParams: dict = field(default_factory=lambda: {
        # iter 0: find high-pt tracks using ML seeds
        0: {
            'minPt': 0.4, 'impactMax': 1.5, 'z': (-1400, 1400), 'r': (None, 350),
            'maxSeedsPerSpM': 3, 'sigmaScattering': 6, 'radLengthPerSeed': 0.06,
            'deltaRTopSP': (10., 200.),
            'deltaRBottomSP': (10., 50.),
            'seedingConfigFile': "geoSelectionForSeeding_3firstML.json",
            'rRangeMiddleSP': [[0.0, 0.0],   # -3000,-1400 (unused)
                               [200., 350.],  # -1400.,-900
                               [80., 350.],  # -900., -700
                               [80., 350.],  # -700., -250
                               [80.,  350.],  # -250.,  250.
                               [80., 350.],  # 250.,  700
                               [80., 350.],  # 700.,  900
                               [200., 350.],  # 900., 1400
                               [0.0, 0.0],   # 1400, 3000 (unused)
                               ],
            # constrain this iteration's seeding (ignored if pre-vertexing iteration is Off)
            'useVertexZ': True,
            'vertexZNSigma': 3.0,        # half-window = N*sigma_z ...
            'vertexZMargin': 2.5,        # additional +/- margin (mm)
        },
        # iter 1: find lower-pt tracks using VD seeds, with strict impact parameter constraint
        1: {
            'minPt': 0.15, 'impactMax': 0.1, 'z': (-1400, 1400), 'r': (None, 27),
            'maxSeedsPerSpM': 4, 'sigmaScattering': 2, 'radLengthPerSeed': 0.02,
            'deltaRTopSP': (3., 200.),
            'deltaRBottomSP': (3., 50.),
            'seedingConfigFile': "geoSelectionForSeeding_VD.json",
            'rRangeMiddleSP': [[0.0, 0.0],   # -3000,-1400 (unused)
                               [200., 350.],  # -1400.,-900
                               [80., 350.],  # -900., -700
                               [80., 250.],  # -700., -250
                               [6.,  25.],   # -250.,  250.
                               [80., 250.],  # 250.,  700
                               [80., 350.],  # 700.,  900
                               [200., 350.],  # 900., 1400
                               [0.0, 0.0],   # 1400, 3000 (unused)
                               ],
            # constrain this iteration's seeding (ignored if pre-vertexing iteration is Off)
            'useVertexZ': True,
            'vertexZNSigma': 3.0,        # half-window = N*sigma_z ...
            'vertexZMargin': 0.5,        # additional +/- margin (mm)
        },
        # iter 2: find even lower-pt tracks using VD seeds, with relaxed impact parameter constraint
        2: {
            'minPt': 0.04, 'impactMax': 0.25, 'z': (-1400, 1400), 'r': (None, 27),
            'maxSeedsPerSpM': 4, 'sigmaScattering': 1.5, 'radLengthPerSeed': 0.01,
            'deltaRTopSP': (3., 200.),
            'deltaRBottomSP': (3., 50.),
            'seedingConfigFile': "geoSelectionForSeeding_VD.json",
            'rRangeMiddleSP': [[0.0, 0.0],   # -3000,-1400 (unused)
                               [200., 350.],  # -1400.,-900
                               [80., 350.],  # -900., -700
                               [80., 250.],  # -700., -250
                               [6.,  25.],   # -250.,  250.
                               [80., 250.],  # 250.,  700
                               [80., 350.],  # 700.,  900
                               [200., 350.],  # 900., 1400
                               [0.0, 0.0],   # 1400, 3000 (unused)
                               ],
            # constrain this iteration's seeding (ignored if pre-vertexing iteration is Off)
            'useVertexZ': True,
            'vertexZNSigma': 3.0,        # half-window = N*sigma_z ...
            'vertexZMargin': 1.0,        # additional +/- margin (mm)
        },
        # iter 3: recover remaining lower-pt tracks using ML seeds, with loose impact parameter constraint
        3: {
            'minPt': 0.1, 'impactMax': 4.5, 'z': (-1400, 1400), 'r': (None, 210),
            'maxSeedsPerSpM': 3, 'sigmaScattering': 5, 'radLengthPerSeed': 0.05,
            'deltaRTopSP': (10., 200.),
            'deltaRBottomSP': (10., 50.),
            'seedingConfigFile': "geoSelectionForSeeding_3firstML.json",
            'rRangeMiddleSP': [[0.0, 0.0],   # -3000,-1400 (unused)
                               [200., 350.],  # -1400.,-900
                               [80., 350.],  # -900., -700
                               [80., 350.],  # -700., -250
                               [80.,  350.],  # -250.,  250.
                               [80., 350.],  # 250.,  700
                               [80., 350.],  # 700.,  900
                               [200., 350.],  # 900., 1400
                               [0.0, 0.0],   # 1400, 3000 (unused)
                               ],
            # constrain this iteration's seeding (ignored if pre-vertexing iteration is Off)
            'useVertexZ': True,
            'vertexZNSigma': 3.0,        # half-window = N*sigma_z ...
            'vertexZMargin': 5.0,        # additional +/- margin (mm)
        },
        # iter 4: very-low-pt recovery with GBTS on the vertex detector + the
        # first three inner-tracker layers (everything a track below ~150 MeV/c
        # can reach at 2 T). Only runs when general.iterations >= 5.
        # The GridTriplet keys below are unused by GBTS but are still read when
        # the SeedFinderConfigArg is assembled; 'seedingConfigFile' is what
        # selects the space points, and must match cfg.gbts layer set.
        4: {
            'seedingAlgorithm': 'Gbts',
            'gbts': {
                'layerSet': 'layersVDIT3',
                'minPt': 0.05,
                'maxOuterRadius': 130.0,
            },
            'minPt': 0.05, 'impactMax': 1.0, 'z': (-1400, 1400), 'r': (None, 130),
            'maxSeedsPerSpM': 4, 'sigmaScattering': 2, 'radLengthPerSeed': 0.02,
            'deltaRTopSP': (3., 60.),
            'deltaRBottomSP': (3., 30.),
            'seedingConfigFile': "geoSelectionForSeeding_VD_IT3.json",
            'rRangeMiddleSP': [[0.0, 0.0],   # -3000,-1400 (unused)
                               [0.0, 0.0],   # -1400.,-900
                               [6., 122.],   # -900., -700
                               [6., 122.],   # -700., -250
                               [6., 122.],   # -250.,  250.
                               [6., 122.],   # 250.,  700
                               [6., 122.],   # 700.,  900
                               [0.0, 0.0],   # 900., 1400
                               [0.0, 0.0],   # 1400, 3000 (unused)
                               ],
            'useVertexZ': False,
        },
    })


@dataclass
class gbtsConfig:
    """GBTS (graph-based track seeding) settings; see docs/GBTS.md.

    Used only by iterations whose seeding.iterationParams entry carries
    'seedingAlgorithm': 'Gbts'. Per-iteration overrides go in that entry's
    'gbts' sub-dict and are layered on top of the defaults here.

    The (volume, layer) pairs are ACTS GeometryIdentifier components of the
    *Gen3* geometry; regenerate them with scripts/dump_gbts_layers.py if the
    geometry changes.
    """

    # ---- layer sets -----------------------------------------------------
    # Vertex detector only (r = 5, 12, 25 mm). Three layers give chains of at
    # most 3 space points, i.e. level 2, so this set needs minSeedLevel = 2.
    layersVD: List[Tuple[int, int]] = field(default_factory=lambda: [
        (40, 1), (40, 3), (40, 5),
    ])
    # Vertex detector + the first three inner-tracker layers (r <= 121 mm).
    # The reach of a ~150 MeV/c track at B = 2 T, and the default low-pt set.
    layersVDIT3: List[Tuple[int, int]] = field(default_factory=lambda: [
        (40, 1), (40, 3), (40, 5),
        (40, 9), (40, 11), (40, 13),
    ])
    # Vertex detector + the first three inner-tracker layers + the three inner
    # FT3 discs per side. The forward-capable high-pt set: it reaches eta ~3.0
    # seedable (>=4 layers), which is all the forward acceptance the full
    # 23-layer set provides (3.05), at 12 layers instead of 23. Above eta ~3 the
    # limiting layers are the vertex detector, not the outer discs or the OT, so
    # those add combinatorics and CKF load without extending the reach.
    # Pairs with geoSelectionForSeeding_VD_3ML_innerDisks.json.
    layersVD3MLDiscs: List[Tuple[int, int]] = field(default_factory=lambda: [
        (40, 1), (40, 3), (40, 5),
        (40, 9), (40, 11), (40, 13),
        (30, 11), (30, 7), (30, 3),
        (60, 3), (60, 7), (60, 11),
    ])

    # Whole barrel + forward tracker, for the "replace the seeding" mode.
    # The barrel TOF (40, 17) and the outer TOF (90, 3) are deliberately left
    # out: they are timing layers, not tracking ones.
    layersFull: List[Tuple[int, int]] = field(default_factory=lambda: [
        (40, 1), (40, 3), (40, 5),
        (40, 9), (40, 11), (40, 13), (40, 15), (40, 19),
        (50, 3), (50, 5), (50, 7),
        (30, 11), (30, 7), (30, 3),
        (20, 11), (20, 9), (20, 7),
        (60, 3), (60, 7), (60, 11),
        (70, 4), (70, 6), (70, 8),
    ])

    # ---- geometry / graph schedule --------------------------------------
    # Which of the sets above an iteration uses, unless it overrides 'layers'.
    layerSet: str = "layersVDIT3"
    # How many layers a generated connection may skip. 1 = consecutive only,
    # 2 = tolerate one missing layer.
    maxLayerSkip: int = 2
    connectionScanZ0Steps: int = 21
    connectionScanEtaSteps: int = 2001
    # Width of the eta bins each GBTS layer is split into. Coarser bins suit
    # the short lever arm of a low-pt pass.
    etaBinWidth: float = 0.3
    # Luminous region the eta-bin table is built against (mm).
    #
    # THE most important GBTS setting for ALICE 3. It gates twice: per doublet
    # (z0 = z1 - r1*tau must lie in the window) and, through
    # GbtsLayer::checkCompatibility, on which eta-bin pairs may be connected at
    # all. Set too narrow, the graph does not just lose tracks, it falls apart:
    # with sigmaVz = 50 mm, a +-60 mm window gave 986 edges / 5 links / 0 seeds,
    # while +-150 mm gave 5506 edges / 4784 links / 99 seeds and ~100% seeding
    # efficiency on the same event. Keep it at ~3 sigma of the vertex z spread
    # (cfg.particleGun.sigmaVz or cfg.pythia.sigmaVz, both 50 mm by default).
    binTableZ0: Tuple[float, float] = (-150.0, 150.0)

    # ---- region of interest ---------------------------------------------
    roiEta: Tuple[float, float] = (-4.0, 4.0)
    # mm; also the seeder's own z0 cut. Same reasoning as binTableZ0 above.
    roiZ: Tuple[float, float] = (-150.0, 150.0)

    # ---- seed finder ----------------------------------------------------
    # The default cut coefficients were tuned at tuningPt = 0.9 GeV and scale
    # by tuningPt / minPt, so lowering minPt widens the phi and curvature
    # windows on its own.
    minPt: float = 0.05            # GeV
    # Chain of minSeedLevel + 1 space points. The ACTS default of 3 asks for
    # four, which a track that only reaches the vertex detector can never give;
    # 2 (a triplet) is what the low-pt layer set can actually deliver.
    minSeedLevel: int = 2
    minDeltaRadius: float = 2.0    # mm; VD layer spacing is 7 mm
    maxAbsTau: float = 27.3        # |eta| < 4
    maxOuterRadius: float = 130.0  # mm; outer edge of the layer set
    nMaxPhiSlice: float = 53.0
    nMaxEdges: int = 2000000
    hitShareThreshold: float = 0.49
    validateTriplets: bool = True
    useAdaptiveCuts: bool = True
    # addTriplets lowers the required chain length by one (i.e. accepts a bare
    # triplet) for |eta| < maxAbsEtaAddTriplets. The ACTS default of 1.5 keeps
    # it central-only, which is useless forward: it has to be opened past
    # edgeMaskMinEta to do anything, and the forward region is precisely where
    # there are too few layers to reach minSeedLevel. See the acceptance table
    # in docs/GBTS.md section 9.
    addTriplets: bool = False
    maxAbsEtaAddTriplets: float = 1.5
    # |eta| above which the edges of an extracted seed are marked consumed, so
    # another candidate cannot reuse them. Controls forward duplicates.
    edgeMaskMinEta: float = 1.5
    d0Max: float = 3.0             # mm
    # The z0-histogram cut on the innermost barrel layer is tuned for ITk
    # occupancy; -1 disables it.
    z0HistogramMaxBarrelOrder: int = -1

    # ---- edge-to-edge linking tolerances --------------------------------
    # ACTS scales its curvature and phi *windows* by tuningPt / minPt, but not
    # these three, which decide whether two edges sharing a node are linked
    # into a triplet. Left absolute they are the binding cut at low pt: at
    # 50 MeV a real track's curvature proxy is ~6e-3/mm, six times cutDCurvMax.
    # With scaleLinkingCutsWithPt the same tuningPt / minPt factor is applied
    # to them here. The values below are the ACTS defaults, before scaling.
    # Demand that the inner node already carries an incoming edge whose tau
    # agrees within tauRatioPrecut BEFORE a new edge is created, on the
    # innermost pixel-barrel layers (up to matchBeforeCreateMaxBarrelOrder).
    # Designed to suppress combinatorial edges where occupancy is worst.
    matchBeforeCreate: bool = False
    tauRatioPrecut: float = 0.009

    scaleLinkingCutsWithPt: bool = True
    tuningPt: float = 0.9          # GeV; what the defaults were tuned at
    tauRatioCut: float = 0.007
    tauRatioCorr: float = 0.006
    cutDPhiMax: float = 0.012
    cutDCurvMax: float = 0.001     # 1/mm
    # Track-level bound the chain-following filter applies after each update.
    # None derives it from minPt: |curvature| ~ 1 / (2 R) with a safety factor.
    filterMaxCurvature: Optional[float] = None   # 1/mm
    filterMaxCurvatureSafety: float = 2.0

    # ---- chain-following filter -----------------------------------------
    # Defaults are ITk pixel resolutions; ALICE 3 has a 20 um pitch, so these
    # can be tightened once the chain is understood.
    filterSigmaX: float = 0.08     # mm
    filterSigmaY: float = 0.25     # mm
    filterMaxDChi2X: float = 5.0
    filterMaxDChi2Y: float = 6.0
    filterAddHit: float = 14.0
    filterRadLen: float = 0.025    # radiation lengths per layer

    # The filter inflates for multiple scattering as
    #     corrMS = sigmaMS * |curvature| / maxCurvature
    # so sigmaMS and maxCurvature are not independent: raising maxCurvature to
    # admit a low-pt track *shrinks* the scattering it is allowed, which is
    # backwards. What is physical is the ratio, so it is what is configured
    # here and sigmaMS is derived from it. The default is the ACTS one,
    # 0.016 / (1e-3 /mm).
    filterMsPerCurvature: float = 16.0   # mm


@dataclass
class trackingConfig:
    minPt: float = 0.05  # GeV
    # used in final ambiguity resolution for merging tracks from different iterations, keep it low to cover also low-pt tracking settings
    nMeasurementsMinForFinalAmbiRes: int = 5
    maxSharedHits: int = 2
    ckfMeasPerSurf: int = 1
    seedDeduplication: bool = True
    stayOnSeed: bool = False
    twoWayCKF: bool = True
    truthMatchingRatio: float = 0.5
    # ROOT output. Per iteration (iter_<n>/ when iterations > 1):
    writeTrackSummary: bool = True           # tracksummary_<ckf tracks>.root
    writeCKFPerformance: bool = True         # performance_finding_<ckf tracks>.root
    writeAmbiTrackSummary: bool = True       # tracksummary_<ambi tracks>_writer.root
    writeAmbiPerformance: bool = True        # performance_finding_<ambi tracks>_writer.root
    writeMatchingDetailsCKF: bool = False
    writeMatchingDetailsAmbi: bool = False
    # After merging the iterations (only when iterations > 1):
    writeMergedTrackSummary: bool = True     # tracksummary_ambi-tracks-merged.root
    writeTrackStates: bool = True            # trackstates_ambi-tracks-merged.root
    writePatternRecoPerformance: bool = True  # performance_merged_seed.root, performance_merged_ambi_tracks.root

    # CKF track state creator (the CKF step that calibrates and chi2-tests the
    # measurements on a surface):
    #   "default" - ACTS: every measurement on the surface is tested. On the VD
    #               (3 full cylinders, ~2k measurements each in pp + pileup)
    #               this dominates the CKF time.
    #   "gridded" - surfaces with >= trackStateCreatorGridMinMeasurements
    #               measurements are binned in (loc0, loc1); only the ones in a
    #               window around the prediction are tested. The window is exact
    #               (see GriddedTrackStateCreator.hpp): same tracks, less time.
    #               Not exact with stayOnSeed = True, which is then refused.
    #               Also implements the per-iteration 'forceSeedMeasurements'.
    trackStateCreator: Literal["default", "gridded"] = "gridded"
    trackStateCreatorGridMinMeasurements: int = 64
    trackStateCreatorGridMeasurementsPerCell: float = 1.0

    # The CKF runs in AliceTrackFindingAlgorithm, actso2's copy of the ACTS
    # examples TrackFindingAlgorithm (same tracks with the same options).
    # The second, inward, pass of the two-way CKF inherits the
    # hole/outlier counts of the first pass, including the trailing holes the
    # final track will not have. A first pass stopped by holes therefore stops
    # the second pass at its first surface. True starts the second pass from
    # the trimmed counts instead. Changes the tracks.
    secondPassTrimmedCounts: bool = False
    # e.g. to match displaced reco tracks to sim-level weak decay daughters
    doTruthMatchingAlsoWithParticlesFromGeant: bool = True

    # 'forceSeedMeasurements' (needs trackStateCreator "gridded"; never applied
    # in the pre-vertexing pass):
    # as ATLAS forceTrackOnSeed, on a surface where the seed has a measurement
    # that measurement is the only candidate and is accepted without chi2 cut.
    # On for the VD-seeded iterations, where the seed hits always pass the
    # chi2 anyway and it saves the wide start-covariance search (PbPb: CKF of
    # iteration 1 -65 %, efficiency up, fakes down). Off for the ML-seeded
    # iterations: there the chi2 on the seed's own hits rejects inconsistent
    # (fake) seeds, and forcing doubles the iteration-0 fake ratio in PbPb.
    # See perf_study/README.md section 9.
    iterationParams: dict = field(default_factory=lambda: {
        0: {'pt': (0.05, None), 'nMeasurementsMin': 7, 'nHoles': 2, 'nOutliers': 2, 'ckfChi2Measurement': 45, 'ckfChi2Outlier': 100,
            'forceSeedMeasurements': False},
        1: {'pt': (0.05, None), 'nMeasurementsMin': 7, 'nHoles': 2, 'nOutliers': 2, 'ckfChi2Measurement': 45, 'ckfChi2Outlier': 100,
            'forceSeedMeasurements': True},
        2: {'pt': (0.05, None), 'nMeasurementsMin': 7, 'nHoles': 2, 'nOutliers': 2, 'ckfChi2Measurement': 45, 'ckfChi2Outlier': 100,
            'forceSeedMeasurements': True},
        3: {'pt': (0.05, None), 'nMeasurementsMin': 7, 'nHoles': 2, 'nOutliers': 2, 'ckfChi2Measurement': 45, 'ckfChi2Outlier': 100,
            'forceSeedMeasurements': False},
        # Iteration 4 is the GBTS very-low-pt pass (only read when
        # general.iterations >= 5). nMeasurementsMin is 5, not 7: a 50-100 MeV
        # pion at 2 T crosses at most ~6 of the six layers in the VD+IT3 set
        # before curling, so a 7-measurement track is not something it can
        # produce, and the 7 cut discards ~half of correctly seeded tracks.
        # NOTE the same relaxation applied to the GridTriplet iterations 2 and 3
        # recovers the same tracks without GBTS at all - see docs/GBTS.md
        # section 7. Re-check the fake rate on pp/PbPb before relying on it.
        4: {'pt': (0.05, None), 'nMeasurementsMin': 5, 'nHoles': 2, 'nOutliers': 2, 'ckfChi2Measurement': 45, 'ckfChi2Outlier': 100,
            'forceSeedMeasurements': False},
    })


@dataclass
class alignmentConfig:
    # Misalign the reconstruction geometry (sim + digi stay nominal), see
    # alice3/performance/alignment.py::addGlobalShifts. Works with full_chain and
    # reco_only (reuse a nominal sim sample); no effect in sim_only.
    testMisalignment: bool = False
    # Global translations in mm, one entry per structure (gen3 ids: see the
    # builder's layer table). layer/sensitive omitted or 0 = whole volume/layer;
    # entries hitting the same surface add up. Example, OT barrel L2 by 100 um in x:
    #   [{'volume': 50, 'layer': 5, 'shift': (0.1, 0., 0.)}]
    globalShifts: List[dict] = field(default_factory=list)


@dataclass
class monitoring:
    addUnbiasedResidualWriters: bool = False


class Config:
    general = General()
    pythia = PythiaConfig()
    detSim = detSimConfig()
    digi = digiConfig()
    particleGun = ParticleGunConfig()
    seeding = seedingConfig()
    gbts = gbtsConfig()
    tracking = trackingConfig()
    alignment = alignmentConfig()
    monitoring = monitoring()
