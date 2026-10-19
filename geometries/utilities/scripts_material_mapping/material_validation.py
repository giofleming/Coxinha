#!/usr/bin/env python3

import argparse
import pathlib
import sys

import acts
from acts import (
    MaterialValidator,
    IntersectionMaterialAssigner,
)

from acts.examples.simulation import (
    addParticleGun,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
)

from acts.examples import (
    MaterialValidation,
    ParticleTrackParamExtractor,
)

from acts.examples.root import RootMaterialTrackWriter

u = acts.UnitConstants


class _MaterialLoader(acts.TrackingGeometryMutableVisitor):
    """Apply a material map onto an already-constructed geometry (Gen3 only).

    Gen1 injects the map through an IMaterialDecorator during construction; the
    Gen3 builder has no such hook, so the map is applied afterwards by walking
    every surface. Counts how many surfaces end up with real material so a
    map/geometry mismatch (0 matched) is caught early.
    """

    def __init__(self, decorator):
        super().__init__()
        self._deco = decorator
        self.n_seen = 0
        self.n_with_material = 0

    def visitSurface(self, surface):
        self.n_seen += 1
        self._deco.decorate(surface)
        if surface.surfaceMaterial is not None:
            self.n_with_material += 1


def runMaterialValidation(
    surfaces,
    s,
    tracksPerEvent=100,
    etaRange=(-4.2, 4.2),
    phiRange=(0.0, 360.0 * u.degree),
    materialTrackCollectionName="material_tracks",
    outputFileBase="material_validation",
    trackingGeometry=None,
):
    rnd = acts.examples.RandomNumbers(seed=42)

    addParticleGun(
        s,
        ParticleConfig(
            num=tracksPerEvent, pdg=acts.PdgParticle.eMuon, randomizeCharge=True
        ),
        EtaConfig(*etaRange),
        PhiConfig(*phiRange),
        rnd=rnd,
    )

    trkParamExtractor = ParticleTrackParamExtractor(
        level=acts.logging.INFO,
        inputParticles="particles_generated",
        outputTrackParameters="params_particles_generated",
    )
    s.addAlgorithm(trkParamExtractor)

    # Assignment setup : Intersection assigner
    materialAssingerConfig = IntersectionMaterialAssigner.Config()
    materialAssingerConfig.surfaces = surfaces
    materialAssinger = IntersectionMaterialAssigner(
        materialAssingerConfig, acts.logging.INFO
    )

    # Validator setup
    materialValidatorConfig = MaterialValidator.Config()
    materialValidatorConfig.materialAssigner = materialAssinger

    # Validation Algorithm
    materialValidationConfig = MaterialValidation.Config()
    materialValidationConfig.inputTrackParameters = "params_particles_generated"
    materialValidationConfig.materialValidator = MaterialValidator(
        materialValidatorConfig, acts.logging.INFO
    )
    materialValidationConfig.outputMaterialTracks = materialTrackCollectionName
    materialValidation = MaterialValidation(materialValidationConfig, acts.logging.INFO)
    s.addAlgorithm(materialValidation)

    # Add the material tracks writer
    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=materialValidationConfig.outputMaterialTracks,
            treeName=materialTrackCollectionName,
            filePath=str(outputFileBase) + ".root",
            storeSurface=True,
            storeVolume=True,
        )
    )

    if trackingGeometry is not None:
        nav = acts.Navigator(trackingGeometry=trackingGeometry)
        stepper = acts.StraightLineStepper()
        propagator = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

        propagationAlgorithm = acts.examples.PropagationAlgorithm(
            propagatorImpl=propagator,
            level=acts.logging.INFO,
            sterileLogger=False,
            inputTrackParameters="params_particles_generated",
            outputSummaryCollection="propagation_summary",
            outputMaterialCollection=materialTrackCollectionName + "_propagated",
        )
        s.addAlgorithm(propagationAlgorithm)

        s.addWriter(
            RootMaterialTrackWriter(
                level=acts.logging.INFO,
                inputMaterialTracks=materialTrackCollectionName + "_propagated",
                treeName=materialTrackCollectionName,
                filePath=str(outputFileBase) + "_propagated.root",
                storeSurface=True,
                storeVolume=True,
            )
        )

    return s


def main():
    p = argparse.ArgumentParser()

    p.add_argument(
        "--geo-dir",
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help="Location of the geometry",
    )
    p.add_argument(
        "-n", "--events", type=int, default=4000, help="Number of events to process"
    )
    p.add_argument(
        "-t", "--tracks", type=int, default=1000, help="Number of tracks per event"
    )
    p.add_argument("-j", "--threads", type=int, default=-1, help="Number of threads")
    p.add_argument(
        "-m",
        "--map",
        type=str,
        default="material_map.json",
        help="Input file for the material map"
    )
    p.add_argument(
        "--eta-range",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        default=(-4.2, 4.2),
        help="Eta range for generated particles",
    )
    p.add_argument(
        "--phi-range",
        nargs=2,
        type=float,
        metavar=("MIN_DEG", "MAX_DEG"),
        default=(0.0, 360.0),
        help="Phi range in degree for generated particles",
    )
    p.add_argument(
        "--material-track-collection",
        type=str,
        default="material_tracks",
        help="Output material track collection name",
    )
    p.add_argument(
        "-o",
        "--output",
        type=str,
        default="propagation-material",
        help="Output file stem (without extension)",
    )
    p.add_argument(
        "-p",
        "--propagate",
        action="store_true",
        help="Enable propagation validation",
    )
    p.add_argument(
        "--geomGeneration",
        type=str,
        default="gen1",
        choices=["gen1", "gen3"],
        help="Geometry generation to build (gen3 applies the map after "
        "construction, see below).",
    )
    p.add_argument(
        "--root-file",
        type=str,
        default="o2sim_geometry.root",
        help="gen3 only: TGeo .root file inside --geo-dir the builder reads",
    )

    args = p.parse_args()

    print("Loading geometry from ... ", str(args.geo_dir))

    if args.geomGeneration == "gen3":
        # Gen3 (Blueprint): the builder has no IMaterialDecorator hook, so the
        # map is applied *after* construction by walking every surface (Gen1
        # injects it during the build instead). Without this the surfaces keep
        # only proto material and every "mapped" plot is empty.
        from AliceActsPythonBindings import buildALICE3GeometryGen3

        rootFile = pathlib.Path(args.geo_dir) / args.root_file
        trackingGeometry = buildALICE3GeometryGen3(
            str(rootFile), acts.logging.INFO, "", True  # withMaterial=True
        )
        if args.map != "":
            decorator = acts.IMaterialDecorator.fromFile(args.map)
            loader = _MaterialLoader(decorator)
            trackingGeometry.apply(loader)
            print(
                f"  {loader.n_with_material}/{loader.n_seen} surfaces carry "
                "material after decoration"
            )
    else:
        materialDecorator = None
        if args.map != "":
            materialDecorator = acts.IMaterialDecorator.fromFile(args.map)

        sys.path.insert(0, str(pathlib.Path(args.geo_dir).resolve()))
        import buildALICE3Geometry

        detector = buildALICE3Geometry.buildALICE3Geometry(
            args.geo_dir, False, False, acts.logging.INFO, materialDecorator
        )
        trackingGeometry = detector.trackingGeometry()

    materialSurfaces = trackingGeometry.extractMaterialSurfaces()

    s = acts.examples.Sequencer(events=args.events, numThreads=args.threads)

    runMaterialValidation(
        surfaces=materialSurfaces,
        s=s,
        tracksPerEvent=args.tracks,
        etaRange=tuple(args.eta_range),
        phiRange=(args.phi_range[0] * u.degree, args.phi_range[1] * u.degree),
        materialTrackCollectionName=args.material_track_collection,
        outputFileBase=args.output,
        trackingGeometry=trackingGeometry if args.propagate else None,
    ).run()


if "__main__" == __name__:
    main()
