#!/usr/bin/env python3

import argparse
from pathlib import Path
import pathlib
import sys
import acts

from acts import (
    Surface,
    MaterialMapper,
    IntersectionMaterialAssigner,
    BinnedSurfaceMaterialAccumulator,
    logging,
    GeometryContext,
)

from acts.json import MaterialMapJsonConverter

from acts.examples import (
    Sequencer,
    WhiteBoard,
    MaterialMapping,
)

from acts.examples.root import (
    RootMaterialTrackReader,
    RootMaterialTrackWriter,
    RootMaterialWriter,
)

from acts.examples.json import (
    JsonMaterialWriter,
    JsonFormat,
)


def runMaterialMapping(
    surfaces: list[Surface],
    inputFile: Path,
    nEvents: int,
    outputFileBase: str,
    outputMapFormats: list[str] = ["json", "root"],
    loglevel: acts.logging.Level = acts.logging.INFO,
    outputMaterialTracks: str = "material_tracks",
    treeName: str = "material_tracks",
):
    # Create a sequencer
    print("Creating the sequencer with 1 thread (inter event information needed)")

    s = Sequencer(numThreads=1, events=nEvents)

    # IO for material tracks reading
    wb = WhiteBoard(acts.logging.INFO)

    # Read material step information from a ROOT TTRee
    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            outputMaterialTracks=outputMaterialTracks,
            treeName=treeName,
            fileList=[str(inputFile)],
            readCachedSurfaceInformation=False,
        )
    )

    # Assignment setup : Intersection assigner
    materialAssingerConfig = IntersectionMaterialAssigner.Config()
    materialAssingerConfig.surfaces = surfaces
    materialAssinger = IntersectionMaterialAssigner(materialAssingerConfig, loglevel)

    # Accumulation setup : Binned surface material accumulator
    materialAccumulatorConfig = BinnedSurfaceMaterialAccumulator.Config()
    materialAccumulatorConfig.materialSurfaces = surfaces
    materialAccumulator = BinnedSurfaceMaterialAccumulator(
        materialAccumulatorConfig, loglevel
    )

    # Mapper setup
    materialMapperConfig = MaterialMapper.Config()
    materialMapperConfig.assignmentFinder = materialAssinger
    materialMapperConfig.surfaceMaterialAccumulator = materialAccumulator
    materialMapper = MaterialMapper(materialMapperConfig, loglevel)

    # Add the map writer(s)
    materialMapWriters = []
    # json map writer
    if "json" in outputMapFormats:
        jmConverterCfg = MaterialMapJsonConverter.Config(
            processSensitives=True,
            processApproaches=True,
            processRepresenting=True,
            processBoundaries=True,
            processVolumes=False,
        )
        # Suffix for the map file is added in the writer depending on the format
        materialMapWriters.append(
            JsonMaterialWriter(
                level=loglevel,
                converterCfg=jmConverterCfg,
                fileName=outputFileBase + "_map",
                writeFormat=JsonFormat.Json,
            )
        )
    if "root" in outputMapFormats:
        materialMapWriters.append(
            RootMaterialWriter(
                level=loglevel,
                filePath=outputFileBase + "_map.root",
            )
        )

    # Mapping Algorithm
    materialMappingConfig = MaterialMapping.Config()
    materialMappingConfig.materialMapper = materialMapper
    materialMappingConfig.inputMaterialTracks = outputMaterialTracks
    materialMappingConfig.mappedMaterialTracks = outputMaterialTracks + "_mapped"
    materialMappingConfig.unmappedMaterialTracks = outputMaterialTracks + "_unmapped"
    materialMappingConfig.materialWriters = materialMapWriters
    materialMapping = MaterialMapping(materialMappingConfig, loglevel)
    s.addAlgorithm(materialMapping)

    # Add the mapped material tracks writer
    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=materialMappingConfig.mappedMaterialTracks,
            filePath=outputFileBase + "_mapped.root",
            storeSurface=True,
            storeVolume=False,
        )
    )

    # Add the unmapped material tracks writer
    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=materialMappingConfig.unmappedMaterialTracks,
            filePath=outputFileBase + "_unmapped.root",
            storeSurface=True,
            storeVolume=False,
        )
    )

    return s


if "__main__" == __name__:
    p = argparse.ArgumentParser()

    p.add_argument(
        "--geo-dir",
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help="Location of the geometry",
    )
    p.add_argument(
        "-n", "--events", type=int, default=1000, help="Number of events to process"
    )
    p.add_argument(
        "-i", "--input", type=str, default="geant4_material_tracks.root", help="Input file with material tracks"
    )

    p.add_argument(
        "-o", "--output", type=str, default="material", help="Output file (core) name"
    )

    p.add_argument(
        "--matconfig", type=str, default="geometry-map.json", help="Material configuration file"
    )

    p.add_argument(
        "--tree-name",
        type=str,
        default="material_tracks",
        help="Input material track tree name",
    )

    p.add_argument(
        "--material_tracks-name",
        type=str,
        default="material_tracks",
        help="Input material track collection name",
    )

    p.add_argument(
        "--geomGeneration",
        type=str,
        default="gen1",
        choices=["gen1", "gen3"],
        help="Geometry generation to build (gen3 ignores --matconfig, see below).",
    )
    p.add_argument(
        "--root-file",
        type=str,
        default="o2sim_geometry.root",
        help="gen3 only: TGeo .root file inside --geo-dir the builder reads",
    )

    args = p.parse_args()
    logLevel = logging.INFO

    print("Loading geometry from ... ", str(args.geo_dir))

    if args.geomGeneration == "gen3":
        # Gen3 (Blueprint): the C++ builder reads a TGeo .root file directly and
        # attaches ProtoSurfaceMaterial to the receiver surfaces itself when
        # called with withMaterial=True. There is no IMaterialDecorator hook at
        # construction, so --matconfig does not apply here.
        from AliceActsPythonBindings import buildALICE3GeometryGen3

        rootFile = pathlib.Path(args.geo_dir) / args.root_file
        trackingGeometry = buildALICE3GeometryGen3(
            str(rootFile), logLevel, "", True  # graphviz="", withMaterial=True
        )
    else:
        sys.path.insert(0, str(pathlib.Path(args.geo_dir).resolve()))
        import buildALICE3Geometry

        matDeco = None
        if args.matconfig != "":
            matDeco = acts.IMaterialDecorator.fromFile(args.matconfig)

        detector = buildALICE3Geometry.buildALICE3Geometry(
            args.geo_dir, False, False, logLevel, matDeco
        )
        trackingGeometry = detector.trackingGeometry()

    materialSurfaces = trackingGeometry.extractMaterialSurfaces()

    runMaterialMapping(
        materialSurfaces,
        inputFile=Path(args.input),
        nEvents=args.events,
        outputFileBase=args.output,
        outputMapFormats=["json", "root"],
        loglevel=logLevel,
        outputMaterialTracks=args.material_tracks_name,
        treeName=args.tree_name,
    ).run()
