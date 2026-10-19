#!/usr/bin/env python3
import os
# from common import getOpenDataDetectorDirectory
# from acts.examples.odd import getOpenDataDetector
from acts.examples import (
    GenericDetector,
    #    AlignedDetector,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    CsvTrackingGeometryWriter,
    ObjTrackingGeometryWriter,
)

from acts.examples.json import (
    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
)

import acts
from acts.json import MaterialMapJsonConverter

# import alice3
import pathlib
import argparse
import sys


def runGeometry(
    trackingGeometry,
    decorators,
    outputDir,
    events=1,
    outputObj=False,
    outputCsv=False,
    outputJson=True,
    outputRoot=False,
):

    for ievt in range(events):
        eventStore = WhiteBoard(name=f"EventStore#{ievt}", level=acts.logging.INFO)
        ialg = 0
        ithread = 0

        context = AlgorithmContext(ialg, ievt, eventStore, ithread)

        for cdr in decorators:
            r = cdr.decorate(context)
            if r != ProcessCode.SUCCESS:
                raise RuntimeError("Failed to decorate event context")

        if outputCsv:
            writer = CsvTrackingGeometryWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=os.path.join(outputDir, "csv"),
                writePerEvent=True,
            )
            writer.write(context)

        if outputObj:
            writer = ObjTrackingGeometryWriter(
                level=acts.logging.INFO, outputDir=os.path.join(outputDir, "obj")
            )
            writer.write(context, trackingGeometry)

        if outputJson:
            writer = JsonSurfacesWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=os.path.join(outputDir, "json"),
                writePerEvent=True,
                writeSensitive=True,
            )
            writer.write(context)

            jmConverterCfg = MaterialMapJsonConverter.Config(
                processSensitives=True,
                processApproaches=True,
                processRepresenting=True,
                processBoundaries=True,
                processVolumes=True,
                processNonMaterial=True,
                context=context.geoContext,
            )

            jmw = JsonMaterialWriter(
                level=acts.logging.VERBOSE,
                converterCfg=jmConverterCfg,
                fileName=os.path.join(outputDir, "geometry-map"),
                writeFormat=JsonFormat.Json,
            )

            jmw.write(trackingGeometry)


parser = argparse.ArgumentParser(
    description="Geometry to json")

parser.add_argument(
    "--geo-dir",
    help="Location of the geometry",
    type=pathlib.Path,
    default=pathlib.Path.cwd()
)

parser.add_argument(
    "--geomGeneration",
    type=str,
    default="gen1",
    choices=["gen1", "gen3"],
    help="Geometry generation to dump to geometry-map.json",
)

parser.add_argument(
    "--root-file",
    type=str,
    default="o2sim_geometry.root",
    help="gen3 only: TGeo .root file inside --geo-dir the builder reads",
)

args = parser.parse_args()


geo_dir = args.geo_dir
print("Loading geometry from ... ", str(geo_dir))
if args.geomGeneration == "gen1":
    sys.path.insert(0, str(pathlib.Path(geo_dir).resolve()))
    import buildALICE3Geometry


if "__main__" == __name__:
    #    detector, trackingGeometry, decorators = AlignedDetector.create()

    #     detector, trackingGeometry, decorators = alice3.buildALICE3Geometry(
    # #    geo_dir, False, False, acts.logging.VERBOSE)
    # #    geo_dir, False, False, acts.logging.DEBUG)
    #     geo_dir, False, False, acts.logging.INFO)
    #     # detector, trackingGeometry, decorators = GenericDetector.create()
    #     # detector, trackingGeometry, decorators = getOpenDataDetector(getOpenDataDetectorDirectory())

    if args.geomGeneration == "gen3":
        # Gen3: withMaterial=True so the passive layers and the synthetic
        # representative surfaces are present and get the SAME geo-ids as in the
        # mapping step - otherwise geometry-map.json would not match the sur_id
        # in material_gen3_mapped.root and the plots would be empty.
        from AliceActsPythonBindings import buildALICE3GeometryGen3
        trackingGeometry = buildALICE3GeometryGen3(
            str(pathlib.Path(geo_dir) / args.root_file), acts.logging.INFO, "", True)
        decorators = []
    else:
        detector = buildALICE3Geometry.buildALICE3Geometry(
            #    geo_dir, False, False, acts.logging.VERBOSE)
            geo_dir, False, False, acts.logging.INFO)
        trackingGeometry = detector.trackingGeometry()
        decorators = detector.contextDecorators()

    runGeometry(trackingGeometry, decorators, outputDir=os.getcwd())
