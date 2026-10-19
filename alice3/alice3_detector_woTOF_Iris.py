#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse

from acts.examples.tgeo import TGeoDetector, Interval

from zipfile import ZipFile

import acts
from acts import UnitConstants as u


def buildALICE3Geometry(
    geo_dir: Path,
    material: bool = False,
    jsonconfig: bool = False,
    logLevel=acts.logging.DEBUG,
    matDeco=None
):

    logger = acts.logging.getLogger("buildALICE3Geometry")

    if material:
        file = geo_dir / "old_geom/iris4/material-map.json"
        zip_file = geo_dir / "old_geom/iris4/material-map-iris4.json.zip"

        # check if JSON exists
        if not file.exists():
            # check if ZIP exists
            if zip_file.exists():
                # unzip
                with ZipFile(zip_file, 'r') as zip_ref:
                    zip_ref.extractall(file.parent)
                    print(f"Extracted {file} from {zip_file}")
            else:

                raise FileNotFoundError(f"Neither {json_file} nor {zip_file} exists!")

        logger.info("Adding material from %s", file.absolute())
        matDeco = acts.IMaterialDecorator.fromFile(
            file,
            level=acts.logging.Level(
                min(acts.logging.INFO.value, logLevel.value)),
        )

    tgeo_fileName = geo_dir / \
        "old_geom/iris4/o2sim_geometry.root"

    if jsonconfig:
        jsonFile = geo_dir / "acts/bin/tgeo-config.json"
        logger.info("Create geometry from %s", jsonFile.absolute())
        return TGeoDetector.create(
            jsonFile=str(jsonFile),
            fileName=str(tgeo_fileName),
            surfaceLogLevel=logLevel,
            layerLogLevel=logLevel,
            volumeLogLevel=logLevel,
            mdecorator=matDeco,
        )

    Volume = TGeoDetector.Config.Volume
    LayerTriplet = TGeoDetector.Config.LayerTriplet
    equidistant = TGeoDetector.Config.BinningType.equidistant
    arbitrary = TGeoDetector.Config.BinningType.arbitrary

    return TGeoDetector(
        fileName=str(tgeo_fileName),
        materialDecorator=matDeco,
        buildBeamPipe=True,
        unitScalor=10.0,  # explicit units
        beamPipeRadius=4.8 * u.mm,
        beamPipeHalflengthZ=1000.0 * u.mm,  # 500.0 * u.mm,
        beamPipeLayerThickness=0.25 * u.mm,  # 0.25 * u.mm,
        beamPipeEnvelopeR=0.01 * u.mm,
        layerEnvelopeR=0.01 * u.mm,
        surfaceLogLevel=logLevel,
        layerLogLevel=logLevel,
        volumeLogLevel=logLevel,
        volumes=[
            Volume(
                name="InnerPixels",
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                layers=LayerTriplet(True),
                subVolumeName=LayerTriplet(
                    negative="PETALCASE*", central="PETALCASE*", positive="PETALCASE*"),
                sensitiveNames=LayerTriplet(
                    negative=["PETALCASE*_DISK*_TRKSensor*"], central=["PETALCASE*_LAYER*_TRKSensor*"], positive=["PETALCASE*_DISK*_TRKSensor*"]),
                sensitiveAxes=LayerTriplet("XYZ"),
                rRange=LayerTriplet(negative=(3 * u.mm, 45 * u.mm),
                                    central=(3 * u.mm, 36 * u.mm),
                                    positive=(3 * u.mm, 45 * u.mm)),

                zRange=LayerTriplet(
                    negative=(-400 * u.mm, -250 * u.mm),
                    central=(-250 * u.mm, 250 * u.mm),
                    positive=(250 * u.mm, 400 * u.mm),
                ),
                splitTolR=LayerTriplet(
                    negative=-1.0, central=3 * u.mm, positive=-1.0),
                splitTolZ=LayerTriplet(
                    negative=10 * u.mm, central=-1.0, positive=10 * u.mm),
                binning0=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                binning1=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                cylinderDiscSplit=False,
                cylinderNZSegments=0,
                cylinderNPhiSegments=0,
                discNRSegments=0,
                discNPhiSegments=0,
                itkModuleSplit=False,
                barrelMap={},
                discMap={},
            ),

            Volume(
                name="OuterPixels",
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                layers=LayerTriplet(True),
                subVolumeName=LayerTriplet(
                    negative="FT3*", central="TRK*", positive="FT3*"),
                sensitiveNames=LayerTriplet(
                    negative=["FT3Sensor*"], central=["TRKSensor*"], positive=["FT3Sensor*"]),
                sensitiveAxes=LayerTriplet("XYZ"),
                rRange=LayerTriplet(negative=(50 * u.mm, 500 * u.mm),
                                    central=(50 * u.mm, 440 * u.mm),
                                    positive=(50 * u.mm, 500 * u.mm)),
                zRange=LayerTriplet(
                    negative=(-1300 * u.mm, -700 * u.mm),
                    central=(-700 * u.mm, 700 * u.mm),
                    positive=(700 * u.mm, 1300 * u.mm),
                ),
                splitTolR=LayerTriplet(
                    negative=-1.0, central=1 * u.mm, positive=-1.0),
                splitTolZ=LayerTriplet(
                    negative=5 * u.mm, central=-1.0, positive=5 * u.mm
                ),
                binning0=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                binning1=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                cylinderDiscSplit=True,
                cylinderNZSegments=6,
                cylinderNPhiSegments=32,
                discNRSegments=6,
                discNPhiSegments=32,
                itkModuleSplit=False,
                barrelMap={},
                discMap={},
            ),
            Volume(
                name="OuterMostBarrel",
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                layers=LayerTriplet(
                    positive=False, central=True, negative=False),
                subVolumeName=LayerTriplet("TRK*"),
                sensitiveNames=LayerTriplet(["TRKSensor*"]),
                sensitiveAxes=LayerTriplet("XYZ"),
                rRange=LayerTriplet(
                    central=(445 * u.mm, 1500 * u.mm),
                ),
                zRange=LayerTriplet(
                    central=(-1400 * u.mm, 1400 * u.mm),
                ),
                splitTolR=LayerTriplet(
                    central=5 * u.mm,
                ),
                splitTolZ=LayerTriplet(-1.0),
                binning0=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                binning1=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                cylinderDiscSplit=True,
                cylinderNZSegments=6,
                cylinderNPhiSegments=32,
                discNRSegments=0,
                discNPhiSegments=0,
                itkModuleSplit=False,
                barrelMap={},
                discMap={},
            ),
            Volume(
                name="OuterMostEndcap",
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                layers=LayerTriplet(
                    positive=True, central=False, negative=True),
                subVolumeName=LayerTriplet("FT3*"),
                sensitiveNames=LayerTriplet(["FT3Sensor*"]),
                sensitiveAxes=LayerTriplet("XYZ"),
                rRange=LayerTriplet(
                    negative=(50 * u.mm, 1500 * u.mm),
                    positive=(50 * u.mm, 1500 * u.mm),
                ),
                zRange=LayerTriplet(
                    negative=(-4200 * u.mm, -1350 * u.mm),
                    positive=(1350 * u.mm, 4200 * u.mm),
                ),
                splitTolR=LayerTriplet(-1.0),
                splitTolZ=LayerTriplet(negative=5 * u.mm, positive=5 * u.mm),
                binning0=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                binning1=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                cylinderDiscSplit=True,
                cylinderNZSegments=0,
                cylinderNPhiSegments=0,
                discNRSegments=6,
                discNPhiSegments=32,
                itkModuleSplit=False,
                barrelMap={},
                discMap={},
            ),
        ],
    )


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Example script to construct the ALICE3 geometry and write it out to CSV and OBJ formats"
    )
    p.add_argument(
        "geo_dir",
        help="Input directory containing the ALICE3 standalone geometry.",
    )
    p.add_argument(
        "--output-dir",
        default=Path.cwd(),
        type=Path,
        help="Directory to write outputs to",
    )
    p.add_argument(
        "--output-csv", action="store_true", help="Write geometry in CSV format."
    )
    p.add_argument(
        "--output-obj", action="store_true", help="Write geometry in OBJ format."
    )
    p.add_argument(
        "--output-json",
        action="store_true",
        help="Write geometry and material in JSON format.",
    )
    p.add_argument(
        "--no-material", action="store_true", help="Decorate material to the geometry"
    )

    args = p.parse_args()
    args.output_dir.mkdir(exist_ok=True, parents=True)

    geo_example_dir = Path(args.geo_dir)
    assert geo_example_dir.exists(), "Detector example input directory missing"

    detector = buildALICE3Geometry(
        geo_example_dir,
        material=not args.no_material,
    )

    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    runALICE3(
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        outputDir=args.output_dir,
        outputCsv=args.output_csv,
        outputObj=args.output_obj,
        outputJson=args.output_json,
    )
