#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse

from acts.examples.tgeo import TGeoDetector, Interval

from acts.examples import (
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
import acts.examples

from acts.json import MaterialMapJsonConverter
from acts import UnitConstants as u


def buildALICE3Geometry(
    geo_dir: Path,
    material: bool = False,
    jsonconfig: bool = False,
    logLevel=acts.logging.DEBUG,
    matDeco=None
):

    logger = acts.logging.getLogger("buildALICE3Geometry")

#    matDeco = None
#    matDeco = acts.IMaterialDecorator.fromFile("geometry-map.json")
    if material:
        #        file = geo_dir / "material-maps.root"
        file = geo_dir / "geom/geom_dec25/material-map.json"
        logger.info("Adding material from %s", file.absolute())
        matDeco = acts.IMaterialDecorator.fromFile(
            file,
            level=acts.logging.Level(
                min(acts.logging.INFO.value, logLevel.value)),
        )

    tgeo_fileName = geo_dir / "geom/geom_dec25/o2sim_geometry.root"
    logger.info("tgeo_filename from %s", tgeo_fileName.absolute())

    if jsonconfig:
        jsonFile = geo_dir / "tgeo-config.json"
        logger.info("Create geometry from %s", jsonFile.absolute())
        # return TGeoDetector.create(
        return TGeoDetector(
            jsonFile=str(jsonFile),
            fileName=str(tgeo_fileName),
            surfaceLogLevel=logLevel,
            layerLogLevel=logLevel,
            volumeLogLevel=logLevel,
            materialDecorator=matDeco,
        )

    Volume = TGeoDetector.Config.Volume
    LayerTriplet = TGeoDetector.Config.LayerTriplet
    equidistant = TGeoDetector.Config.BinningType.equidistant
    arbitrary = TGeoDetector.Config.BinningType.arbitrary

    # return TGeoDetector.create(
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
                binTolerancePhi=(0.25 * u.mm, 0.25 * u.mm),
                layers=LayerTriplet(True),
                subVolumeName=LayerTriplet(
                    negative="FT3*", central="TRK*", positive="FT3*"),
                sensitiveNames=LayerTriplet(
                    negative=["FT3Sensor*"], central=["TRKSensor*"], positive=["FT3Sensor*"]),
                sensitiveAxes=LayerTriplet("XYZ"),
                rRange=LayerTriplet(negative=(5 * u.mm, 45 * u.mm),
                                    central=(5.1 * u.mm, 36 * u.mm),
                                    positive=(5 * u.mm, 45 * u.mm)),

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
                cylinderNZSegments=6,
                cylinderNPhiSegments=32,
                discNRSegments=1,
                discNPhiSegments=32,
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
                    negative="FT3*", central="*TRK*", positive="FT3*"),
                sensitiveNames=LayerTriplet(
                    negative=["FT3Sensor*"], central=["*TRKSensor*"], positive=["FT3Sensor*"]),
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
                cylinderDiscSplit=False,
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
                cylinderDiscSplit=False,
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
                cylinderDiscSplit=False,
                cylinderNZSegments=0,
                cylinderNPhiSegments=0,
                discNRSegments=6,
                discNPhiSegments=32,
                itkModuleSplit=False,
                barrelMap={},
                discMap={},
            ),


            Volume(
                name="oTOF",
                binToleranceR=(1 * u.mm, 1 * u.mm),
                binToleranceZ=(1 * u.mm, 1 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                layers=LayerTriplet(
                    positive=False, central=True, negative=False),
                subVolumeName=LayerTriplet("OTOFLayer*"),
                sensitiveNames=LayerTriplet(["OTOFSensor*"]),
                sensitiveAxes=LayerTriplet("XYZ"),
                rRange=LayerTriplet(
                    central=(845 * u.mm, 1205 * u.mm),
                ),
                zRange=LayerTriplet(
                    central=(-3400 * u.mm, 3400 * u.mm),
                ),
                splitTolR=LayerTriplet(
                    central=1 * u.mm,
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
                cylinderDiscSplit=False,
                cylinderNZSegments=4,
                cylinderNPhiSegments=32,
                discNRSegments=0,
                discNPhiSegments=0,
                itkModuleSplit=False,
                barrelMap={},
                discMap={},
            ),


            Volume(
                name="EndcapsTOF",
                binToleranceR=(1 * u.mm, 1 * u.mm),
                binToleranceZ=(1 * u.mm, 1 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                layers=LayerTriplet(
                    positive=True, central=False, negative=True),
                # subVolumeName=LayerTriplet("BTOFLayer*"),
                # sensitiveNames=LayerTriplet(["BTOFSensor*"]),
                subVolumeName=LayerTriplet(
                    negative="BTOFLayer*", positive="FTOFLayer*"),
                sensitiveNames=LayerTriplet(
                    negative=["BTOFSensor*"], positive=["FTOFSensor*"]),

                sensitiveAxes=LayerTriplet("XYZ"),
                rRange=LayerTriplet(
                    central=(20 * u.mm, 1010 * u.mm),
                ),
                zRange=LayerTriplet(
                    negative=(-3750 * u.mm, -3650 * u.mm),
                    positive=(3650 * u.mm, 3750 * u.mm),
                ),
                splitTolR=LayerTriplet(
                    central=1 * u.mm,
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
                cylinderDiscSplit=False,
                cylinderNZSegments=4,
                cylinderNPhiSegments=32,
                discNRSegments=0,
                discNPhiSegments=0,
                itkModuleSplit=False,
                barrelMap={},
                discMap={},
            ),

        ],
    )
