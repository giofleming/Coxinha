import json
import sys

import acts
import acts.examples

from acts.examples import (
    StructureSelector
)

from acts.examples.alignment import (
    AlignmentDecorator,
    GeoIdAlignmentStore,
)


def selectAlignableSurfaces(trackingGeometry, geoId, selector=None):
    """Alignable surfaces of the structure identified by geoId.

    Same selection as acts.examples.StructureSelector (lowest non-zero level
    of geoId: volume, layer or single module), restricted to surfaces backed
    by a detector element (Surface.isAlignable). Only those follow an
    AlignmentContext. This drops the gen3 passive cylinders/discs and the
    material representatives, which carry sensitive IDs but no placement.

    Returns a dict {GeometryIdentifier: Surface}. Pass a prebuilt
    StructureSelector to avoid re-visiting the geometry on repeated calls.
    """
    if selector is None:
        selector = StructureSelector(trackingGeometry)
    return {s.geometryId: s for s in selector.selectSurfaces(geoId) if s.isAlignable}


_SHIFT_KEYS = {"volume", "layer", "sensitive", "shift"}


def addGlobalShifts(sequencer,
                    trackingGeometry,
                    globalShifts,
                    outputDir=None,
                    target=AlignmentDecorator.Target.eReco):
    """Misalign structures by global translations, all in one alignment store.

    globalShifts is a list of dicts, one per structure:
        {'volume': 50, 'layer': 5, 'shift': (tx, ty, tz)}   # mm, global frame
    'layer' / 'sensitive' omitted or 0 select the whole volume / layer (see
    selectAlignableSurfaces). The shifted transform is T' = Translation(shift) * T.
    Entries selecting the same surface add up (translations commute).

    A single AlignmentDecorator carries all the shifts: each decorator replaces
    the geometry context, so one decorator per structure would keep only the last.

    With the default target eReco, simulation and digitization see the nominal
    geometry and reconstruction the shifted one. This is what Geant4 needs, as
    its geometry is built from the GDML and ignores the ACTS context. It is
    equivalent to the real detector being displaced by -shift.

    If outputDir is given, the applied entries are written to misalignment.json.
    """
    if not globalShifts:
        raise ValueError("addGlobalShifts: no entries in globalShifts")

    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
    selector = StructureSelector(trackingGeometry)

    nominal = {}   # geoId -> nominal Transform3
    offsets = {}   # geoId -> summed global shift [tx, ty, tz]
    record = []
    for entry in globalShifts:
        unknown = set(entry) - _SHIFT_KEYS
        if unknown or "volume" not in entry or "shift" not in entry:
            raise ValueError(f"addGlobalShifts: bad entry {entry}, expected keys {sorted(_SHIFT_KEYS)} "
                             "with 'volume' and 'shift' required")
        geoId = acts.GeometryIdentifier(volume=entry["volume"],
                                        layer=entry.get("layer", 0),
                                        sensitive=entry.get("sensitive", 0))
        shift = [float(x) for x in entry["shift"]]
        if len(shift) != 3:
            raise ValueError(f"addGlobalShifts: 'shift' must be (tx, ty, tz), got {entry['shift']}")

        surfaces = selectAlignableSurfaces(trackingGeometry, geoId, selector)
        if not surfaces:
            raise RuntimeError(f"addGlobalShifts: no alignable surfaces selected by {geoId}: "
                               "the misalignment would be a no-op")
        for gid, srf in surfaces.items():
            if gid not in nominal:
                nominal[gid] = srf.localToGlobalTransform(gctx)
                offsets[gid] = [0., 0., 0.]
            offsets[gid] = [o + d for o, d in zip(offsets[gid], shift)]

        print(f"== INFO == misalignment: {geoId} shifted by {shift} mm ({len(surfaces)} surfaces)")
        record.append({"volume": geoId.volume, "layer": geoId.layer, "sensitive": geoId.sensitive,
                       "shift_mm": shift, "nSurfaces": len(surfaces)})

    misaligned = {}
    for gid, T in nominal.items():
        t = T.translation
        misaligned[gid] = acts.Transform3(acts.Vector3(t[0] + offsets[gid][0],
                                                       t[1] + offsets[gid][1],
                                                       t[2] + offsets[gid][2]),
                                          T.rotation)

    decoConfig = AlignmentDecorator.Config()
    decoConfig.target = target
    # one store valid for every event of the job
    decoConfig.iovStores = [((0, sys.maxsize), GeoIdAlignmentStore(misaligned))]
    sequencer.addContextDecorator(AlignmentDecorator(decoConfig, acts.logging.INFO))
    print(f"== INFO == misalignment: {len(misaligned)} surfaces in the alignment store, target {target}")

    if outputDir is not None:
        with open(f"{outputDir}/misalignment.json", "w") as f:
            json.dump({"target": str(target), "frame": "global", "globalShifts": record}, f, indent=2)
