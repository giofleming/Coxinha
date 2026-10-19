#!/usr/bin/env python3
"""
add_disk_eos_cards.py

Add end-of-stave (EOS) copper cards at the OUTER-RADIUS ends of the ALICE 3 FT3
disk staves (analogous to the barrel TRKStave*_EOSCard) and write a new GDML.

Each disk stave (FT3_Stave_*, an extruded-triangular-prism "xtru" subtraction)
gets one copper box card at its outer-radius tip, with:

  - length (default 12 cm) along the GLOBAL z axis (beam),
  - width  (default  4 cm) along the stave width (card centred on the width),
  - thickness (default 4% X/X0 of copper = 0.5744 mm), along the stave axis, so
    the card is "rotated 90 deg" wrt the stave: its long axis points along z.

Card z-placement:
  The card is NOT centred symmetrically on the stave. Instead its IP-facing edge
  (along z) is aligned with the FRONT of the sensor (the stave face toward z=0),
  and the card body extends AWAY from the interaction point (downstream).

Because a displaced 12 cm card no longer fits inside the ~12.5 cm-thick disk
mother tube, the mother ENVELOPE is enlarged downstream only: its solid is
replaced by a union(original tube + short downstream extension tube). This keeps
the mother placement and ALL existing daughters (sensors) untouched, does not
grow toward the barrel (IP side unchanged), and contains the cards.

Approach: minimal XML edit (stdlib xml.etree), no ROOT. The rotation convention is
auto-detected (GDML/Geant4 use the passive/inverse convention) by requiring
every stave to sit inside its mother tube.

Usage:
  python3 add_disk_eos_cards.py [-i in.gdml] [-o out.gdml]
      [--thickness CM] [--length CM] [--width CM]
"""
import argparse
import numpy as np
import xml.etree.ElementTree as etree   # stdlib only, no external dependency

# ----------------------------------------------------------------- constants
X0_CU_CM = 1.436          # copper radiation length [cm] (PDG 12.86 g/cm2 / 8.96)
FRAC_X0 = 0.04            # default thickness = 4% X/X0
CARD_MATERIAL = "TRK_COPPER"
CARD_BOX = "FT3_EOSdiskCard_box"
CARD_VOL = "FT3_EOSdiskCard_Cu"
DISK_MOTHERS = ["FT3Layer0_0", "FT3Layer0_3", "FT3Layer1_0", "FT3Layer1_3"]
Z_MARGIN_CM = 0.20        # small downstream clearance when sizing the envelope


def localname(e):
    # tag is "{ns}local" or "local"
    tag = e.tag
    return tag.rsplit("}", 1)[-1] if isinstance(tag, str) else tag


def rot_matrix(rx, ry, rz, unit, inverse=True):
    if unit == "deg":
        rx, ry, rz = np.radians([rx, ry, rz])
    cx, sx = np.cos(rx), np.sin(rx)
    cy, sy = np.cos(ry), np.sin(ry)
    cz, sz = np.cos(rz), np.sin(rz)
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    R = Rz @ Ry @ Rx
    return R.T if inverse else R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default="o2sim_geometry_noDiskCards.gdml")
    ap.add_argument("-o", "--output",
                    default="o2sim_geometry_with_disk_EOScards.gdml")
    ap.add_argument("--thicknessX", type=float, default=FRAC_X0,
                    help="card thickness in X/X0 units (default %(default)s = 4%% X0)")
    ap.add_argument("--length", type=float, default=12.0,
                    help="card length along beam-z [cm] (default 12)")
    ap.add_argument("--width", type=float, default=4.0,
                    help="card width along stave width [cm] (default 4)")
    args = ap.parse_args()

    t_cm = args.thicknessX * X0_CU_CM       # X/X0 -> cm via copper X0
    L_cm, W_cm = args.length, args.width
    print("=" * 70)
    print("Copper radiation length / card thickness")
    print("  TRK_COPPER in GDML: Z=29, A=63.546 g/mol, rho=8.96 g/cm3 "
          "(= standard copper)")
    print(f"  X0 = 12.86 g/cm2 / 8.96 g/cm3 = {X0_CU_CM:.4f} cm  "
          f"(identical to standard Cu; difference ~ 0)")
    print(f"  card thickness    : {args.thicknessX:.4f} X/X0 -> "
          f"{t_cm:.5f} cm ({t_cm*10:.4f} mm)")
    print(f"  card size (full)  : length(beam-z)={L_cm} cm, width={W_cm} cm")
    print("=" * 70)

    tree = etree.parse(args.input)
    root = tree.getroot()
    NS = root.tag[:root.tag.index("}") + 1] if root.tag.startswith("{") else ""

    define = root.find(NS + "define")
    solids = root.find(NS + "solids")
    structure = root.find(NS + "structure")

    positions = {}
    rotations = {}
    for e in define:
        ln = localname(e)
        if ln == "position":
            f = 1.0 if e.get("unit", "mm") == "cm" else 0.1
            positions[e.get("name")] = np.array(
                [float(e.get("x", 0)), float(e.get("y", 0)),
                 float(e.get("z", 0))]) * f
        elif ln == "rotation":
            rotations[e.get("name")] = (
                float(e.get("x", 0)), float(e.get("y", 0)),
                float(e.get("z", 0)), e.get("unit", "deg"))

    xtru = {}
    tube = {}    # name -> (element, dict rmin/rmax/z)
    subtr = {}
    for e in solids:
        ln = localname(e)
        nm = e.get("name")
        if ln == "xtru":
            f = 1.0 if e.get("lunit", "mm") == "cm" else 0.1
            xs, ys, zs = [], [], []
            for v in e:
                if localname(v) == "twoDimVertex":
                    xs.append(float(v.get("x")) * f)
                    ys.append(float(v.get("y")) * f)
                elif localname(v) == "section":
                    zs.append(float(v.get("zPosition")) * f)
            xtru[nm] = (min(xs), max(xs), min(ys), max(ys), min(zs), max(zs))
        elif ln == "tube":
            f = 1.0 if e.get("lunit", "mm") == "cm" else 0.1
            # GDML tube z is the FULL length -> store the HALF length
            tube[nm] = (e, dict(rmin=float(e.get("rmin", 0)) * f,
                                rmax=float(e.get("rmax")) * f,
                                z=float(e.get("z")) * f / 2.0))
        elif ln == "subtraction":
            subtr[nm] = e.find(NS + "first").get("ref")

    vol_solid = {}
    vol_elem = {}
    for v in structure:
        if localname(v) not in ("volume", "assembly"):
            continue
        vol_elem[v.get("name")] = v
        sr = v.find(NS + "solidref")
        vol_solid[v.get("name")] = sr.get("ref") if sr is not None else None

    # mother global-z (from the placement in the parent assembly, e.g. TRKV)
    mother_globalz = {}
    for parent in structure:
        if localname(parent) not in ("volume", "assembly"):
            continue
        for p in parent.findall(NS + "physvol"):
            vr = p.find(NS + "volumeref").get("ref")
            if vr in DISK_MOTHERS:
                pr = p.find(NS + "positionref")
                mother_globalz[vr] = positions[pr.get("ref")][2] if pr is not None else 0.0

    def stave_bbox(stavevol):
        s = vol_solid.get(stavevol)
        if s not in subtr:
            return None
        return xtru.get(subtr[s])

    # ---- auto-detect rotation convention -----------------------------------
    def fits_all(inverse):
        for mother in DISK_MOTHERS:
            tb = tube[vol_solid[mother]][1]
            for p in vol_elem[mother].findall(NS + "physvol"):
                bb = stave_bbox(p.find(NS + "volumeref").get("ref"))
                if bb is None:
                    continue
                pos = positions[p.find(NS + "positionref").get("ref")]
                rr = p.find(NS + "rotationref")
                ang = rotations[rr.get("ref")] if rr is not None else (0, 0, 0, "deg")
                R = rot_matrix(*ang, inverse=inverse)
                xmn, xmx, ymn, ymx, zmn, zmx = bb
                c = np.array([[x, y, z] for x in (xmn, xmx)
                              for y in (ymn, ymx) for z in (zmn, zmx)])
                g = (R @ c.T).T + pos
                r = np.hypot(g[:, 0], g[:, 1])
                if not (r.min() >= tb["rmin"] - 0.1 and r.max() <= tb["rmax"] + 0.1
                        and np.abs(g[:, 2]).max() <= tb["z"] + 0.1):
                    return False
        return True

    inverse = True if fits_all(True) else (False if fits_all(False) else True)
    print(f"Rotation convention: {'passive/inverse' if inverse else 'active'} "
          f"(auto-detected)\n")

    # ---- shared copper card box + volume -----------------------------------
    box = etree.SubElement(solids, NS + "box")
    box.set("name", CARD_BOX)
    box.set("x", f"{W_cm}")       # width  -> stave width (local x)
    box.set("y", f"{L_cm}")       # length -> beam z (local y -> mother z)
    box.set("z", f"{t_cm:.6f}")   # thickness -> stave axis (local z)
    box.set("lunit", "cm")

    cardvol = etree.Element(NS + "volume")
    cardvol.set("name", CARD_VOL)
    etree.SubElement(cardvol, NS + "materialref").set("ref", CARD_MATERIAL)
    etree.SubElement(cardvol, NS + "solidref").set("ref", CARD_BOX)
    structure.insert(0, cardvol)

    hx, hy, hz = W_cm / 2, L_cm / 2, t_cm / 2
    box_corners = np.array([[sx * hx, sy * hy, sz * hz]
                            for sx in (-1, 1) for sy in (-1, 1)
                            for sz in (-1, 1)])

    ncards = 0
    nnudged = 0
    report = []
    for mother in DISK_MOTHERS:
        tb_elem, tb = tube[vol_solid[mother]]
        melem = vol_elem[mother]
        gz = mother_globalz.get(mother, 0.0)
        d = 1.0 if gz >= 0 else -1.0            # downstream (away from IP) sign
        H0 = tb["z"]                             # original half-z (local)

        stave_pvs = [p for p in melem.findall(NS + "physvol")
                     if stave_bbox(p.find(NS + "volumeref").get("ref")) is not None]

        sf_lo, sf_hi = 1e9, -1e9
        card_zlo, card_zhi = 1e9, -1e9
        ds_edge = 0.0                            # most downstream card edge (local z)
        nudges = []
        for i, p in enumerate(stave_pvs):
            vr = p.find(NS + "volumeref").get("ref")
            pr = p.find(NS + "positionref").get("ref")
            rr = p.find(NS + "rotationref")
            rrname = rr.get("ref") if rr is not None else None
            pos = positions[pr]
            ang = rotations[rrname] if rrname else (0, 0, 0, "deg")
            R = rot_matrix(*ang, inverse=inverse)
            xmn, xmx, ymn, ymx, zmn, zmx = stave_bbox(vr)

            # stave footprint in mother frame -> z range and outer-tip radius end
            corners = np.array([[x, y, z] for x in (xmn, xmx)
                                for y in (ymn, ymx) for z in (zmn, zmx)])
            g = (R @ corners.T).T + pos
            sz_lo, sz_hi = g[:, 2].min(), g[:, 2].max()
            # sensor FRONT (toward IP) = stave face in the -d direction
            sensor_front = sz_lo if d > 0 else sz_hi
            sf_lo, sf_hi = min(sf_lo, sensor_front), max(sf_hi, sensor_front)
            # card centre z: IP-edge at sensor_front, body extends downstream
            centre_z = sensor_front + d * hy
            ds_edge = (max(ds_edge, sensor_front + d * L_cm) if d > 0
                       else min(ds_edge, sensor_front + d * L_cm))

            # radial: outer tip = stave end (local z zmn/zmx) with larger r
            def end_r(ze):
                gg = R @ np.array([0.0, 0.5 * (ymn + ymx), ze]) + pos
                return np.hypot(gg[0], gg[1])
            z_out = zmx if end_r(zmx) >= end_r(zmn) else zmn
            sgn = 1.0 if z_out == zmx else -1.0
            z_thick = z_out + sgn * hz          # just outside the tip (thickness)

            # y_c so that the card centre lands at mother-z = centre_z
            if abs(R[2, 1]) > 1e-6:
                y_c = (centre_z - R[2, 2] * z_thick - pos[2]) / R[2, 1]
            else:
                y_c = 0.5 * (ymn + ymx)
            centre = R @ np.array([0.0, y_c, z_thick]) + pos

            # radial clamp within mother rmax
            wc = (R @ box_corners.T).T + centre
            rmx = np.hypot(wc[:, 0], wc[:, 1]).max()
            if rmx > tb["rmax"]:
                rc = np.hypot(centre[0], centre[1])
                if rc > 1e-6:
                    centre = centre - np.array([centre[0], centre[1], 0.0]) / rc * (rmx - tb["rmax"])
                    nudges.append(rmx - tb["rmax"])
                    nnudged += 1
            wc = (R @ box_corners.T).T + centre
            card_zlo = min(card_zlo, wc[:, 2].min())
            card_zhi = max(card_zhi, wc[:, 2].max())

            posname = f"{CARD_VOL}_{mother}_{i}_pos"
            pe = etree.SubElement(define, NS + "position")
            pe.set("name", posname)
            for a, val in zip("xyz", centre):
                pe.set(a, f"{val:.10f}")
            pe.set("unit", "cm")

            pv = etree.SubElement(melem, NS + "physvol")
            pv.set("name", f"{CARD_VOL}_{mother}_pv{i}")
            etree.SubElement(pv, NS + "volumeref").set("ref", CARD_VOL)
            etree.SubElement(pv, NS + "positionref").set("ref", posname)
            if rrname is not None:
                etree.SubElement(pv, NS + "rotationref").set("ref", rrname)
            ncards += 1

        # ---- enlarge the mother envelope DOWNSTREAM only (union) -----------
        # need coverage out to |ds_edge| (+margin); original covers +-H0.
        need = abs(ds_edge) + Z_MARGIN_CM
        new_half_total = None
        if need > H0:
            ext_half = 0.5 * (need - H0)
            ext_centre = d * (H0 + ext_half)     # local z of extension centre
            ext_tube_name = f"{mother}_EOSext_tube"
            ext_pos_name = f"{mother}_EOSext_pos"
            et = etree.SubElement(solids, NS + "tube")
            et.set("name", ext_tube_name)
            et.set("rmin", f"{tb['rmin']}")
            et.set("rmax", f"{tb['rmax']}")
            et.set("z", f"{2.0 * ext_half:.6f}")   # GDML tube z is FULL length
            et.set("startphi", "0")
            et.set("deltaphi", "360")
            et.set("aunit", "deg")
            et.set("lunit", "cm")
            ep = etree.SubElement(define, NS + "position")
            ep.set("name", ext_pos_name)
            ep.set("x", "0")
            ep.set("y", "0")
            ep.set("z", f"{ext_centre:.6f}")
            ep.set("unit", "cm")
            uni = etree.SubElement(solids, NS + "union")
            uni.set("name", f"{mother}_env")
            etree.SubElement(uni, NS + "first").set("ref", vol_solid[mother])
            etree.SubElement(uni, NS + "second").set("ref", ext_tube_name)
            etree.SubElement(uni, NS + "positionref").set("ref", ext_pos_name)
            melem.find(NS + "solidref").set("ref", f"{mother}_env")
            new_half_total = H0 + 2 * ext_half   # downstream reach from centre

        report.append((mother, gz, d, len(stave_pvs), sf_lo, sf_hi,
                       card_zlo, card_zhi, H0, need, new_half_total,
                       len(nudges), max(nudges) if nudges else 0.0))

    print("Per-disk summary (lengths in cm, local mother z):")
    print(f"  {'mother':12s} {'globalz':>8s} {'d':>3s} {'#cards':>6s} "
          f"{'sensFrontZ':>16s} {'cardZ':>16s} {'origHalf':>8s} "
          f"{'needHalf':>8s} {'#nudge':>6s}")
    for (m, gz, d, n, sflo, sfhi, czlo, czhi, H0, need, nht, nn, mn) in report:
        print(f"  {m:12s} {gz:8.2f} {int(d):3d} {n:6d} "
              f"[{sflo:6.2f},{sfhi:6.2f}] [{czlo:6.2f},{czhi:6.2f}] "
              f"{H0:8.2f} {need:8.2f} {nn:6d}")

    tree.write(args.output, xml_declaration=True, encoding="UTF-8")
    print(f"\nAdded {ncards} disk EOS cards ({nnudged} radially nudged within "
          f"rmax). Mother envelopes extended downstream via union.")
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
