from pathlib import Path

p = Path("o2sim_geometry.gdml")
text = p.read_text(encoding="utf-8")

# 1. Shrink the cave volume to the tracker volume
# old = '<box name="cave" x="100" y="100" z="800" lunit="cm"/>'
old = '<box name="cave" x="4000" y="4000" z="6000" lunit="cm"/>'
new = '<box name="cave" x="350" y="350" z="1100" lunit="cm"/>'
assert old in text, f"Could not find: {old}"
text = text.replace(old, new)
print(f"cave solid: {old} --> {new}")

# 2. Add new barrel tube solid before </solids>
old_solids_close = '</solids>'
new_tube = '  <tube name="shBarrelTRK" rmin="0" rmax="135" z="1090" startphi="0" deltaphi="360" aunit="deg" lunit="cm"/>\n</solids>'
assert old_solids_close in text, f"Could not find: {old_solids_close}"
assert 'shBarrelTRK' not in text, "shBarrelTRK already exists in the file!"
text = text.replace(old_solids_close, new_tube, 1)  # replace only first occurrence
print(f"Added shBarrelTRK tube solid before </solids>")

# 3. Swap solidref in barrel volume only
old_solidref = '<solidref ref="shCaveTR"/>'
new_solidref = '<solidref ref="shBarrelTRK"/>'
assert old_solidref in text, f"Could not find: {old_solidref}"
assert text.count(old_solidref) == 1, f"Found multiple occurrences of {old_solidref} - unsafe to replace!"
text = text.replace(old_solidref, new_solidref)
print(f"barrel solidref: {old_solidref} --> {new_solidref}")

# 4. Write back
p.write_text(text, encoding="utf-8")

# 5. Verify
assert 'shBarrelTRK' in text, "shBarrelTRK not found after replacement!"
assert new_solidref in text, "shBarrelTRK solidref not found after replacement!"
print("All replacements verified successfully.")
