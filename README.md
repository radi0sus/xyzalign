# xyzalign
A Python 3 script for aligning, rotating and translating atomic coordinates in xyz files. The main purpose of the script is to align ligand atoms along the axes of the cartesian coordinate system, which is useful for computational chemistry (especially in the field of coordination chemistry). Single atoms or the centroid of several atoms (useful for cyclopentadienyl ligands for example) can be aligned. The modified xyz file will be saved. Alternatively, the modified coordinates will be printed to the console. 

## External modules
`pandas`, `numpy`

## Quick start
 Start the script with:
```console
python3 xyzalign.py filename.xyz
```
to open the XYZ. It will calculate the centroid or geometric center of all atoms and make this the origin at x=0, y=0 and z=0. All other coordinates will be recalculated with respect to the new origin. The file will be saved as `filename-mod.xyz`.

```console
python3 xyzalign.py filename.xyz -o 1 -x 2 -y 3 -z 4
```
Open `filename.xyz`, set the origin (x=y=z=0 or 0,0,0) to atom 1 or the first atom in the xyz file (`-o 1`), align atom 2 to x (1,0,0) (`-x 2`) , align atom 3 to y (0,1,0) (`-y 3`) and align atom 4 to z (0,0,1) (`-z 4`). Since the atoms are probably not perfectly aligned to x, y and z the script is trying to align them as close as possible to the three axes (x has some priority. The resultant xyz file will be saved as `filename-mod.xyz`.

```console
python3 xyzalign.py ferrocene.xyz -o 1 -z 2 3 4 5 6
```
Open `ferrocene.xyz`, set the origin (x=y=z=0 or 0,0,0) to atom 1 (`-o 1`), which is Fe and align the centroid of the atoms 2, 3, 4, 5, 6 (carbon atoms of one cyclopentadienyl  ring) to z. Save the modified xyz file as `ferrocene-mod.xyz`.
