# xyzalign
A Python 3 script for aligning, rotating and translating atomic coordinates in xyz files. The main purpose of the script is to align ligand atoms along the axes of the cartesian coordinate system which is useful for computational chemistry (especially in the field of coordination chemistry). Single atoms or the centroid of several atoms (useful for cyclopentadienyl ligands for example) can be aligned. The modified xyz file will be saved. Alternatively, the modified coordinates will be printed to the console. 

## External modules
`pandas`, `numpy`

## Quick start
 Start the script with:
```console
python3 xyzalign.py filename.xyz
```
to open the XYZ. It will calculate the centroid or geometric center of all atoms and make this the origin at x=0, y=0 and z=0. All other coordinates will be recalculated with respect to the new origin. The file will be saved as `filename-mod.xyz`.
