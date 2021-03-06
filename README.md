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
Open `filename.xyz`, set the origin (x=y=z=0 or 0,0,0) to atom 1 or the first atom in the xyz file (`-o 1`), align atom 2 to x (1,0,0) (`-x 2`) , align atom 3 to y (0,1,0) (`-y 3`) and align atom 4 to z (0,0,1) (`-z 4`). Since the atoms are probably not perfectly aligned to x, y and z the script is trying to align them as close as possible to the three axes (x has some priority). The resultant xyz file will be saved as `filename-mod.xyz`.

```console
python3 xyzalign.py ferrocene.xyz -o 1 -z 2 3 4 5 6
```
Open `ferrocene.xyz`, set the origin (x=y=z=0 or 0,0,0) to atom 1 (`-o 1`), which is Fe and align the centroid of the atoms 2, 3, 4, 5, 6 (carbon atoms of one cyclopentadienyl  ring) to the z axis (0,0,1). Save the modified xyz file as `ferrocene-mod.xyz`.

```console
python3 xyzalign.py fe2s2.xyz -o 1 2 3 4 -r 0 0 90 -s
```
Open `fe2s2.xyz`, set the origin (x=y=z=0 or 0,0,0) to the centroid of the first 4 atoms (`-o 1 2 3 4`). Counterclockwise rotate the molecule around the z-axis by 90° (`-r 0 0 90`). Do not save the xyz file, print the content of the xyz file to the console (`-s`).

## Command-line options
- `filename` , required: filename, e.g. `my_xyz.xyz`, first two lines will be ignored, file format must be `element x y z`, cartesian coordinates, (units in Å)
- `-o` `atom(s)`, optional:  define the origin of the molecule by one or more atoms, e.g. `-o 1` or `-o 1 2 3`. If no atom is defined (i.e. the `-o` option is not given), the centroid of all atoms in the xy file will be calculated. The origin will have the coordinates x=0, y=0 and z=0. All other coordinates will be recalculated with respect to the new origin.
-  `-x` `atom(s)`, optional: define the atom or atoms that should be aligned with the x-axis (1,0,0), e.g. `-x 1` or `-x 1 2 3`. If more than one atom is given, the centroid of these atoms is aligned to the x-axis.
-  `-y` `atom(s)`, optional: define the atom or atoms that should be aligned with the y-axis (0,1,0), e.g. `-y 1` or `-y 1 2 3`. If more than one atom is given, the centroid of these atoms is aligned to the y-axis.
-  `-z` `atom(s)`, optional: define the atom or atoms that should be aligned with the z-axis (0,0,1), e.g. `-z 1` or `-z 1 2 3`. If more than one atom is given, the centroid of these atoms is aligned to the z-axis.
-  `-r` `angle-x angle-y angle-z`, optional: counterclockwise rotate the coordinates around the x-, y- and z-axes, e.g. `-r 45.1 0 90` or `-r 0 0 90`. Angles are in degrees.
-  `-t` `x y z`, optional: translate coordinates in x, y and z, e.g. `-t 2 1 0` (translate coordinates by 2 Å or units of the coordinates in x, 1 Å in y and 0 Å in z)
-  `-m` `x1 y1 z1 x2 y2 z2 x3 y3 z3`, optional: rotation / transformation matrix, e.g. `-m -1 0 0 0 -1 0 0 0 -1` for inversion of the coordinates.
-  `-s`, optional: print the content of the altered xyz file to the console, do not save the file.
-  `-v`, optional: verbose mode. Print more information, like the origin of the input xyz file or rotation matrices.

## Remarks
- Only the standard XYZ file format is supported. 
- If the script is opened with no options, the molecule will be centered and the modified coordinates will be saved as `-mod.xyz`.
- If saving is not necessary or not possible, the `-s` option should be used.
- A perfect alignment of atoms of 'real world' molecules to the three axes of a cartesian coordinate system is probably impossible. The script tries to align the atoms as close as possible to the axes.
- The selected atoms are preferably placed along the positive xyz-axes of the right-handed coordinate system. However, this is not always possible. Depending on the molecular geometry, the atoms are then sometimes also located on negative axes after the transformation. 
- The rotation of positive rotation angles is counterclockwise about the axis of rotation.
- To makes calculations easier, the origin is always at 0,0,0. To transfer the coordinates back to the starting origin, use the `-v` option, note the original origin and translate the coordinates back with the `-t` option.
- If counting of atoms should start at zero, comment the line `xyz_df.index +=1` in the script.

## Examples
All images were created with [xyz2tab](https://github.com/radi0sus/xyz2tab).
### Example 1:
Not aligned molecule:
<p align="center">
<img width="400" alt="nipor1" src="/examples/nipor1.png">
</p>
Origin set to atom 1 (metal atom in the center, red),...:

```console
python3 xyzalign.py filename.xyz -o 1 
```
<p align="center">
<img width="400" alt="nipor2" src="/examples/nipor2.png">
</p>
...atoms 2 and 3 (nitrogen atoms, blue) aligned with the x- and y-axes,...:

```console
python3 xyzalign.py filename.xyz -o 1 -x 2 -y 3
```

<p align="center">
<img width="400" alt="nipor3" src="/examples/nipor3.png">
</p>
...45° counterclockwise rotation around the z-axis,...:

```console
python3 xyzalign.py filename.xyz -o 1 -x 2 -y 3 -r 0 0 45
```

<p align="center">
<img width="400" alt="nipor4" src="/examples/nipor4.png">
</p>
...1.6 Å translation in z.:

```console
python3 xyzalign.py filename.xyz -o 1 -x 2 -y 3 -r 0 0 45 -t 0 0 1.6
```

<p align="center">
<img width="400" alt="nipor4" src="/examples/nipor5.png">
</p>

<p align="center">
<img width="800" alt="nipor4" src="/examples/show-use1.gif">
</p>

### Example 2:
Not aligned molecule:
<p align="center">
<img width="400" alt="nipor1" src="/examples/fer1.png">
</p>
Set origin to atom 1 (metal atom in the center, red) and align the centroid of atoms 2, 3, 4, 5, and 6 (carbon atoms, black) to the z-axis:

```console
python3 xyzalign.py filename.xyz -o 1 -z 2 3 4 5 6
```

<p align="center">
<img width="400" alt="nipor4" src="/examples/fer2.png">
</p>
