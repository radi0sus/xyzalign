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
- To makes calculations easier, the origin is always at 0,0,0. To transfer the coordinates back to the starting origin, use the `-v` option, note the original origin and translate the coordinates back with the `-t` option.
- If counting of atoms should start at zero, comment the line `xyz_df.index +=1` in the script.

## Examples

```console
python3 xyzalign.py test.xyz
```

```
7						             7
test					           test
Ni  1  1  1				     Ni  0  0  0
Cl  3  1  1				     Cl  2  0  0
Br  1  3  1		===>		 Br  0  2  0
I   1  1  3				     I   0  0  2	
He -1  1  1				     He -2  0  0	
Ne  1 -1  1				     Ne  0 -2  0
Ar  1  1 -1				     Ar  0  0 -2
```
