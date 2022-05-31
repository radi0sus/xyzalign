#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#xyzalign

#import sys                                                     #stdout
import argparse                                                 #argument parser
import pandas as pd                                             #pandas data frames
import numpy as np                                              #calculations
import os.path                                                  #filename split 

#options for printig pandas tables
pd.set_option('display.max_rows', None)
pd.set_option('display.float_format', lambda x: '%.4f' % x)

#rotation matrix from two vectors
def rotmat_from_vec(vec1, vec2):
	#normalize
	a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
	v = np.cross(a, b)
	c = np.dot(a, b)
	#to avoid division by zero
	if c == -1:
		c = -0.9999999999999999
	k = 1.0 / (1.0 + c)
	return np.array([[v[0] * v[0] * k + c, v[1] * v[0] * k - v[2], v[2] * v[0] * k + v[1]],
	[v[0] * v[1] * k + v[2], v[1] * v[1] * k + c, v[2] * v[1] * k - v[0]],
	[v[0] * v[2] * k - v[1], v[1] * v[2] * k + v[0], v[2] * v[2] * k + c]])

#apply rotation matrix to coordinates
def align_xyz(vec1,vec2,coord):
	rotmatrix = rotmat_from_vec(vec1,vec2)
	return np.dot(coord,rotmatrix.T)

#get rotation matrix from angles
def rotmat_from_ang(theta, R = np.zeros((3,3))):
	theta = np.array(theta)*np.pi/180.
	cx, cy, cz = np.cos(theta)
	sx, sy, sz = np.sin(theta)
#	R.flat = (cx*cz, -cy*sz+sy*sx*cz,  sy*sz+cy*sx*cy, 
#		cx*sz,  cy*cz+sy*sx*sz,  -sy*cz+cy*sx*sy, 
#		-sx  ,  sy*cx         ,  cy*cx         )
	R.flat = (cy*cz, -cx*sz+sx*sy*cz,  sx*sz+cx*sy*cx, 
			  cy*sz,  cx*cz+sx*sy*sz,  -sx*cz+cx*sy*sx, 
			 -sy   ,  sx*cy         ,  cx*cy          )
	return R

parser = argparse.ArgumentParser(
	prog='xyzalign.py', 
	formatter_class=argparse.RawDescriptionHelpFormatter,
	description=('''\
align, rotate, translate xyz coordinates
output will be saved as filename-mod.xyz
atom number 1 in the xyz-file is 1
'''))

#filename is required
parser.add_argument('filename', 
	help = "filename, xyz; e.g. mymolecule.xyz")

#origin for transformations / rotations
parser.add_argument('-o', '--origin',
	nargs='+',
	default=12345,
	type=int,
	help='define the origin (at 0,0,0) by one or more atoms,\
		if no atom is defined, the centroid of all atoms will be calculated, \
		e.g. -o 1 or -o 1 2 3')

#atom in x direction
parser.add_argument('-x','--x',
	nargs='+',
	type=int,
	help='define the atom or atoms in the x-direction, \
		e.g. -x 1 or -x 2 3 4')

#atom in y direction
parser.add_argument('-y','--y',
	nargs='+',
	type=int,
	help='define the atom or atoms in the y-direction, \
		e.g. -y 1 or -y 2 3 4')

#atom in z direction
parser.add_argument('-z','--z',
	nargs='+',
	type=int,
	help='define the atom or atoms in the z-direction, \
		e.g. -z 1 or -z 2 3 4')

#rotate
parser.add_argument('-r','--rotate',
	nargs=3,
	type=float,
	help='rotation about the x-, y- and z-axis, \
		order is x y z, angles in degrees, \
		e.g. -r 45.11 90 0 or -r 0 -180 0')

#transformation matrix
parser.add_argument('-m','--matrix',
	nargs=9,
	type=float,
	help='use a rotation matrix x1 y1 z1 x2 y2 z2 x3 y3 z3 \
		to align the xyz coordinates, \
		e.g. -m -1 0 0 0 -1 0 0 0 -1')

#translate
parser.add_argument('-t','--translate',
	nargs=3,
	type=float,
	help='translate in x-, y- and z-direction,\
		order is x y z, distances in Å,\
		e.g. -t 2.11 0 3 or -t 0 -1.81 0')

#verbose output
parser.add_argument('-v','--verbose',
	default=0,
	action='store_true',
	help='verbose output')

#print to stdout
parser.add_argument('-s','--stdout',
	default=0,
	action='store_true',
	help='output to stdout, no file be saved')

#parse arguments
args = parser.parse_args()

#read xyz into data frame
#skip first two rows of the xyz file
#only XMol xyz is supportet, atom(as element) x y z, e.g. C 1.58890 -1.44870 -0.47000
try:
	#get the two header lines and store it for later
	head = open(args.filename).readlines()[0:2]
	if head[-1] == "\n":
		head[-1] = " \n"
	head= ''.join(head).rstrip('\n')
	
	#xyz --> data frame
	xyz_df = pd.read_csv(args.filename, 
			delim_whitespace=True, 
			skiprows=2, 
			names=["element", "x", "y", "z"])
	#index +1, first atom is atom number 1
	xyz_df.index +=1
#file not found
except IOError:
	print(f"'{args.filename}'" + " not found")
	sys.exit(1)

#center molecule
#center on one or more atoms or all atoms
if args.origin:
	if args.origin == 12345:
		#if no atoms are given, 12345 is for no atoms
		#generate centroid of all atoms 
		origin=np.mean(xyz_df[['x','y','z']],axis=0)
	else:
		#generate centroid from selected atoms
		try:
			origin_atoms_df = xyz_df.loc[args.origin]
			origin=np.mean(origin_atoms_df[['x','y','z']],axis=0)
		except KeyError: 
			#exit in case of atoms are not in the xyz-file
			print('Warning! Atom not in xyz file. Exit.')
			sys.exit(1)

#copy xyz_df (data frame) to aligned_xyz_df 
aligned_xyz_df = xyz_df
#subtract centroid from xyz coordinates, 'zero' the molecule
aligned_xyz_df[['x','y','z']]=aligned_xyz_df[['x','y','z']].apply(lambda x: x-origin, axis=1)

#print detailed information
if args.verbose:
	print('')
	print('selected atom(s) for origin:')
	try:
		print(origin_atoms_df)
	except NameError:
		print('all atoms')
	print('')
	print('origin: %.4f %.4f %.4f' % (origin['x'],origin['y'],origin['z']))
	print('')

#align atom or centroid of atoms to the x-axis (1 0 0)
#e.g. -x 1 or -x 1 2 3
if args.x:
	try:
		#select input atoms from data frame 
		#move selected atoms in a new data frame
		x_atoms_df = aligned_xyz_df.loc[args.x]
	except KeyError: 
		#exit in case of atoms are not in the xyz-file
		print('Warning! Atom not in xyz file. Exit.')
		sys.exit(1)
		
	#generate centroid from selected atoms (or use coordinates of one atom)
	atoms_x=np.mean(x_atoms_df[['x','y','z']],axis=0)
	#align the atom or centroid vetor to the x-axis (1 0 0)
	aligned_xyz_df[['x','y','z']]=align_xyz((atoms_x['x'],atoms_x['y'],atoms_x['z']), (1,0,0), aligned_xyz_df[['x','y','z']])
	
	#print detailed information
	if args.verbose:
		with np.printoptions(precision=4, suppress=True, linewidth=100): 
			print('selected atom(s) for x direction:')
			print(x_atoms_df)
			print('')
			print('rotation matrix x: ', *rotmat_from_vec((atoms_x['x'],atoms_x['y'],atoms_x['z']), (1,0,0)).T)
			print('')
			
#align atom or centroid of atoms to the y-axis (0 1 0)
#e.g. -y 1 or -y 1 2 3
if args.y:
	try:
		#select input atoms from data frame 
		#move selected atoms in a new data frame
		y_atoms_df = aligned_xyz_df.loc[args.y]
	except KeyError: 
		#exit in case of atoms are not in the xyz-file
		print('Warning! Atom not in xyz file. Exit.')
		sys.exit(1)
	
	#generate centroid from selected atoms (or use coordinates of one atom)
	atoms_y=np.mean(y_atoms_df[['x','y','z']],axis=0)
	#align the atom or centroid vetor to the y-axis (0 1 0)
	aligned_xyz_df[['x','y','z']]=align_xyz((atoms_y['x'],atoms_y['y'],atoms_y['z']), (0,1,0), aligned_xyz_df[['x','y','z']])
	
	#print detailed information
	if args.verbose:
		with np.printoptions(precision=4, suppress=True, linewidth=100): 
			print('selected atom(s) for y direction:')
			print(y_atoms_df)
			print('')
			print('rotation matrix y: ', *rotmat_from_vec((atoms_y['x'],atoms_y['y'],atoms_y['z']), (0,1,0)).T)
			print('')

#align atom or centroid of atoms to the z-axis (0 0 1)
#e.g. -z 1 or -z 1 2 3
if args.z:
	try:
		#select input atoms from data frame 
		#move selected atoms in a new data frame
		z_atoms_df = aligned_xyz_df.loc[args.z]
	except KeyError: 
		#exit in case of atoms are not in the xyz-file
		print('Warning! Atom not in xyz file. Exit.')
		sys.exit(1)
		
	#generate centroid from selected atoms (or use coordinates of one atom)
	atoms_z=np.mean(z_atoms_df[['x','y','z']],axis=0)
	#align the atom or centroid vetor to the z-axis (0 0 1)
	aligned_xyz_df[['x','y','z']]=align_xyz((atoms_z['x'],atoms_z['y'],atoms_z['z']), (0,0,1), aligned_xyz_df[['x','y','z']])
	
	#print detailed information
	if args.verbose:
		with np.printoptions(precision=4, suppress=True, linewidth=100): 
			print('selected atom(s) for z direction:')
			print(z_atoms_df)
			print('')
			print('rotation matrix z: ', *rotmat_from_vec((atoms_z['x'],atoms_z['y'],atoms_z['z']), (0,0,1)).T)
			print('')

#align atom or centroid of atoms to the xyz-axes (1 1 1)
#take atoms from -x, -y, -z input
if args.x and args.y and args.z:
	
	#select input atoms from data frame 
	x_atoms_df = aligned_xyz_df.loc[args.x]
	y_atoms_df = aligned_xyz_df.loc[args.y]
	z_atoms_df = aligned_xyz_df.loc[args.z]
	
	#move selected atoms in a new data frame
	xyz_atoms_df = pd.concat([x_atoms_df, y_atoms_df,z_atoms_df])
	#generate centroid from selected atoms 
	atoms_xyz=np.mean(xyz_atoms_df[['x','y','z']],axis=0)
	#align the centroid vetor to the xyz-axes (1 1 1)
	aligned_xyz_df[['x','y','z']]=align_xyz((atoms_xyz['x'],atoms_xyz['y'],atoms_xyz['z']), (1,1,1), aligned_xyz_df[['x','y','z']])
	
	#print detailed information
	if args.verbose:
		with np.printoptions(precision=4, suppress=True, linewidth=100): 
			print('rotation matrix xyz: ', *rotmat_from_vec((atoms_xyz['x'],atoms_xyz['y'],atoms_xyz['z']), (1,1,1)).T)
			print('')
			
#align atom or centroid of atoms to the xy-axes (1 1 0)
#take atoms from -x, -y input
if args.x and args.y:
	
	#select input atoms from data frame 
	x_atoms_df = aligned_xyz_df.loc[args.x]
	y_atoms_df = aligned_xyz_df.loc[args.y]
	
	#move selected atoms in a new data frame
	xy_atoms_df = pd.concat([x_atoms_df, y_atoms_df])
	#generate centroid from selected atoms 
	atoms_xy=np.mean(xy_atoms_df[['x','y','z']],axis=0)
	#align the centroid vetor to the xy-axes (1 1 0)
	aligned_xyz_df[['x','y','z']]=align_xyz((atoms_xy['x'],atoms_xy['y'],atoms_xy['z']), (1,1,0), aligned_xyz_df[['x','y','z']])
	
	#print detailed information
	if args.verbose:
		with np.printoptions(precision=4, suppress=True, linewidth=100): 
			print('rotation matrix xy: ', *rotmat_from_vec((atoms_xy['x'],atoms_xy['y'],atoms_xy['z']), (1,1,0)).T)
			print('')

#re-align atom or centroid of atoms to the y-axis (0 1 0)
#take atoms from -y input
if args.y:
	y_atoms_df = aligned_xyz_df.loc[args.y]
	atoms_y=np.mean(y_atoms_df[['x','y','z']],axis=0)
	aligned_xyz_df[['x','y','z']]=align_xyz((atoms_y['x'],atoms_y['y'],atoms_y['z']), (0,1,0), aligned_xyz_df[['x','y','z']])
	
	#print detailed information
	if args.verbose:
		with np.printoptions(precision=4, suppress=True, linewidth=100): 
			print('rotation matrix y 2nd: ', *rotmat_from_vec((atoms_y['x'],atoms_y['y'],atoms_y['z']), (0,1,0)))
			print('')

#re-align atom or centroid of atoms to the x-axis (1 0 0)
#take atoms from -x input
if args.x:
		x_atoms_df = aligned_xyz_df.loc[args.x]
		atoms_x=np.mean(x_atoms_df[['x','y','z']],axis=0)
		aligned_xyz_df[['x','y','z']]=align_xyz((atoms_x['x'],atoms_x['y'],atoms_x['z']), (1,0,0), aligned_xyz_df[['x','y','z']])
		
		#print detailed information
		if args.verbose:
			with np.printoptions(precision=4, suppress=True, linewidth=100): 
				print('rotation matrix x 2nd: ', *rotmat_from_vec((atoms_x['x'],atoms_x['y'],atoms_x['z']), (1,0,0)).T)
				print('')

#rotate counterclockwise about the x-, y- and z-axes
#input is x y z, input is in degrees
#e.g. -r 0 90 0 - rotate 90° about the y-axis
#e.g. -r 45.11 90 0 - rotate 45° about the x-axis and 90° about the y-axis
if args.rotate:
	rotmatrix=rotmat_from_ang(args.rotate)
	aligned_xyz_df[['x','y','z']]=np.dot(aligned_xyz_df[['x','y','z']],rotmatrix.T)
	
	#print detailed information
	if args.verbose:
		with np.printoptions(precision=4, suppress=True, linewidth=100): 
			print('rotation matrix angles:\n ', *rotmatrix.T)
			
#process the rotation/transformation matrix on coordinates
#input is x1 y1 z1, x2 y2 z2, x3 y3 z3
#e.g. -m -1 0 0 0 -1 0 0 0 -1 - invert coordinates
if args.matrix:
	rot_matrix=np.asarray(args.matrix).reshape([3,3])
	aligned_xyz_df[['x','y','z']]=np.dot(aligned_xyz_df[['x','y','z']],rot_matrix)

#translate coordinates
#input is x y z, input is in Å or in units of the input coordinates
if args.translate:
	aligned_xyz_df[['x','y','z']]=np.add(aligned_xyz_df[['x','y','z']],args.translate)

#print data frame after all transformations / rotations
if args.verbose:
	print('transformed coordinates:')
	print(aligned_xyz_df)
	print('')

#print xyz-file to stdout, do not save the file
if args.stdout:
	pd.set_option('display.float_format', lambda x: '%12.8f' % x)
	print(head)
	print(aligned_xyz_df.to_string(index=False, header=False))

#save the modified xyz-file
#if not printed to stdout
if not args.stdout:
	try:
		my_numpy = aligned_xyz_df.to_numpy()
		np.savetxt(os.path.splitext(args.filename)[0] +'-mod.xyz', my_numpy,fmt='%-2s  %12.8f  %12.8f  %12.8f', delimiter='', header=head, comments='')
	#write errot -> exit here
	except IOError:
		print("Write error. Exit.")
		sys.exit(1)
