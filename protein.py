#!/usr/bin/env python

"""
PROTEIN SIMULATION
Atomistic protein in water.
"""

from amx import *

make_step()
write_mdp()
interpret_start_structure()
remove_hetero_atoms(
	structure='start-structure.pdb',
	out='start-structure-trim.pdb')
gmx('pdb2gmx',
	posre='vacuum-posre.itp',
	structure='start-structure-trim.pdb',
	out='vacuum-alone.gro',
	water=settings.water,
	ff=settings.force_field,
	top='system.top',
	log='pdb2gmx')
copy_file('system.top','vacuum.top')
extract_itp('vacuum.top')
write_top('vacuum.top')
gmx('editconf',
	structure='vacuum-alone.gro',
	out='vacuum.gro',
	center=True,d='%.2f'%settings.water_buffer,
	log='editconf-vacuum-room')
minimize('vacuum',method='steep')
solvate_protein(
	structure='vacuum-minimized.gro',
	top='vacuum.top')
minimize('solvate')
counterions(
	structure='solvate-minimized.gro',
	top='solvate',
	ff_includes='ions')
minimize('counterions')
write_structure_pdb(
	pdb='start-structure.pdb',
	structure='counterions.gro')
write_top('system.top')
write_continue_script()
equilibrate()
