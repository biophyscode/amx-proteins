#!/usr/bin/env python

"""
PROTEIN SIMULATION
Topology only!
"""

from amx import *

init()
make_step(settings.step)
if state.pdb_source: get_pdb(state.pdb_source)
else: get_start_structure(state.start_structure)
remove_hetero_atoms(
	structure='start-structure.pdb',
	out='start-structure-trim.pdb')
gmx('pdb2gmx',
	base='vacuum',
	structure='start-structure-trim.pdb',
	gro='vacuum-alone',
	water=settings.water,
	ff=settings.force_field,
	log='pdb2gmx')
extract_itp('system.top')
#---extract_itp always produces protein.itp
state.protein_prepared = {'gro':state.here+'vacuum-alone.gro','itp':state.here+'protein.itp'}
