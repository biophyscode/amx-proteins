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
itp = extract_itp('system.top')
#---extract_itp always produces protein.itp
last_call = gmx_get_last_call('pdb2gmx')
prepped_files = {'itp':state.here+itp}
for flag,key in [('-i','posre'),('-o','gro')]: 
	prepped_files[key] = state.here+last_call['flags'][flag]
state.protein_prepared = prepped_files
