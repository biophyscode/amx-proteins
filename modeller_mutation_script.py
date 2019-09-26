#!/usr/bin/env python

"""
MODELLER homology modeling script
This script is called by automacs to prepare a protein structure.
"""

from modeller import *
from modeller.automodel import *
import os,sys,json

settings = {}
# optional settings
settings_fn = 'modeller_settings.json'
if os.path.isfile(settings_fn):
	with open(settings_fn) as fp:
		settings.update(**json.load(fp))

# generic names
# historical note: we used a json file to pass more information
target_name = settings.get('target','target')
template_name = settings.get('template','template')
refinement = settings.get('refinement','fast')
n_models = settings.get('n_models',2)
segment_ids = settings.get('segment_ids','A')
# ensure the resulting structures start from the correct residue number
renumber_residues = settings.get('renumber_residues',[1])

class mymodel(automodel):
	def special_patches(self, aln):
		self.rename_segments(segment_ids=segment_ids,
			renumber_residues=renumber_residues)
		
env = environ()
env.io.hetatm = False
env.io.water = False
env.io.atom_files_directory = ['./']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

a = mymodel(env,
	alnfile='%s.ali'%target_name,
	knowns=template_name,
	assess_methods=(assess.DOPE),
	sequence='%s'%target_name)

if refinement=='fast':
	a.library_schedule = autosched.fast
	a.max_var_iterations = 300
	a.md_level = refine.fast
if refinement=='slow':
	a.library_schedule = autosched.slow
	a.max_var_iterations = 300
	a.md_level = refine.slow

a.starting_model = 1
a.ending_model = n_models

a.make()
