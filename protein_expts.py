{

'protein':{
#####
####
###
##
#
'tags':['aamd','protein'],
'script':'protein.py',
'params':'parameters.py',
'extensions':[],
'settings':"""

step: protein                       # name of the folder is s01-protein
force field: charmm27               # which gromacs-standard force-field to use (see pdb2gmx list)
water: tip3p                        # which water model (another question from pdb2gmx)
equilibration: nvt-short,nvt,npt    # which equilibration step to use (must have `input-name-in.mdp` below)
pdb source: 1yrf                    # PDB code for download. overrides the start structure
start structure: None               # path to PDB structure or None to use a single PDB in inputs
protein water gap: 3.0              # Angstroms distance around the protein to remove water
water buffer: 1.2                   # distance (nm) of solvent to the box wall 
solvent: spc216                     # starting solvent box (use spc216 from gromacs share)
ionic strength: 0.150               # desired molar ionic strength
cation: NA                          # name of the cation for neutralizing the system
anion: CL                           # name of the anion for neutralizing the system

#---INTEGRATOR PARAMETERS generated via parameters.py
mdp_specs:| {
	'group':'aamd',
	'mdps':{
		'input-em-steep-in.mdp':['minimize'],
		'input-em-cg-in.mdp':['minimize',{'integrator':'cg'}],
		'input-md-nvt-eq-in.mdp':['nvt-protein','nvt-protein',{'nsteps':10000}],
		'input-md-nvt-short-eq-in.mdp':['nvt-protein-short',{'nsteps':10000}],
		'input-md-npt-eq-in.mdp':['npt-protein',{'nsteps':10000}],
		'input-md-in.mdp':{'nsteps':100000},
		},
	}
"""},

'generate_charmm_landscape':{
#####
####
###
##
#
#
'tags':['aamd','landscape'],
'params':'@bilayers/parameters.py',
'extensions':[
	'@charmm/landscape.py'],
'quick':"""

from amx import *
init()
write_charmm_landscape()

""",
'settings':"""
landscape at: @charmm/landscape.json
lipids: ['DOPC','DOPS','PI2P','CHL1','DOPE','POPC']
""",},

'bilayer_protein_topology_only':{
#####
####
###
##
#
'tags':['aamd','protein'],
'script':'protein-topology-only.py',
'params':None,
'extensions':[],
'settings':"""

step: protein                       # name of the folder is s01-protein
force field: charmm27               # which gromacs-standard force-field to use (see pdb2gmx list)
water: tip3p                        # which water model (another question from pdb2gmx)

"""},

}