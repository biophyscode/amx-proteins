#!/usr/bin/env python

"""
PROTEIN MUTATION
Make a point mutation in a protein.
Development notes:
- pending multiple chain handling
"""

from amx import *
import re,shutil,glob,textwrap
from ortho import Handler,bash

### HELPERS

def get_start_structure(path):
	"""
	Get a start structure or auto-detect a single PDB in the inputs folder.
	"""
	if path: 
		altpath = os.path.join(globals().get(
			'expt',{}).get('cwd_source','./'),path)
	else: altpath = None
	if path and os.path.isfile(path): fn = path
	elif altpath and os.path.isfile(altpath): fn = altpath
	else: 
		fns = glob.glob('inputs/*.pdb')
		if len(fns)>1: raise Exception('multiple PDBs in inputs')
		elif len(fns)==0: raise Exception('no PDBs in inputs')
		else: fn = fns[0]
	shutil.copy(fn,os.path.join(state.here,''))
	shutil.copyfile(fn,os.path.join(state.here,'start-structure.pdb'))

### SETTINGS

target_name = 'target'
template_pdb = 'start-structure.pdb'

### MAIN

# track mutation state
if state.mutation: 
	raise Exception('preexisting mutation state: %s'%str(state.mutation))
state.mutation = {}

make_step(settings.step)

# check for PDB code or path
is_pdb = (settings.start_structure!=None 
	and re.match('^[A-Za-z0-9]{4}$',settings.start_structure.strip())!=None)
is_path = (settings.start_structure==None 
	or os.path.isfile(settings.start_structure))
if is_path + is_pdb != 1:
	raise Exception(
		'is_pdb = %s, is_path = %s, settings.start_structure = %s'%(
		is_path,is_pdb,settings.start_structure))
# collect the structure
if is_path or settings.start_structure==None:
	#! fix the above
	get_start_structure(settings.start_structure)
elif is_pdb:
	# clean the PDB code here
	settings.start_structure = settings.start_structure.strip().upper()
	get_pdb(settings.start_structure)

# get the sequence
import Bio
import Bio.PDB
from Bio import SeqIO
# hide warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
# via: https://stackoverflow.com/questions/9134795
import warnings
warnings.simplefilter('ignore',PDBConstructionWarning)
target = state.here+'start-structure.pdb'
records = list(SeqIO.parse(target,"pdb-atom"))
if len(records)!=1:
	raise Exception('development note: configured for a single chain only')
record = records[0]
seq = str(record.seq)
# get the start position from the PDB
pdb_start_position = record.annotations.get('start',1)
pdb_chain = record.annotations.get('chain','')

#! dev: note the single-chain requirement must be relaxed

# apply sequence adjustments
if settings.sequence_check:
	class SeqShifter(Handler):
		def confirm(self,seq,subseq,start,chain='',target=None,actual=None):
			candidates = [seq]
			if target: candidates += [target]
			if actual: candidates += [actual]
			for candidate in candidates:
				found = re.findall(subseq,candidate)
				if len(found)!=1:
					raise Exception(
						'failed to find subsequence %s in sequence: %s'%(
							candidate,seq))
			# expected start for the PDB marked in the ali file written below
			global pdb_start_resid
			if (pdb_start_position != 1 
				and pdb_start_position != start):
				raise Exception(('start position from the PDB is %d but '
					'the sequence_check setting says it should be %d and we '
					'only allow 1 or the correct value to avoid losing '
					'information')%(pdb_start_position,start))
			elif pdb_start_position == 1:
				# tell modeller the start position is 1
				pdb_start_resid = 1
			elif pdb_start_position == start:
				pdb_start_resid = start
			# override the target sequence with a full sequence if provided
			if target: seq = target
			return dict(start_position=start,seq=seq,chain=chain)
	# get sequence data
	state.mutation.update(**SeqShifter(
		seq=seq,**settings.sequence_check).solve)
else: 
	pdb_start_resid = pdb_start_position
	state.mutation.update(
		start_position=settings.sequence_check.start,
		chain=pdb_chain,seq=seq)

#! dev: extract the full information from the PDB header if missing residues
# override the target and PDB sequence
# note that we also use the full sequence in the SeqShifter above
# the ability to override the full sequence allows extensions and gap filling
seq = settings.sequence_check.get('target',seq)
pdb_seq = settings.sequence_check.get('actual',seq)

# apply the mutations
seq_target = list(seq)
if not settings.mutations: 
	raise Exception('we require mutations in settings')
for chain,mutations in settings.mutations.items():
	#! only handles point mutations
	for mutation in mutations:
		regex_point_mutation = \
			r'^(?P<source>[A-Za-z])(?P<spot>\d+)(?P<target>[A-Za-z]+)$'
		mutate = re.match(regex_point_mutation,mutation).groupdict()
		mutate['spot'] = int(mutate['spot'])
		pos_rel = mutate['spot']-state.mutation['start_position']
		# check the residue at the current position
		found = seq_target[pos_rel]
		if found != mutate['source']:
			raise Exception(
				'found %s at position %d and not expected %s'%(
					found,mutate['spot'],mutate['source']))
		# note that all changes are cumulative
		print('status applying mutation %s'%mutation)
		seq_target[pos_rel] = mutate['target']
# set the sequences
seq_target = ''.join(seq_target)
seq_template = state.mutation['seq']

# prepare the ali file
wrap_ali = lambda x,w=75: '\n'.join(textwrap.wrap(x,width=w))
template_ali = (">P1;%(name)s\n%(kind)s:%(pdb)s:"
"%(start_resid)s:%(start_chain)s:%(end_resid)s:%(end_chain)s:"
"protein_name:source:resolution:rfactor\n%(seq)s*\n")
with open(state.here+'%s.ali'%target_name,'w') as fp: 
	#! assume one chain
	# write the template
	fp.write(template_ali%dict(seq=wrap_ali(pdb_seq,w=75),
		pdb=re.match(r'^(.+)\.pdb',template_pdb).group(1),
		name='template',kind='structure',
		start_chain=state.mutation['chain'],
		start_resid=pdb_start_resid,
		end_resid=state.mutation['start_position']+len(seq),
		end_chain=''))
	# write the target
	fp.write(template_ali%dict(seq=wrap_ali(seq_template,w=75),
		pdb='XXXX',name='target',kind='sequence',
		start_chain=state.mutation['chain'],
		start_resid=state.mutation['start_position'],
		end_resid=state.mutation['start_position']+len(seq_template),
		end_chain=''))

# export settings
settings_out = dict(
	n_models=2,
	segment_ids=state.mutation['chain'],
	# names from the ali file
	target_name='target',
	template_name='template',
	refinement=settings.get('refinement','fast'),
	renumber_residues=state.mutation['start_position'])
settings_fn = 'modeller_settings.json'
with open(settings_fn,'w') as fp:
	json.dump(settings_out,fp)

# run modeller
bash('python %s'%os.path.join(os.getcwd(),settings.modeller_script),
	cwd=state.here,log=state.here+'log-modeller')

# select the best structure
with open(state.here+'log-modeller') as fp: 
	log_text = fp.read()
regex_get_best = \
	r">> Summary of successfully produced models:\n(.*?)\n(?:-+)\n(.*?)(?:\Z)"
result = re.search(regex_get_best,log_text,flags=re.M+re.DOTALL).groups()
#! check this header
result_header = result[0] 
results = result[1].splitlines()
if not result:
	raise Exception('failed to get best result. see log: %s'%(
		state.here+'log-modeller'))
dopes = [i.split()[2] for i in result[1:]]
best_ind = [ii for ii,i in enumerate(dopes) if i==min(dopes)][0]
best_fn = results[best_ind].split()[0]
# best structure is symlinked
os.symlink(state.here+best_fn,state.here+'structure-out.pdb')
print('status homology model is complete and saved to structure-out.pdb')