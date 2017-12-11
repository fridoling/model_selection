## Description goes here

from SloppyCell.ReactionNetworks import *
import numpy as np
import matplotlib.pyplot as plt
import itertools as it

## load full network and set dynamic vars
net_full = IO.from_SBML_file('../model/model.xml', id='net_full')
dyn_vars =  ['paurb', 'pmelt_bub1','pmelt_bub1_pp2a', 'pmelt', 'rvsf', 'rvsf_pp1']
nonvars = list(set(net_full.dynamicVars.keys()) - set(dyn_vars))
for key in nonvars:
	net_full.set_var_constant(key, is_constant=True)
	net_full.set_var_optimizable(key, is_optimizable=False)
for var in ['pp1_tot', 'bub1_tot', 'pp2a_tot', 'kln1_tot', 'aurb_tot']:
	net_full.set_var_optimizable(var, False)

## generate all possible models
## ids are strings of 1/0 indicating presence/absence of reactions
phospho_pars = [
	['kdp_rvsf_pp1', 'kdp_rvsf_pp2a'],
	['kp_aurb_bub1','kp_aurb_aurb'],
	['kdp_aurb_pp1','kdp_aurb_pp2a']]
nets = {}
alf = ['01','10','11']
net_iterator = it.product(alf,repeat=3)
for n_i in net_iterator:
	id = ''.join(n_i)
	net_it = net_full.copy(new_id='net_'+id)
	for i in range(len(n_i)):
		pars = phospho_pars[i]
		for j in range(len(n_i[i])):
			val = float(n_i[i][j])
			net_it.set_var_ic(pars[j], val)
			net_it.set_var_optimizable(pars[j], bool(val))
	nets[id] = net_it


