from SloppyCell.ReactionNetworks import *
import numpy as np

class wiring:
	"""Class that stores info about a wiring and its performance"""
	def __init__(self, net, id):
		self.basic_net = net
		self.id = id
		self.nets = dict()
		self.expts = dict()
		self.opt = dict()
						    
	def add_net(self, net):
		self.nets[net.id] = net

	def add_expt(self, expt):
		self.expts[expt.name] = expt
									    
	def add_opt(self, conditions, id):
		nets = []
		expts = []
		params = []
		for c in conditions:
			if c not in self.nets:
				raise NameError('no net for given condition')
			nets.append(self.nets[c])
			params+=self.nets[c].optimizableVars.keys()
			if c not in self.expts:
				raise NameError('no expt for given condition')
			expts.append(self.expts[c])
		m = Model(expts,nets)
		for p in list(set(params)):    
			res = Residuals.PriorInLog(p+'_prior', p, np.log(1), np.log(np.sqrt(100)))
			m.AddResidual(res)
		self.opt[id] = optimization(m,id)

	def add_condition(self, mutant, intervention, data, sf=None):
		id = mutant.id+'_'+intervention.id
		net = self.basic_net.copy(new_id=id)
		for param, val in mutant.settings.items():
			net.set_var_ic(param, val)
			net.set_var_optimizable(param, is_optimizable=False)
		if len(intervention.settings)>0:
			for var, val in intervention.settings.items():
				net.set_var_constant(var, is_constant=False)
			net.add_event(id=id+'_event', trigger=intervention.trigger, event_assignments=intervention.settings)
		self.add_net(net)
		expt = Experiment(id)
		expt.set_data({id: data})
		if sf is not None:
			expt.set_fixed_sf({key: sf for key in data.keys()})
		self.add_expt(expt)


class mutant:
	"""Description goes here"""
	def __init__(self, id):
		self.settings = dict()
		self.id = id

	def add_setting(self, param, val):
		self.settings[param] = val


class intervention:
	"""Description goes here"""
	def __init__(self, id):
		self.settings = dict()
		self.trigger = ''
		self.id = id

	def add_setting(self, var, val):
		self.settings[var] = val
	
	def add_trigger(self, trigger):
		self.trigger = trigger


class optimization:
	"""Description goes here"""
	def __init__(self, m, id):
		self.id = id
		self.model = m.copy()
		self.params_init = m.params.copy()
		self.cost_init = m.cost(m.params)
		self.params_opt = []
		self.cost_opt = np.inf

	def set_params_init(self, params):
		self.params_init = params.copy()
		self.cost_init = self.model.cost(params)
	
	def run(self, maxiter=200, disp=False):
		p = Optimization.fmin_log_params(self.model,
			self.params_init, disp=disp)
		self.cost_opt = self.model.cost(p)
		self.params_opt = p

mutants = dict()
mutants['wt'] = mutant('wt')

mutants['dkard'] = mutant('dkard')
mutants['dkard'].add_setting('pmelt_bub1_pp2a', 0)
mutants['dkard'].add_setting('kb_pmelt_bub1_pp2a', 0)
mutants['dkard'].add_setting('kd_pmelt_bub1_pp2a', 0)

interventions = dict()
interventions['ic'] = intervention('ic')

interventions['mps1'] = intervention('mps1')
interventions['mps1'].add_setting('mps1', 0)
interventions['mps1'].add_trigger('gt(time, 1000)')

interventions['aurb'] = intervention('aurb')
interventions['aurb'].add_setting('paurb', 0)
interventions['aurb'].add_setting('aurb_tot', 0)
interventions['aurb'].add_trigger('gt(time, 1000)')

interventions['both'] = intervention('both')
interventions['both'].add_setting('paurb', 0)
interventions['both'].add_setting('aurb_tot', 0)
interventions['both'].add_trigger('gt(time, 1000)')
interventions['both'].add_setting('mps1', 0)
interventions['both'].add_trigger('gt(time, 1000)')

sf = {
	'wt': {
		'ic': 1,
		'mps1': None,
		'aurb': None,
		'both': None
		},
	'dkard': {
		'ic': 1,
		'mps1': None,
		'aurb': None,
		'both': None
		}
	}

data = {
	'wt': {
		'ic':{
			'prvsf': {
					100: (0.5, 0.2),
					1000: (0.5, 0.2)},              
			'pmelt_tot' : {
					100: (0.5, 0.2),
					1100: (0.5, 0.2)}
			},
		'mps1':{
			'prvsf': {
					100: (1, 0.2),
					1000: (1, 0.2)},              
			'pmelt_tot' : {
					100: (1, 0.2),
					1100: (0, 0.2)}, 
			},
		'aurb':{
			'pmelt_tot' : {
					100: (1, 0.2),
					1100: (0, 0.2)}, 
			},
		'both':{
			'pmelt_tot' : {
					100: (1, 0.2),
					1100: (0, 0.2)}
			}
		},
	'dkard': {
		  'ic':{
			'prvsf': {
					100: (0.9, 0.2),
					1000: (0.9, 0.2)},              
			'pmelt_tot' : {
					100: (0.9, 0.2),
					1100: (0.9, 0.2)}
			  },
		  'mps1':{
    			'pmelt_tot' : {
					100: (1, 0.2),
					1100: (1, 0.2)
				}
			  },
		  'aurb':{
    			'prvsf' : {
					100: (1, 0.2),
	                		1100: (0, 0.2)},                              
			  },
		  'both':{
			 'prvsf' : {
				 	100: (1, 0.2),
					1100: (0, 0.2)}
			  }
		  }

	}


