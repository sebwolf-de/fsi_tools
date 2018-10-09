import json

class ParametersHandler:
    def __init__(self,filename):
        with open(filename) as f:
            self.params = json.load(f)

        self.solver_type = self.params["solver_type"]
        self.base = self.params["delta_time_base"]
        self.esponente = self.params["delta_time_negative_esponent"]
        self.kappa = self.params["structure_stiffness_kappa"]
        self.reynolds = self.params["reynolds_number"]

        self.n_delta_x = self.params["fluid_triangulation_intervals"]
        self.n_delta_s = self.params["structure_triangulation_intervals"]

        self.dt = self.base*10**(-self.esponente)
        self.n_times = self.params["time_intervals_to_be_simulated"]
        self.no_print_intervals = self.params["time_intervals_in_between_printed_results"]
        self.stampa = []
        for i in range(self.n_times):
            if (i%self.no_print_intervals==0):
                self.stampa.append(True)
            else:
                self.stampa.append(False)
        self.mesh_prefix = self.params["mesh_prefix"]

        self.mesh_name = self.mesh_prefix+str(int(self.n_delta_s))

        bool_conversion = {"true_string" : True, "false_string" : False}
        self.equilibrium_at_zero = bool_conversion[self.params["equilibrium_at_zero"]]

        if self.solver_type == "ns_":
            sp =  self.solver_type+self.mesh_name+'_'
            sp += 'dt'+str(int(self.base))+'em'+str(int(self.esponente))
            sp += '_hx'+str(int(self.n_delta_x))
            sp += '_re'+str(int(self.reynolds))
        else:
            sp = self.solver_type+self.mesh_name+'_'
            sp += 'dt'+str(int(self.base))+'em'+str(int(self.esponente))
            sp += '_hx'+str(int(self.n_delta_x))+'_hs'+str(int(self.n_delta_s))
            sp += '_k'+str(int(self.kappa))
            sp += '_re'+str(int(self.reynolds))
            sp += '_eq_at_zero_'+str(self.equilibrium_at_zero)

        self.sim_prefix = sp

        self.time_index_digits = self.params["time_index_digits"]
        self.results_directory = self.params["results_directory"]

    def simulation_info(self):
        print '-----------------------------------'
        print 'started simulation: '
        print self.sim_prefix
        print 'dt = '+str(int(self.base))+'em'+str(int(self.esponente))
        print 'hx = '+str(int(self.n_delta_x))
        print 'hs = '+str(int(self.n_delta_s))
        print 'k  = '+str(int(self.kappa))
        print 're = '+str(int(self.reynolds))
        print 'eq_at_zero = '+str(self.equilibrium_at_zero)
        print '-----------------------------------'
        return

    def dump_to_json(self,filename):
        with open(filename, 'w') as outfile:
            json.dump(self.params, outfile)
