import json

class ParametersHandler:
    def __init__(self,filename):
        with open(filename) as f:
            params = json.load(f)

        self.base = params["delta_time_base"]
        self.esponente = params["delta_time_negative_esponent"]
        self.kappa = params["structure_stiffness_kappa"]
        self.reynolds = params["reynolds_number"]

        self.n_delta_x = params["fluid_triangulation_intervals"]
        self.n_delta_s = params["structure_triangulation_intervals"]

        self.dt = self.base*10**(-self.esponente)
        self.n_times = params["time_intervals_to_be_simulated"]
        self.no_print_intervals = params["time_intervals_in_between_printed_results"]
        self.stampa = []
        for i in range(self.n_times):
            if (i%self.no_print_intervals==0):
                self.stampa.append(True)
            else:
                self.stampa.append(False)
        self.mesh_prefix = params["mesh_prefix"]

        self.mesh_name = self.mesh_prefix+str(int(self.n_delta_s))

        bool_conversion = {"true_string" : True, "false_string" : False}
        self.equilibrium_at_zero = bool_conversion[params["equilibrium_at_zero"]]

        sp = 'dlm_'+self.mesh_name+'_'
        sp += 'dt'+str(int(self.base))+'em'+str(int(self.esponente))
        sp += '_hx'+str(int(self.n_delta_x))+'_hs'+str(int(self.n_delta_s))
        sp += '_k'+str(int(self.kappa))
        sp += '_re'+str(int(self.reynolds))
        sp += '_eq_at_zero_'+str(self.equilibrium_at_zero)

        self.sim_prefix = sp

        self.time_index_digits = params["time_index_digits"]
        self.results_directory = params["results_directory"]

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
