import json

class ParametersHandler:
    def __init__(self,filename):
        with open(filename) as f:
            self.params = json.load(f)

        self.solver_type = self.params.get("solver_type")
        self.time_integration = self.params.get("time_integration")
        self.tolerance = self.params.get("nonlin_tolerance")
        if self.tolerance == None:
            self.tolerance = 1e-6
        self.base = self.params.get("delta_time_base")
        self.esponente = self.params.get("delta_time_negative_esponent")

        self.kappa = self.params.get("kappa")
        self.rho_fluid = self.params.get("rho_fluid")
        self.rho_structure = self.params.get("rho_structure")
        self.nu = self.params.get("nu")

        self.n_delta_x = self.params.get("fluid_triangulation_intervals")
        self.n_delta_s = self.params.get("structure_triangulation_intervals")

        self.dt = self.base*10**(-self.esponente)
        self.n_times = self.params.get("time_intervals_to_be_simulated")
        self.no_print_intervals = self.params.get("time_intervals_in_between_printed_results")
        self.stampa = []
        for i in range(self.n_times):
            if (i%self.no_print_intervals==0):
                self.stampa.append(True)
            else:
                self.stampa.append(False)
        self.mesh_prefix = self.params.get("mesh_prefix")

        self.mesh_name = self.mesh_prefix+str(int(self.n_delta_x))

        if self.solver_type == "ns_":
            sp =  self.solver_type+self.mesh_name+'_'
            sp += self.time_integration+'_'
            sp += 'dt'+str(int(self.base))+'em'+str(int(self.esponente))
            sp += '_hx'+str(int(self.n_delta_x))
            sp += '_nu'+str(self.nu)
        elif self.solver_type == "dlm_":
            sp = self.solver_type+self.mesh_name+'_'
            sp += self.time_integration+'_'
            sp += 'dt'+str(int(self.base))+'em'+str(int(self.esponente))
            sp += '_hx'+str(int(self.n_delta_x))+'_hs'+str(int(self.n_delta_s))
            sp += '_k'+str(self.kappa)
            sp += '_nu'+str(self.nu)

        self.sim_prefix = sp

        self.time_index_digits = self.params.get("time_index_digits")
        self.results_directory = self.params.get("results_directory")

    def simulation_info(self):
        print('-----------------------------------')
        print('started simulation: ')
        print(self.sim_prefix)
        print('dt = '+ str(self.dt))
        print('hx = '+ str(int(self.n_delta_x)))
        if self.solver_type == 'dlm_':
            print('hs = '+str(int(self.n_delta_s)))
            print('k  = '+str(self.kappa))
        print('nu = '+str(self.nu))
        print('tol = ' + str(self.tolerance))
        print('-----------------------------------')
        return

    def dump_to_json(self,filename):
        with open(filename, 'w') as outfile:
            json.dump(self.params, outfile)
