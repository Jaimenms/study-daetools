__doc__ = """
This file concentrates all the test cases for the network simulation. So consider te functions beginning with the
name case_ and returning the data dictionary.

"""

def case_01():

    network_data = {
        'name': 'biofilm_test_01',
        'kind': 'biofilm',
        'module': 'models.biofilm',
        'class': 'Biofilm',
        'submodels': {
        },
        'specifications': {
            'Tbf': 273.15 + 37,
            'vbf': 0.54,
        },
        'parameters': {
            'rhomf': 980,
            'mfi': 1e-6,
        },
        'initial_guess': {
            'mf': 1e-6,
        }
    }


    def define_variables(self):

        # Variable types
        mass_biofilm_t = daeVariableType("mass_biofilm_t", (kg ** (1)) * (m ** (-2)), 1e-10, 1e+04, 1e-06, 1e-07)
        self.mf = daeVariable("mf", mass_biofilm_t, self, "Mass of Biofilm", [self.x, ])
        self.Tbf = daeVariable("Tbf", temperature_t, self, "Temperature of Biofilm", [self.x, ])


    def define_parameters(self):

        self.lagt = Constant(3600*24 * s)

        self.rhomf = daeParameter("rhomf",(kg ** (1)) * (m ** (-3)), self, "Density of the biofilm")
        self.mfi = daeParameter("mfi", (kg ** (1)) * (m ** (-2)), self, "Initial biofilm density")



    return network_data