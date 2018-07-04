__doc__ = """
This file concentrates all the test cases for the network simulation. So consider te functions beginning with the
name case_ and returning the data dictionary.

"""

from water_properties import density, heat_capacity, conductivity, viscosity


def case_pipe():

    p1 = 400000 - 1669
    p2 = 400000 - 1669 - 65 - 1845.8 - 65 - 1062
    t1 = 305.15
    t2 = 328.35

    p1t = p1 - 65
    p2t = p2 + 65

    q = 0.10425

    Pm = 0.5 * (p1t + p2t)
    Tm = 0.5 * (t1 + t2)

    network_data = {
        'name': 'network_test_00',
        'kind': 'network',
        'submodels': {
            'node_A': {
                'kind': 'node',
                'module': 'models.source',
                'class': 'Source',
                'specifications': {
                    'P': p1,
                },
                'parameters': {
                    'Text': t1,
                    'Pext': p1,
                    'x': 0.0,
                    'y': 0.0,
                    'z': 0.0
                },
                'initial_guess':{
                    'w': q,
                    'T': t1,
                    'P': p1,
                }
            },
            'node_B': {
                'kind': 'node',
                'module': 'models.sink',
                'class': 'Sink',
                'specifications': {
                    'P': p2,
                },
                'parameters': {
                    'Text': t2,
                    'Pext': p2,
                    'x': 0.0,
                    'y': 1.0,
                    'z': 0.0
                },
                'initial_guess':{
                    'w': q,
                    'T': t2,
                    'P': p2,
                }
            },
            'pipe_01': {
                'kind': 'edge',
                'module': 'models.pipe',
                'class': 'Pipe',
                'from': 'node_A',
                'to': 'node_B',
                "domains": {
                    "x": {
                        "initial": 0.0,
                        "final": 1.0,
                        "N": 10
                    },
                },
                'states': {
                },
                'parameters': {
                    'Di': 0.01575,
                    'tetha': 0.0,
                    'L': 6.1,
                    'ep': 45e-6,
                    'Klb': 0.5*0.9,
                    'Kub': 0.5*0.9,
                    'Npipes': 1.0,
                },
                'specifications':{

                },
                'initial_guess': {
                    'Re': 13793,
                    'D': 0.01575,
                    'v': 0.54,
                    'k': q,
                    'T': 300.0,
                    'P':  {
                      "initial": p1t,
                      "final": p2t
                    },
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True)
                },
            },
        },
    }


    return network_data



def case_fixed_external_temperature_pipe():

    p1 = 400000 - 1669
    p2 = 400000 - 1669 - 65 - 1845.8 - 65 - 1062
    t1 = 305.15
    t2 = 328.35

    p1t = p1 - 65
    p2t = p2 + 65

    q = 0.10425

    Pm = 0.5 * (p1t + p2t)
    Tm = 0.5 * (t1 + t2)

    network_data = {
        'name': 'network_test_01',
        'kind': 'network',
        'submodels': {
            'node_A': {
                'kind': 'node',
                'module': 'models.source',
                'class': 'Source',
                'specifications': {
                    'P': p1,
                },
                'parameters': {
                    'Text': t1,
                    'Pext': p1,
                    'x': 0.0,
                    'y': 0.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': t1,
                    'P': p1,
                }
            },
            'node_B': {
                'kind': 'node',
                'module': 'models.sink',
                'class': 'Sink',
                'specifications': {
                    'P': p2,
                },
                'parameters': {
                    'Text': t2,
                    'Pext': p2,
                    'x': 0.0,
                    'y': 1.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': t2,
                    'P': p2,
                }
            },
            'pipe_01': {
                'kind': 'edge',
                'module': 'models.fixed_external_temperature_pipe',
                'class': 'FixedExternalTemperaturePipe',
                'from': 'node_A',
                'to': 'node_B',
                "domains": {
                    "x": {
                        "initial": 0.0,
                        "final": 1.0,
                        "N": 10
                    },
                },
                'states': {
                },
                'parameters': {
                    'Di': 0.01575,
                    'Do': 0.01905,
                    'kwall': 16.0,
                    'tetha': 0.0,
                    'L': 6.1,
                    'ep': 45e-6,
                    'Klb': 0.5 * 0.9,
                    'Kub': 0.5 * 0.9,
                    'Npipes': 1.0,
                    'Text': 328.,
                },
                'specifications': {

                },
                'initial_guess': {
                    'Re': 13793,
                    'D': 0.01575,
                    'v': 0.54,
                    'k': q,
                    'T':  {
                      "initial": t1,
                      "final": t2
                    },
                    'P':  {
                      "initial": p1t,
                      "final": p2t
                    },
                    'fD': 0.0330,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                    'Qout':{
                      "initial": -3526.5,
                      "final":-605.2
                    },
                    'To': 328.,
                    'Ti': {
                      "initial": 323.07,
                      "final": 331.42
                    },
                    'hint': 3976.69,
                    'Resistance': 0.007941,
                },
            },
        },
    }

    return network_data


def case_fixed_external_convection_pipe():

    p1 = 400000 - 1669
    p2 = 400000 - 1669 - 65 - 1845.8 - 65 - 1062
    t1 = 305.15
    t2 = 328.35

    p1t = p1 - 65
    p2t = p2 + 65

    q = 0.10425

    Pm = 0.5 * (p1t + p2t)
    Tm = 0.5 * (t1 + t2)

    network_data = {
        'name': 'network_test_01',
        'kind': 'network',
        'submodels': {
            'node_A': {
                'kind': 'node',
                'module': 'models.source',
                'class': 'Source',
                'specifications': {
                    'P': p1,
                },
                'parameters': {
                    'Text': t1,
                    'Pext': p1,
                    'x': 0.0,
                    'y': 0.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': t1,
                    'P': p1,
                }
            },
            'node_B': {
                'kind': 'node',
                'module': 'models.sink',
                'class': 'Sink',
                'specifications': {
                    'P': p2,
                },
                'parameters': {
                    'Text': t2,
                    'Pext': p2,
                    'x': 0.0,
                    'y': 1.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': t2,
                    'P': p2,
                }
            },
            'pipe_01': {
                'kind': 'edge',
                'module': 'models.fixed_external_convection_pipe',
                'class': 'FixedExternalConvectionPipe',
                'from': 'node_A',
                'to': 'node_B',
                "domains": {
                    "x": {
                        "initial": 0.0,
                        "final": 1.0,
                        "N": 10
                    },
                },
                'states': {
                },
                'parameters': {
                    'Di': 0.01575,
                    'Do': 0.01905,
                    'kwall': 16.0,
                    'tetha': 0.0,
                    'L': 6.1,
                    'ep': 45e-6,
                    'Klb': 0.5 * 0.9,
                    'Kub': 0.5 * 0.9,
                    'Npipes': 1.0,
                    'Text': 273.15 + 60.,
                    'hext': 17284.56,
                },
                'specifications': {

                },
                'initial_guess': {
                    'Re': 13793,
                    'D': 0.01575,
                    'v': 0.54,
                    'k': q,
                    'T':  {
                      "initial": t1,
                      "final": t2
                    },
                    'P':  {
                      "initial": p1t,
                      "final": p2t
                    },
                    'fD': 0.0330,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                    'Qout': {
                        "initial": -3526.5,
                        "final": -605.2
                    },
                    'To':  {
                      "initial": 329.75,
                      "final": 332.57
                    },
                    'Ti': {
                        "initial": 323.07,
                        "final": 331.42
                    },
                    'hint': 3976.69,
                    'Resistance': 0.007941,
                },
            },
        },
    }

    return network_data


def case_external_film_condensation_pipe():

    p1 = 400000 - 1669
    p2 = 400000 - 1669 - 65 - 1845.8 - 65 - 1062
    t1 = 305.15
    t2 = 328.35

    p1t = p1 - 65
    p2t = p2 + 65

    q = 0.10425

    Pm = 0.5 * (p1t + p2t)
    Tm = 0.5 * (t1 + t2)

    network_data = {
        'name': 'network_test_01',
        'class': 'Network',
        'submodels': {
            'node_A': {
                'kind': 'node',
                'module': 'models.source',
                'class': 'Source',
                'specifications': {
                    'P': p1,
                },
                'parameters': {
                    'Text': t1,
                    'Pext': p1,
                    'x': 0.0,
                    'y': 0.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': t1,
                    'P': p1,
                }
            },
            'node_B': {
                'kind': 'node',
                'module': 'models.sink',
                'class': 'Sink',
                'specifications': {
                    'P': p2,
                },
                'parameters': {
                    'Text': t2,
                    'Pext': p2,
                    'x': 0.0,
                    'y': 1.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': t2,
                    'P': p2,
                }
            },
            'pipe_01': {
                'kind': 'edge',
                'module': 'models.external_film_condensation_pipe',
                'class': 'ExternalFilmCondensationPipe',
                'from': 'node_A',
                'to': 'node_B',
                "domains": {
                    "x": {
                        "initial": 0.0,
                        "final": 1.0,
                        "N": 10
                    },
                },
                'states': {
                },
                'parameters': {
                    'Di': 0.01575,
                    'Do': 0.01905,
                    'kwall': 16.0,
                    'tetha': 0.0,
                    'L': 6.1,
                    'ep': 45e-6,
                    'Klb': 0.5 * 0.9,
                    'Kub': 0.5 * 0.9,
                    'Npipes': 1.0,
                    'rhov':  0.13044 ,
                    'hvap': 2357680,
                    'Tsat': 60 + 273.15,
                    'Psat': 19950,
                    'fNtub': 1.0,
                },
                'specifications': {

                },
                'initial_guess': {
                    'Re': 13793,
                    'D': 0.01575,
                    'v': 0.54,
                    'k': q,
                    'T':  {
                      "initial": t1,
                      "final": t2
                    },
                    'P':  {
                      "initial": p1t,
                      "final": p2t
                    },
                    'fD': 0.0330,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                    'Qout': {
                        "initial": -3526.5,
                        "final": -605.2
                    },
                    'To': {
                        "initial": 329.75,
                        "final": 332.57
                    },
                    'Ti': {
                        "initial": 323.07,
                        "final": 331.42
                    },
                    'hint': 3976.69,
                    'Resistance': 0.007941,
                },
            },
        },
    }

    return network_data


def case_biofilmed_pipe():

    p1 = 400000 - 1669
    p2 = 400000 - 1669 - 65 - 1845.8 - 65 - 1062
    t1 = 305.15
    t2 = 328.35

    p1t = p1 - 65
    p2t = p2 + 65

    q = 0.10425

    Pm = 0.5 * (p1t + p2t)
    Tm = 0.5 * (t1 + t2)

    network_data = {
        'name': 'network_test_01',
        'kind': 'network',
        'submodels': {
            'node_A': {
                'kind': 'node',
                'module': 'models.source',
                'class': 'Source',
                'specifications': {
                    'P': p1,
                },
                'parameters': {
                    'Text': t1,
                    'Pext': p1,
                    'x': 0.0,
                    'y': 0.0,
                    'z': 0.0
                },
                'initial_guess':{
                    'w': q,
                    'T': t1,
                    'P': p1,
                }
            },
            'node_B': {
                'kind': 'node',
                'module': 'models.sink',
                'class': 'Sink',
                'specifications': {
                    'P': p2,
                },
                'parameters': {
                    'Text': t2,
                    'Pext': p2,
                    'x': 0.0,
                    'y': 1.0,
                    'z': 0.0
                },
                'initial_guess':{
                    'w': q,
                    'T': t2,
                    'P': p2,
                }
            },
            'pipe_01': {
                'kind': 'edge',
                'module': 'models.biofilmed_pipe',
                'class': 'BiofilmedPipe',
                'from': 'node_A',
                'to': 'node_B',
                "domains": {
                    "x": {
                        "initial": 0.0,
                        "final": 1.0,
                        "N": 10
                    },
                },
                'states': {
                },
                'parameters': {
                    'Di': 0.01575,
                    'tetha': 0.0,
                    'L': 6.1,
                    'ep': 45e-6,
                    'Klb': 0.5*0.9,
                    'Kub': 0.5*0.9,
                    'Npipes': 1.0,
                    'rhomf': 980.,
                    'mfi': 1e-6,
                },
                'specifications': {
                },
                'initial_guess': {
                    'Re': 13793,
                    'D': 0.01575,
                    'v': 0.54,
                    'k': q,
                    'T':  {
                      "initial": t1,
                      "final": t2
                    },
                    'P':  {
                      "initial": p1t,
                      "final": p2t
                    },
                    'fD': 0.0330,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                    'Tbf': t1,
                    'vbf': 0.54,
                },
            },
        },
    }

    return network_data

def case_biofilmed_fixed_external_convection_pipe():
    p1 = 400000 - 1669
    p2 = 400000 - 1669 - 65 - 1845.8 - 65 - 1062
    t1 = 305.15
    t2 = 328.35

    p1t = p1 - 65
    p2t = p2 + 65

    q = 0.10425

    Pm = 0.5 * (p1t + p2t)
    Tm = 0.5 * (t1 + t2)

    network_data = {
        'name': 'network_test_01',
        'kind': 'network',
        'submodels': {
            'node_A': {
                'kind': 'node',
                'module': 'models.source',
                'class': 'Source',
                'specifications': {
                    'P': p1,
                },
                'parameters': {
                    'Text': t1,
                    'Pext': p1,
                    'x': 0.0,
                    'y': 0.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': t1,
                    'P': p1,
                }
            },
            'node_B': {
                'kind': 'node',
                'module': 'models.sink',
                'class': 'Sink',
                'specifications': {
                    'P': p2,
                },
                'parameters': {
                    'Text': t2,
                    'Pext': p2,
                    'x': 0.0,
                    'y': 1.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': t2,
                    'P': p2,
                }
            },
            'pipe_01': {
                'kind': 'edge',
                'module': 'models.biofilmed_fixed_external_convection_pipe',
                'class': 'BiofilmedFixedExternalConvectionPipe',
                'from': 'node_A',
                'to': 'node_B',
                "domains": {
                    "x": {
                        "initial": 0.0,
                        "final": 1.0,
                        "N": 10
                    },
                },
                'states': {
                },
                'parameters': {
                    'Di': 0.01575,
                    'Do': 0.01905,
                    'kwall': 16.0,
                    'tetha': 0.0,
                    'L': 6.1,
                    'ep': 45e-6,
                    'Klb': 0.5 * 0.9,
                    'Kub': 0.5 * 0.9,
                    'Npipes': 1.0,
                    'Text': 273.15 + 60.,
                    'hext': 17284.56,
                    'mfi': 1e-6,
                    'rhomf': 980.,
                },
                'specifications': {

                },
                'initial_guess': {
                    'Re': 13793,
                    'D': 0.01575,
                    'v': 0.54,
                    'k': q,
                    'T':  {
                      "initial": t1,
                      "final": t2
                    },
                    'P':  {
                      "initial": p1t,
                      "final": p2t
                    },
                    'fD': 0.0330,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                    'Qout': {
                        "initial": -3526.5,
                        "final": -605.2
                    },
                    'To': {
                        "initial": 329.75,
                        "final": 332.57
                    },
                    'Ti': {
                        "initial": 323.07,
                        "final": 331.42
                    },
                    'hint': 3976.69,
                    'Resistance': 0.007941,

                    'Tbf': t1,
                    'vbf': 0.54,
                },
            },
        },
    }

    return network_data


def case_biofilmed_external_film_cond_pipe():


    p1 = 400000 - 1669
    p2 = 400000 - 1669 - 65 - 1845.8 - 65 - 1062
    t1 = 305.15
    t2 = 328.35

    p1t = p1 - 65
    p2t = p2 + 65

    q = 0.10425

    Pm = 0.5 * (p1t + p2t)
    Tm = 0.5 * (t1 + t2)

    network_data = {
        'name': 'network_test_01',
        'class': 'Network',
        'submodels': {
            'node_A': {
                'kind': 'node',
                'module': 'models.source',
                'class': 'Source',
                'specifications': {
                    'P': p1,
                },
                'parameters': {
                    'Text': t1,
                    'Pext': p1,
                    'x': 0.0,
                    'y': 0.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': t1,
                    'P': p1,
                }
            },
            'node_B': {
                'kind': 'node',
                'module': 'models.sink',
                'class': 'Sink',
                'specifications': {
                    'P': p2,
                },
                'parameters': {
                    'Text': t2,
                    'Pext': p2,
                    'x': 0.0,
                    'y': 1.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': t2,
                    'P': p2,
                }
            },
            'pipe_01': {
                'kind': 'edge',
                'module': 'models.biofilmed_external_film_cond_pipe',
                'class': 'BiofilmedExternalFilmCondPipe',
                'from': 'node_A',
                'to': 'node_B',
                "domains": {
                    "x": {
                        "initial": 0.0,
                        "final": 1.0,
                        "N": 10
                    },
                },
                'states': {
                },
                'parameters': {
                    'Di': 0.01575,
                    'Do': 0.01905,
                    'kwall': 16.0,
                    'tetha': 0.0,
                    'L': 6.1,
                    'ep': 45e-6,
                    'Klb': 0.5 * 0.9,
                    'Kub': 0.5 * 0.9,
                    'Npipes': 1.0,
                    'rhov':  0.13044 ,
                    'hvap': 2357680,
                    'Tsat': 60 + 273.15,
                    'Psat': 19950,
                    'fNtub': 1.0,
                    'rhomf': 980.,
                    'mfi': 1e-6,
                },
                'specifications': {

                },
                'initial_guess': {
                    'Re': 13793,
                    'D': 0.01575,
                    'v': 0.54,
                    'k': q,
                    'T':  {
                      "initial": t1,
                      "final": t2
                    },
                    'P':  {
                      "initial": p1t,
                      "final": p2t
                    },
                    'fD': 0.0330,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                    'Qout': {
                        "initial": -3526.5,
                        "final": -605.2
                    },
                    'To': {
                        "initial": 329.75,
                        "final": 332.57
                    },
                    'Ti': {
                        "initial": 323.07,
                        "final": 331.42
                    },
                    'hint': 3976.69,
                    'Resistance': 0.007941,
                    'Tbf': {
                        "initial": t1 + 1,
                        "final": t2 + 1
                    },
                    'vbf': 0.54,
                },
            },
        },
    }

    return network_data


def case_external_film_condensation_tube_arrange():

    p1 = 400000 - 1669
    p2 = 400000 - 1669 - 65 - 1845.8 - 65 - 1062
    t1 = 305.15
    t2 = 328.35

    p1t = p1 - 65
    p2t = p2 + 65

    q = 0.10425

    Pm = 0.5 * (p1t + p2t)
    Tm = 0.5 * (t1 + t2)

    network_data = {
        'name': 'network_test_01',
        'class': 'Network',
        'submodels': {
            'node_A': {
                'kind': 'node',
                'module': 'models.source',
                'class': 'Source',
                'specifications': {
                    'P': p1,
                },
                'parameters': {
                    'Text': t1,
                    'Pext': p1,
                    'x': 0.0,
                    'y': 0.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': t1,
                    'P': p1,
                }
            },
            'node_B': {
                'kind': 'node',
                'module': 'models.sink',
                'class': 'Sink',
                'specifications': {
                    'P': p2,
                },
                'parameters': {
                    'Text': t2,
                    'Pext': p2,
                    'x': 0.0,
                    'y': 1.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': t2,
                    'P': p2,
                }
            },
            'pipe_01': {
                'kind': 'edge',
                'module': 'models.external_film_condensation_tube_arrange',
                'class': 'ExternalFilmCondensationTubeArrange',
                'from': 'node_A',
                'to': 'node_B',
                "domains": {
                    "x": {
                        "initial": 0.0,
                        "final": 1.0,
                        "N": 10
                    },
                    "y": {
                        "N": 2,
                    },
                },
                'states': {
                },
                'parameters': {
                    'Di': 0.01575,
                    'Do': 0.01905,
                    'kwall': 16.0,
                    'tetha': 0.0,
                    'L': 6.1,
                    'ep': 45e-6,
                    'Klb': 0.5 * 0.9,
                    'Kub': 0.5 * 0.9,
                    'Npipes': [1.0, 1.0],
                    'rhov':  0.13044 ,
                    'hvap': 2357680,
                    'Tsat': 60 + 273.15,
                    'Psat': 19950,
                    'fNtub': [1.0, 1.0],
                },
                'specifications': {

                },
                'initial_guess': {
                    'Re': 13793,
                    'D': 0.01575,
                    'v': 0.54,
                    'k': q,
                    'T':  {
                      "initial": t1,
                      "final": t2
                    },
                    'P':  {
                      "initial": p1t,
                      "final": p2t
                    },
                    'fD': 0.0330,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                    'Qout': {
                        "initial": -3526.5,
                        "final": -605.2
                    },
                    'To': {
                        "initial": 329.75,
                        "final": 332.57
                    },
                    'Ti': {
                        "initial": 323.07,
                        "final": 331.42
                    },
                    'hint': 3976.69,
                    'Resistance': 0.007941,
                },
            },
        },
    }

    return network_data





def case_biofilmed_external_film_cond_tube_arrange():


    p1 = 400000 - 1669
    p2 = 400000 - 1669 - 65 - 1845.8 - 65 - 1062
    t1 = 305.15
    t2 = 328.35

    p1t = p1 - 65
    p2t = p2 + 65

    q = 0.10425

    Pm = 0.5 * (p1t + p2t)
    Tm = 0.5 * (t1 + t2)

    network_data = {
        'name': 'network_test_01',
        'class': 'Network',
        'submodels': {
            'node_A': {
                'kind': 'node',
                'module': 'models.source',
                'class': 'Source',
                'specifications': {
                    'P': p1,
                },
                'parameters': {
                    'Text': t1,
                    'Pext': p1,
                    'x': 0.0,
                    'y': 0.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': t1,
                    'P': p1,
                }
            },
            'node_B': {
                'kind': 'node',
                'module': 'models.sink',
                'class': 'Sink',
                'specifications': {
                    'P': p2,
                },
                'parameters': {
                    'Text': t2,
                    'Pext': p2,
                    'x': 0.0,
                    'y': 1.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': t2,
                    'P': p2,
                }
            },
            'pipe_01': {
                'kind': 'edge',
                'module': 'models.biofilmed_external_film_cond_tube_arrange',
                'class': 'BiofilmedExternalFilmCondensationTubeArrange',
                'from': 'node_A',
                'to': 'node_B',
                "domains": {
                    "x": {
                        "initial": 0.0,
                        "final": 1.0,
                        "N": 10
                    },
                    "y": {
                      "N": 2
                    }
                },
                'states': {
                },
                'parameters': {
                    'Di': 0.01575,
                    'Do': 0.01905,
                    'kwall': 16.0,
                    'tetha': 0.0,
                    'L': 6.1,
                    'ep': 45e-6,
                    'Klb': 0.5 * 0.9,
                    'Kub': 0.5 * 0.9,
                    'Npipes': [1, 1],
                    'rhov':  0.13044 ,
                    'hvap': 2357680,
                    'Tsat': 60 + 273.15,
                    'Psat': 19950,
                    'fNtub': [1.0, 1.0],
                    'rhomf': 980.,
                    'mfi': 1e-6,
                },
                'specifications': {

                },
                'initial_guess': {
                    'Re': 13793,
                    'D': 0.01575,
                    'v': 0.54,
                    'k': q,
                    'T':  {
                      "initial": t1,
                      "final": t2
                    },
                    'P':  {
                      "initial": p1t,
                      "final": p2t
                    },
                    'fD': 0.0330,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                    'Qout': {
                        "initial": -3526.5,
                        "final": -605.2
                    },
                    'To': {
                        "initial": 329.75,
                        "final": 332.57
                    },
                    'Ti': {
                        "initial": 323.07,
                        "final": 331.42
                    },
                    'hint': 3976.69,
                    'Resistance': 0.007941,
                    'Tbf': {
                        "initial": t1 + 1,
                        "final": t2 + 1
                    },
                    'vbf': 0.54,
                },
            },
        },
    }

    return network_data
