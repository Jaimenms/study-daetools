__doc__ = """
This file concentrates all the test cases for the network simulation. So consider te functions beginning with the
name case_ and returning the data dictionary.

"""

from water_properties import density, heat_capacity, conductivity, viscosity

def case_blank():

    network_data = {}

    return network_data


def case_000():
    p1 = 9.95 * 1e5
    p2 = 9.75 * 1e5
    q = 10.0

    t1 = 300.0
    t2 = 301.0
    Pm = 0.5 * (p1 + p2)
    Tm = 0.5 * (t1 + t2)

    network_data = {
        'name': 'network_test_01',
        'kind': 'network',
        'module': 'models.network',
        'class': 'Network',
        'submodels': {
            'node_A': {
                'kind': 'node',
                'module': 'models.node',
                'class': 'Source',
                'specifications': {
                    'P': p1,
                },
                'parameters': {
                    'Text': t1,
                    'cpext': heat_capacity(t1, p1, simplified=True),
                    'x': 0.0,
                    'y': 0.0,
                    'z': 0.0
                },
                'initial_guess': {
                    'w': q,
                    'T': 300.0,
                    'P': p1,
                }
            },
        }
    }

    return network_data



def case_00():

    p1 = 9.95*1e5
    p2 = 9.75*1e5
    q = 10.0

    t1 = 300.0
    t2 = 301.0
    Pm = 0.5 * (p1 + p2)
    Tm = 0.5 * (t1 + t2)

    network_data = {
        'name': 'network_test_01',
        'kind': 'network',
        'module': 'models.network',
        'class': 'Network',
        'submodels': {
            'node_A': {
                'kind': 'node',
                'module': 'models.node',
                'class': 'Source',
                'specifications': {
                    'P': p1,
                },
                'parameters': {
                    'Text': t1,
                    'cpext': heat_capacity(t1, p1, simplified=True),
                    'x': 0.0,
                    'y': 0.0,
                    'z': 0.0
                },
                'initial_guess':{
                    'w': q,
                    'T': 300.0,
                    'P': p1,
                }
            },
            'node_B': {
                'kind': 'node',
                'module': 'models.node',
                'class': 'Sink',
                'specifications': {
                    'P': p2,
                },
                'parameters': {
                    'Text': t2,
                    'cpext': heat_capacity(t2, p2, simplified=True),
                    'x': 0.0,
                    'y': 1.0,
                    'z': 0.0
                },
                'initial_guess':{
                    'w': q,
                    'T': 300.0,
                    'P': p2,
                }
            }
        },
    }


    return network_data


def case_01():

    p1 = 9.95*1e5
    p2 = 9.75*1e5
    q = 10.0

    t1 = 300.0
    t2 = 301.0
    Pm = 0.5 * (p1 + p2)
    Tm = 0.5 * (t1 + t2)

    network_data = {
        'name': 'network_test_01',
        'kind': 'network',
        'module': 'models.network',
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
                'initial_guess':{
                    'w': q,
                    'T': 300.0,
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
                    'T': 300.0,
                    'P': p2,
                }
            },
            'pipe_01': {
                'kind': 'edge',
                'module': 'models.pipe',
                'class': 'Pipe',
                'from': 'node_A',
                'to': 'node_B',
                'domains': {
                    'N': 10,
                },
                'states': {
                    'stnWaterPropertiesSetup': 'Fixed', # Variable or Fixed
                },
                'parameters': {
                    'Di': 1.0 * 4.026*0.0254,
                    'Do': 1.1 * 4.026*0.0254,
                    'tetha': 0.0,
                    'L': 100.0,
                    'ep': 0.0018*0.0254,
                    'kwall': 51.9,
                    'Text': Tm,
                    'hext' : 2.5,
                    'rhomf': 980.,
                },
                'specifications':{

                },
                'initial_guess': {
                    'mf': 0.0,
                    'Re': 187515,
                    'D': 1.0 * 4.026*0.0254,
                    'v': 1.5,
                    'k': q,
                    'T': 300.0,
                    'Tw': 300.0,
                    'P': [p1, p2],
                    'fD': 0.018,
                    'Qout': 0.0,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                },
            },
        },
    }


    return network_data

