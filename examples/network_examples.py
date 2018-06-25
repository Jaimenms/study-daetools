__doc__ = """
This file concentrates all the test cases for the network simulation. So consider te functions beginning with the
name case_ and returning the data dictionary.

"""

from water_properties import density, heat_capacity, conductivity, viscosity


def case_pipe():

    p1 = 9.95*1e5
    p2 = 985422.9428151655
    q = 8.181359

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
                    'tetha': 0.0,
                    'L': 100.0,
                    'ep': 0.0018*0.0254,
                    'Klb': 0.2,
                    'Kub': 0.2, 'Npipes': 1.0,
                },
                'specifications':{

                },
                'initial_guess': {
                    'Re': 187515,
                    'D': 1.0 * 4.026*0.0254,
                    'v': 1.5,
                    'k': q,
                    'T': 300.0,
                    'P': [p1, p2],
                    'fD': 0.018,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                },
            },
        },
    }


    return network_data


def case_fixed_external_convection_pipe():

    p1 = 9.95*1e5
    p2 = 985422.9428151655
    q = 8.18375730036397

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
                'module': 'models.fixed_external_convection_pipe',
                'class': 'FixedExternalConvectionPipe',
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
                    'Text': 273 + 60.,
                    'hext' : 6873,
                    'Klb': 0.2,
                    'Kub': 0.2, 'Npipes': 1.0,
                },
                'specifications':{

                },
                'initial_guess': {
                    'Re': 118786.29183227483,
                    'D': 1.0 * 4.026*0.0254,
                    'v': 1.0,
                    'k': q,
                    'T': 300.0,
                    'Ti': 322,
                    'To': 328,
                    'P': [p1, p2],
                    'fD': 0.019647235607534896,
                    'Qout': -20218.08370826369,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                    'hint': 1398.8091902015963,

                },
            },
        },
    }


    return network_data


def case_fixed_external_temperature_pipe():

    p1 = 9.95*1e5
    p2 = 985422.9428151655
    q = 8.181359

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
                'module': 'models.fixed_external_temperature_pipe',
                'class': 'FixedExternalTemperaturePipe',
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
                    'Text': 273 + 60.,
                    'Klb': 0.2,
                    'Kub': 0.2, 'Npipes': 1.0,
                },
                'specifications':{

                },
                'initial_guess': {
                    'Re': 187515,
                    'D': 1.0 * 4.026*0.0254,
                    'v': 1.0,
                    'k': q,
                    'T': 300.0,
                    'Ti': 326.13475139707714,
                    'To': 333.0,
                    'P': [p1, p2],
                    'fD': 0.019647,
                    'Qout': -23477.07466549454,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                    'hint': 1398.8091902015963,
                },
            },
        },
    }


    return network_data




def case_external_film_condensation_pipe():

    p1 = 9.95*1e5
    p2 = 985422.9428151655
    q = 8.181359
    Ntub = 1

    t1 = 300.0
    t2 = 301.0
    Pm = 0.5 * (p1 + p2)
    Tm = 0.5 * (t1 + t2)
    Psat = 0.2*1e5

    import water_at_saturation_properties

    Tsat = water_at_saturation_properties.calculate_saturation_temperature(Psat)

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
                'module': 'models.external_film_condensation_pipe',
                'class': 'ExternalFilmCondensationPipe',
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
                    'rhov': 1/7.650,
                    'hvap': 2358400,
                    'Tsat': Tsat,
                    'Psat': Psat,
                    'fNtub': Ntub ** 0.75 - (Ntub - 1) ** 0.75,
                    'Klb': 0.2,
                    'Kub': 0.2, 'Npipes': 1.0,

                },
                'specifications':{

                },
                'initial_guess': {
                    'Re': 187515,
                    'D': 1.0 * 4.026*0.0254,
                    'v': 1.0,
                    'k': q,
                    'T': 300.0,
                    'Ti': 322.5068326854783,
                    'To': 328.41907571019976,
                    'P': [p1, p2],
                    'fD': 0.018,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                    'hint': 1425.,
                    'hext': 10000.,
                    'Qout': -21000.
                },
            },
        },
    }


    return network_data


def case_biofilmed_fixed_external_convection_pipe():

    p1 = 9.95 * 1e5
    p2 = 985422.9428151655
    q = 8.18375730036397

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
                'initial_guess': {
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
                'initial_guess': {
                    'w': q,
                    'T': 300.0,
                    'P': p2,
                }
            },
            'pipe_01': {
                'kind': 'edge',
                'module': 'models.biofilmed_fixed_external_convection_pipe',
                'class': 'BiofilmedFixedExternalConvectionPipe',
                'from': 'node_A',
                'to': 'node_B',
                'domains': {
                    'N': 10,
                },
                'states': {
                    'stnWaterPropertiesSetup': 'Fixed',  # Variable or Fixed
                },
                'parameters': {
                    'Di': 1.0 * 4.026 * 0.0254,
                    'Do': 1.1 * 4.026 * 0.0254,
                    'tetha': 0.0,
                    'L': 100.0,
                    'ep': 0.0018 * 0.0254,
                    'kwall': 51.9,
                    'Text': 273 + 60.,
                    'hext': 6873,
                    'rhomf': 980.,
                    'mfi': 1e-9,
                    'Klb': 0.2,
                    'Kub': 0.2, 'Npipes': 1.0,

                },
                'specifications': {

                },
                'initial_guess': {
                    'Re': 118786.29183227483,
                    'D': 1.0 * 4.026 * 0.0254,
                    'v': 1.0,
                    'k': q,
                    'T': 300.0,
                    'Ti': 322,
                    'To': 328,
                    'P': [p1, p2],
                    'fD': 0.019647235607534896,
                    'Qout': -20218.08370826369,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                    'hint': 1398.8091902015963,
                    'mf': 1e-9,
                },
                'initial_conditions': {
                    #'mf': 0.0,
                },
            },
        },
    }

    return network_data


def case_biofilmed_pipe():

    p1 = 9.95 * 1e5
    p2 = 985422.9428151655
    q = 8.18375730036397

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
                'initial_guess': {
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
                'initial_guess': {
                    'w': q,
                    'T': 300.0,
                    'P': p2,
                }
            },
            'pipe_01': {
                'kind': 'edge',
                'module': 'models.biofilmed_pipe',
                'class': 'BiofilmedPipe',
                'from': 'node_A',
                'to': 'node_B',
                'domains': {
                    'N': 10,
                },
                'states': {
                    'stnWaterPropertiesSetup': 'Fixed',  # Variable or Fixed
                },
                'parameters': {
                    'Di': 1.0 * 4.026 * 0.0254,
                    'tetha': 0.0,
                    'L': 100.0,
                    'ep': 0.0018 * 0.0254,
                    'rhomf': 980.,
                    'mfi': 1e-9,
                    'Klb': 0.2,
                    'Kub': 0.2, 'Npipes': 1.0,

                },
                'specifications': {

                },
                'initial_guess': {
                    'Re': 118786.29183227483,
                    'D': 1.0 * 4.026 * 0.0254,
                    'v': 1.0,
                    'k': q,
                    'T': 300.0,
                    'P': [p1, p2],
                    'fD': 0.019647235607534896,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                    'mf': 1e-9,
                },
                'initial_conditions': {
                    #'mf': 0.0,
                },
            },
        },
    }

    return network_data


def case_biofilmed_external_film_cond_pipe():

    p1 = 9.95*1e5
    p2 = 985422.9428151655
    q = 8.181359
    Ntub = 1

    t1 = 300.0
    t2 = 301.0
    Pm = 0.5 * (p1 + p2)
    Tm = 0.5 * (t1 + t2)
    Psat = 0.2*1e5

    import water_at_saturation_properties

    Tsat = water_at_saturation_properties.calculate_saturation_temperature(Psat)

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
                'module': 'models.biofilmed_external_film_cond_pipe',
                'class': 'BiofilmedExternalFilmCondensationPipe',
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
                    'rhov': 1/7.650,
                    'hvap': 2358400,
                    'Tsat': Tsat,
                    'Psat': Psat,
                    'fNtub': Ntub ** 0.75 - (Ntub - 1) ** 0.75,
                    'rhomf': 980.,
                    'mfi': 1e-9,
                    'Klb': 0.2,
                    'Kub': 0.2, 'Npipes': 1.0,

                },
                'specifications':{

                },
                'initial_guess': {
                    'Re': 187515,
                    'D': 1.0 * 4.026*0.0254,
                    'v': 1.0,
                    'k': q,
                    'T': 300.0,
                    'Ti': 322.5068326854783,
                    'To': 328.41907571019976,
                    'P': [p1, p2],
                    'fD': 0.018,
                    'kappa': conductivity(Tm, Pm, simplified=True),
                    'cp': heat_capacity(Tm, Pm, simplified=True),
                    'rho': density(Tm, Pm, simplified=True),
                    'mu': viscosity(Tm, Pm, simplified=True),
                    'hint': 1425.,
                    'hext': 10000.,
                    'Qout': -21000.,
                    'mf': 1e-9,
                },
            },
        },
    }


    return network_data
