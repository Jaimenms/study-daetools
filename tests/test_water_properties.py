import water_properties


def test_density():

    T = 300

    P = 100000

    rho = water_properties.density(T, P)

    print("----", rho)

    assert False


def test_enthalpy():

    T = 300

    P = 100000

    h = water_properties.enthalpy(T, P)

    print("----", h)

    assert False


def test_heat_capacity():

    T = 300

    P = 100000

    cp = water_properties.heat_capacity(T, P)

    print("----", cp)

    assert False


def test_conductivity():

    T = 300

    P = 100000

    kappa = water_properties.conductivity(T, P)

    print("----", kappa)

    assert False


def test_viscosity():

    T = 300

    P = 100000

    mu = water_properties.viscosity(T, P)

    print("----", mu)

    assert False