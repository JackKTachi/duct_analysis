import numpy as np
from pyspedas.erg.ground.camera.omti_attitude_params import omti_attitude_params
from pyspedas.utilities.time_double import time_double


def ex_omti_attitude_params():

    params = omti_attitude_params()
    print(params)

    params = omti_attitude_params(date='2018-09-12T00:00:00', site='nai')
    print(params)

    params = omti_attitude_params('2015-05-18T00:00:00', 'eur')
    print(params)

    params = omti_attitude_params('2009-01-11/00:00:00', 'TRS')
    print(params)

    print()
    # Following 5 examples should give the same results

    date_double = time_double('2005-01-01 00:00')
    params = omti_attitude_params(date_double, 'Rsb')
    print(params)

    params = omti_attitude_params(int(date_double), 'Rsb')
    print(params)

    params = omti_attitude_params(np.intc(date_double), 'Rsb')
    print(params)

    params = omti_attitude_params(np.float64(date_double), 'Rsb')
    print(params)

    date_double = time_double() # current date and time (> 2005-01-01)
    params = omti_attitude_params(date_double, 'Rsb')
    print(params)


if __name__ == '__main__':
    ex_omti_attitude_params()
