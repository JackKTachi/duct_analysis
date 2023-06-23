

import numpy as np
import matplotlib.pyplot as plt
import plasma_params as pp
import calc_RI as RI


def plot_RI(N, B, f, Nrate, Brate):
    thetaBp, muBp = RI.calc_RI(N, B, f, Nrate=1, Brate=1+Brate)
    thetaNp, muNp = RI.calc_RI(N, B, f, Nrate=1+Nrate, Brate = 1)
    thetaBm, muBm = RI.calc_RI(N, B, f, Nrate=1, Brate=1-Brate)
    thetaNm, muNm = RI.calc_RI(N, B, f, Nrate=1-Nrate, Brate = 1)
    theta, mu = RI.calc_RI(N, B, f, Nrate=1, Brate = 1)

    fig = plt.figure(figsize=(8,8))

    ax = fig.add_subplot(111, projection='polar')
    ax.plot(thetaBp, muBp, label='B+')
    ax.plot(thetaNp, muNp, label='N+')
    ax.plot(thetaBm, muBm, label='B-')
    ax.plot(thetaNm, muNm, label='N-')
    ax.plot(theta, mu, label='standard')


    ax.set_thetamin(-90)
    ax.set_thetamax(90)
    ax.set_rmin(0)
    ax.set_rmax(20)
    ax.set_theta_zero_location('N')
    ax.legend()
    plt.show()    

