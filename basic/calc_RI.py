import numpy as np
import matplotlib.pyplot as plt
import plasma_params as pp

def calc_RI(N, B, f):
    theta = 0
    fc = (pp.Q*B/pp.ME)/2/np.pi
    fp = ((pp.Q**2*N)/(pp.EPS*pp.ME))**(1/2)/2/np.pi
    mu = fp*(fc*f*(np.cos(theta)-(f/fc)))**(-1/2)
    return mu

def calc_RI_plot(N, B, f):
    theta = np.linspace(-np.pi/2, np.pi/2, 4001)
    fc = (pp.Q*B/pp.ME)/2/np.pi
    fp = ((pp.Q**2*N)/(pp.EPS*pp.ME))**(1/2)/2/np.pi
    mu = fp*(fc*f*(np.cos(theta)-(f/fc)))**(-1/2)
    return theta, mu

def calc_RI_2(N1, B1, f, Nrate, Brate):

    N1 = N1*Nrate

    B1 = B1*Brate

    
    theta = np.linspace(-np.pi/2, np.pi/2, 501)
    fc = (pp.Q*B1/pp.ME)/2/np.pi
    fp = ((pp.Q**2*N1)/(pp.EPS*pp.ME))**(1/2)/2/np.pi
    mu = fp*(fc*f*(np.cos(theta)-(f/fc)))**(-1/2)
    x = mu*np.sin(theta)
    y = mu*np.cos(theta)
    return x,y

def calc_RI_3(PoC, f, B):
    theta = np.linspace(-np.pi/2, np.pi/2, 500)
    fc = (pp.Q*B/pp.ME)/2/np.pi
    mu = np.outer(PoC, (f*(fc**-2)*(fc*np.cos(theta)-f))**(-1/2)) 
    x = mu*np.sin(theta)
    y = mu*np.cos(theta)
    return x,y

def calc_RI_4(PoC, f):
    theta = 0
    mu = np.outer(PoC, (f*((f/0.2)**-2)*((f/0.2)*np.cos(theta)-f))**(-1/2)) 
    
    x = mu*np.sin(theta)
    y = mu*np.cos(theta)
    return x,y