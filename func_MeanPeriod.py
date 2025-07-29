
import sys
import csv
import numpy as np
import math
import time

import sobol_seq

from scipy import integrate
#from scipy.integrate import quad
#from scipy.special import i0
#from scipy.stats import rice





def NoticeMessages():
    print('<< Notice >>')
    print('The list argument "param_sets" must contain the following 10 elements in its list.')
    print('the 10 elements -> 1 : skipping time for ppp')
    print('                   2 : skipping time for tcp')
    print('                   3 : moving velocity of a user')
    print('                   4 : intensity for ppp')
    print('                   5 : parent intensity for tcp')
    print('                   6 : daughter intensity for tcp')
    print('                   7 : daighter variance for tcp')
    print('                   8 : transmitting power for ppp (macro base stations)')
    print('                   9 : transmitting power for tcp (small base stations)')
    print('                   10: path-loss exponent')


def MeanPeriod(param_sets, DivNum=100):

    ### parameter input
    #s1, s2, veloc, lam1, lam2p, mb, sigma, P0, P1, beta = param_sets
    s1, s2, veloc, lam1, lam2p, mb, rd, P0, P1, beta = param_sets


    def P_bar(tier_i, tier_j):
        def P(tier_i):
            if tier_i == 1:
                return P0
            elif tier_i == 2:
                return P1
            else:
                print('Error;def P(tier_i): tier_i must be either 1 or 2.')
                sys.exit()
        return (P(tier_j)/P(tier_i))**(1/beta)

    #def func_fd_M(x, z):
    #    return rice.pdf(x, z/sigma, scale=sigma)
    #def func_Fd_M(r, z): 
    #    return rice.cdf(r, z/sigma, scale=sigma)
    def func_fd_M(r, z):
        if r<=max(rd - z, 0):
            return 2*r/rd**2
        elif abs(rd - z)<=r and r<=rd + z: 
            return 1/np.pi*np.arccos( round((r**2 + z**2 - rd**2)/(2*r*z), 8) )*2*r/rd**2
        else:
            return 0.0
    def func_Fd_M(r, z):
        def integd_x(L_x, z):
            return np.array([ x*np.arccos( round((x**2 + z**2 - rd**2)/(2*x*z), 8) ) for x in L_x])

        if min(r, abs(rd - z)) == min(r, rd + z):
            return (min(r, max(rd - z, 0))**2)/rd**2
        else:
            L_x = np.linspace(min(r, abs(rd - z)), min(r, rd + z), DivNum+1)
            return (min(r, max(rd - z, 0))**2 + 2/np.pi*integrate.trapz(integd_x(L_x, z), L_x))/rd**2
    

    def func_eta(tier_i, r):    ## same as func_A in main_cDRHOR.py
        def integd_z_cv(L_z_cv, r):
            def integd_z(z, r):
                return z*(1 - np.exp( -mb*func_Fd_M(P_bar(tier_i, 2)*r, z) ))
            L_z = np.tan(np.pi/2*L_z_cv)
            return np.array([ np.pi/2*(1 + z**2)*integd_z(z, r) for z in L_z ])
        L_z_cv = np.linspace(0, 1, DivNum+1)
        return np.exp( -2*np.pi*lam2p * integrate.trapz(integd_z_cv(L_z_cv, r), L_z_cv) )

    def func_mu(r):    ## same as func_B in main_cDRHOR.py
        def integd_z_cv(L_z_cv, r):
            def integd_z(z, r):
                return z*func_fd_M(r, z)*np.exp( -mb*func_Fd_M(r, z) )
            L_z = np.tan(np.pi/2*L_z_cv)
            return np.array([ np.pi/2*(1 + z**2)*integd_z(z, r) for z in L_z ])

        L_z_cv = np.linspace(10**(-8), 1, DivNum+1)
        return 2*np.pi*lam2p * integrate.trapz(integd_z_cv(L_z_cv, r), L_z_cv)


    def func_Aprob_ppp():
        def integd_r_cv(L_r_cv):
            def integd_r(r):
                return r*np.exp(-np.pi*lam1*r**2) * func_eta(1, r)
            L_r = np.tan(np.pi/2*L_r_cv)
            return np.array([ np.pi/2*(1 + r**2)*integd_r(r) for r in L_r ])

        L_r_cv = np.linspace(10**(-8), 1, DivNum+1)
        return 2*np.pi*lam1*integrate.trapz(integd_r_cv(L_r_cv), L_r_cv)

    def func_Aprob_pcp():
        def integd_r_cv(L_r_cv):
            def integd_r(r):
                return np.exp( -np.pi*lam1*P_bar(2, 1)**2*r**2 ) * func_eta(2, r) * func_mu(r)
            L_r = np.tan(np.pi/2*L_r_cv)
            return np.array([ np.pi/2*(1 + r**2)*integd_r(r) for r in L_r ])
        L_r_cv = np.linspace(10**(-8), 1, DivNum+1)
        return mb*integrate.trapz(integd_r_cv(L_r_cv), L_r_cv)


    Aprob_ppp = func_Aprob_ppp()
    Aprob_pcp = func_Aprob_pcp()

    ### test start ###
#    print('--- func_MeanPeriod_v2 ---')
#    print('Aprob_ppp = {}, Aprob_pcp = {}, Aprob_ppp + Aprob_pcp = {}'.format(Aprob_ppp, Aprob_pcp, Aprob_ppp + Aprob_pcp))
    ### test end ###


    return s1*Aprob_ppp + s2*Aprob_pcp
  
  
