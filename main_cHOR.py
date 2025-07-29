
import sys
import csv
import random as rand
import numpy as np
import math
import time

from scipy import integrate


import sobol_seq

from scipy.integrate import quad
from scipy.special import i0
from scipy.stats import rice

## module import ##
#import func_MeanPeriod_v1 as f_MP
import func_MeanPeriod as f_MPv2


## description of the arguments:
    ## param_sets : the parameter sets of the system parameters of our model.


def output_HOR(param_sets):

    ### parameter input
    s1, s2, veloc, lam1, lam2p, mb, rd, P0, P1, beta = param_sets

#    DivNum1 = 100
#    DivNum2 = 16
#    DivNum_x = 100

    DivNum1 = 200
    DivNum2 = 32
    DivNum_x = 200

    def wcos(r, l, theta):
        return np.sqrt(r**2 + l**2 - 2*r*l*np.cos(theta))
    def R_til(x, h, l):
        if not x:
            x = 10**(-8)
        if not l:
            l = 10**(-8)
        return np.arccos(np.clip((x**2+l**2-h**2)/(2*x*l),-1,1))
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
            L_x = np.linspace(min(r, abs(rd - z)), min(r, rd + z), DivNum1+1)
            return (min(r, max(rd - z, 0))**2 + 2/np.pi*integrate.trapz(integd_x(L_x, z), L_x))/rd**2

    def func_G_M(rr, h, l, z, phi):
        def integd_x1(L_x1, h, l, z):
            return np.array([ x*R_til(x, h, l) for x in L_x1 ])
        def integd_x2(L_x2, h, l, z, phi):
            def func_H_M(x, psi, z, phi):
                R = np.arccos( round((z**2 + x**2 - rd**2)/(2*z*x), 8) )
                return max(phi + psi + R - 2*np.pi, 0) + min(phi + psi, R) - np.clip(phi - psi, -R, R)
            return np.array([ x*func_H_M(x, R_til(x, h, l), z, phi) for x in L_x2 ])

        L_x1 = np.linspace(min(rr, max(rd - z, 0)), min(h + l, max(rd - z, 0)), DivNum_x+1)
        L_x2_hl = np.linspace(min(h + l, abs(rd - z)), min(h + l, rd + z), DivNum_x+1)
        L_x2_rr = np.linspace(min(rr, abs(rd - z)), min(rr, rd + z), DivNum_x+1)
        if min(rr, max(rd - z, 0)) == min(h + l, max(rd - z, 0)):
            term1 = 0.0
        else:
            term1 = integrate.trapz(integd_x1(L_x1, h, l, z), L_x1)
        if min(h + l, abs(rd - z)) == min(h + l, rd + z):
            term2 = 0.0
        else:
            term2 = integrate.trapz(integd_x2(L_x2_hl, h, l, z, phi), L_x2_hl)
        if min(rr, abs(rd - z)) == min(rr, rd + z):
            term3 = 0.0
        else:
            term3 = integrate.trapz(integd_x2(L_x2_rr, h, l, z, phi), L_x2_rr)
        return (term1 + term2 - term3)/(np.pi*rd**2)

    #### out_HOR starts: the function for caluculating HOR (the handover rate) ####

    def func_A(tier_i, r):
        def integd_z_cv(L_z_cv, r):
            def integd_z(z):
                return z*(1 - np.exp( -mb*func_Fd_M(P_bar(tier_i, 2)*r, z) ))
            L_z = np.tan(np.pi/2*L_z_cv)
            return np.array([ np.pi/2*(1 + z**2)*integd_z(z) for z in L_z ])
        L_z_cv = np.linspace(0, 1, DivNum1+1)
        return np.exp( -2*np.pi*lam2p * integrate.trapz(integd_z_cv(L_z_cv, r), L_z_cv) )

    def func_B(r):
        def integd_z_cv(L_z_cv, r):
            def integd_z(z):
                return z*func_fd_M(r, z)*np.exp( -mb*func_Fd_M(r, z) )
            L_z = np.tan(np.pi/2*L_z_cv)
            return np.array([ np.pi/2*(1 + z**2)*integd_z(z) for z in L_z ])
        L_z_cv = np.linspace(10**(-8), 1, DivNum1+1)
        return 2*np.pi*lam2p * integrate.trapz(integd_z_cv(L_z_cv, r), L_z_cv)

    def func_C(r, l, theta):
        def integd_z_cv(L_z_cv, r, l, theta):
            def integd_z(z, r, l, theta):
                def integd_phi(L_phi, z, r, l, theta):
                    return np.array([ np.exp( -mb*func_G_M(r, wcos(r, l, theta), l, z, phi) ) for phi in L_phi ])
                L_phi = np.linspace(0, np.pi, DivNum2+1)
                return z*func_fd_M(r, z)*np.exp(-mb*func_Fd_M(r, z))*integrate.trapz(integd_phi(L_phi, z, r, l, theta), L_phi)
            L_z = np.tan(np.pi/2*L_z_cv)
            return np.array([ np.pi/2*(1 + z**2)*integd_z(z, r, l, theta) for z in L_z ])
        L_z_cv = np.linspace(10**(-8), 1, DivNum1+1)
        return 2*lam2p*integrate.trapz(integd_z_cv(L_z_cv, r, l, theta), L_z_cv)

    def func_D(tier_i, r, l, theta):

        def area_D1(tier_i, r, l, theta):
            def term1(tier_i, r, l, theta):
                def func_eta(r, l, theta):
                    return wcos(r, l, theta)**2*np.arccos( round((r*np.cos(theta) - l)/wcos(r, l, theta), 8) ) - r**2*theta + r*l*np.sin(theta)
                return func_eta(P_bar(tier_i, 1)*r, l, R_til(P_bar(tier_i, 1)*r, P_bar(tier_i, 1)*wcos(r, l, theta), l))
            def term2(tier_i, r, l, theta):
                return np.pi*max(P_bar(tier_i, 1)**2*wcos(r, l, theta)**2 - (l + P_bar(tier_i, 1)*r)**2, 0)
            return term1(tier_i, r, l, theta) + term2(tier_i, r, l, theta)

        def integd_z_cv(L_z_cv, r, l, theta):
            def integd_z(z, r, l, theta):
                def integd_phi(L_phi, z, r, l, theta):
                    return np.array([ 1 - np.exp(-mb*func_G_M(P_bar(tier_i, 2)*r, P_bar(tier_i, 2)*wcos(r, l, theta), l, z, phi)) for phi in L_phi ])
                L_phi = np.linspace(0, np.pi, DivNum2+1)
                return z*np.exp(-mb*func_Fd_M(P_bar(tier_i, 2)*r, z))*integrate.trapz(integd_phi(L_phi, z, r, l, theta), L_phi)
            L_z = np.tan(np.pi/2*L_z_cv)
            return np.array([ np.pi/2*(1 + z**2)*integd_z(z, r, l, theta) for z in L_z ])

        L_z_cv = np.linspace(0, 1, DivNum1+1)
        return np.exp(-lam1*(area_D1(tier_i, r, l, theta)) -2*lam2p*integrate.trapz(integd_z_cv(L_z_cv, r, l, theta), L_z_cv))

    def HOR_term1():

        def integd_r_cv(L_r_cv):
            def integd_r(r):
                def integd_theta(L_theta, r):
                    return np.array([ func_D(1, r, s1*veloc, theta) for theta in L_theta ])
                L_theta = np.linspace(0, np.pi, DivNum2+1)
                return r*np.exp(-lam1*np.pi*r**2) * func_A(1, r) * (1 - integrate.trapz(integd_theta(L_theta, r), L_theta)/np.pi)
            L_r = np.tan(np.pi/2*L_r_cv)
            return np.array([ np.pi/2*(1 + r**2)*integd_r(r) for r in L_r ])
        L_r_cv = np.linspace(10**(-8), 1, DivNum1+1)
        #return 2*np.pi*lam1/s1 * integrate.trapz(integd_r_cv(L_r_cv), L_r_cv)
        return 2*np.pi*lam1 * integrate.trapz(integd_r_cv(L_r_cv), L_r_cv)

    def HOR_term2():

        def integd_r_cv(L_r_cv):
            def integd_r(r):
                def integd_theta(L_theta, r):
                    return np.array([  func_C(r, s2*veloc, theta)*func_D(2, r, s2*veloc, theta) for theta in L_theta ])
                L_theta = np.linspace(0, np.pi, DivNum2+1)
                return np.exp(-lam1*np.pi*P_bar(2, 1)**2*r**2) * func_A(2, r) * (func_B(r) - integrate.trapz(integd_theta(L_theta, r), L_theta)/np.pi)
            L_r = np.tan(np.pi/2*L_r_cv)
            return np.array([ np.pi/2*(1 + r**2)*integd_r(r) for r in L_r ])
        L_r_cv = np.linspace(10**(-8), 1, DivNum1+1)
        #return mb/s2 * integrate.trapz(integd_r_cv(L_r_cv), L_r_cv)
        return mb * integrate.trapz(integd_r_cv(L_r_cv), L_r_cv)


    if s1 and s2:
        start = time.time()
        HORt1 = HOR_term1()
        end = time.time()
        etime1 = end - start

        start = time.time()
        HORt2 = HOR_term2()
        end = time.time()
        etime2 = end - start

#        mean_period = f_MP.MeanPeriod(param_sets, flag_sobol=True)
        mean_period = f_MPv2.MeanPeriod(param_sets, DivNum=DivNum1)
        
        return [s1, s2, HORt1/mean_period, HORt2/mean_period, etime1, etime2]    ## return a line of row
    else: 
        return [s1, s2, None, None, None, None]    ## The case either s1 or s2 is zero is not included in this calculus.




    
#### main function ####

def NoticeMessages():
    print('<< Notice >>')
    print('11 arguments are required to be written in the input csv.')
    print('the 11 args -> 1 : skipping time for ppp')
    print('               2 : skipping time for tcp')
    print('               3 : moving velocity of a user')
    print('               4 : intensity for ppp')
    print('               5 : parent intensity for tcp')
    print('               6 : daughter intensity for tcp')
    print('               7 : daighter variance for tcp')
    print('               8 : transmitting power for ppp (macro base stations)')
    print('               9 : transmitting power for tcp (small base stations)')
    print('               10: path-loss exponent')
    print('               11: curve type <<must be either of "straight"/"circle"/"spiral">>  ')

if __name__ == "__main__":

    argvs = sys.argv

    if not len(argvs[1:]) == 1:
        print('InputError: 1 arguments need to be input.')
        print('                1st: csv filename')
        sys.exit()

    InputCsvFilename = argvs[1]
    if not InputCsvFilename[-4:] == '.csv':
        print('InputError: Please set a csv file (containing 11 parameters) for the 1st argument.')
        NoticeMessages()
        sys.exit()
    else:
        f = open(InputCsvFilename, 'r')
        csvreader = csv.reader(f)
        header = next(csvreader)
        matrix = [v for v in csvreader]
        f.close()


    #outputcsvs_from_inputcsv(matrix)

    def make_parameters_set(matrix_row):
        error_flag = False
        try:
            paras10_L = list(map(float, matrix_row[:10]))
        except ValueError:
            print('\nInput Error: 11 arguments are required to be written in each row of the input csv.')
            print('skip message: the parameter set of the {}th row is omitted due to missing the requirements...\n'.format(i+1))
            paras10_L = None; error_flag = True

        if matrix_row[10].lower() in ['straight', 'circle', 'spiral']:
            curve_type = matrix_row[10].lower()
        else:
            print('\nInput Error: the 11th argument must be either of "straight"/"circle"/"spiral" .')
            print('skip message: the parameter set of the {}th row is omitted due to missing the requirements...\n'.format(i+1))
            curve_type = None; error_flag = True
        return paras10_L, curve_type, error_flag


    result_header = ['skipping time1 (macro)', 'skipping time2 (small)', 'handover rate1 (macro)', 'handover rate2 (small)', 'elapsed time1 (macro)', 'elapsed time2 (small)']
    result_rows_L = []
    NoticeMessagesFlag = False
    for i in range(len(matrix)):
        ParamSets, CurveType, ErrorFlag = make_parameters_set(matrix[i])

        if not ErrorFlag:
            result_row = output_HOR(ParamSets)
            ### test start ###
            print('i = {}, {}'.format(i, result_row))
            ### test end ###
        else:
            NoticeMessagesFlag = True
            result_row = ['', '']
            continue
        result_rows_L.append(result_row)
    if NoticeMessagesFlag:
        NoticeMessages()


    OutputCsvFilename = 'result_calc-HOR-exact-MCP_' + InputCsvFilename
    f = open(OutputCsvFilename, 'w')
    csvwriter = csv.writer(f)
    csvwriter.writerow(result_header)
    csvwriter.writerows(result_rows_L)
    f.close()
