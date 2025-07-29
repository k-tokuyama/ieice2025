
import sys
import csv
import numpy as np
import time

from scipy import integrate
from scipy import stats

## module import ##
#import func_MeanPeriod as f_MP
import func_MeanPeriod as f_MPv2


### Caution ###
# The option PCP_id = 'TCP' does not work. Only PCP_id = 'MCP' is available.




## description of the arguments:
    ## param_sets : the parameter sets of the system parameters of our model.
        ## s1/s2    : the skipping time for the 1st/2nd (macro/small) tier.
        ## veloc    : the velocity of the moving UE.
        ## lam1     : the intensity of the PPP in the 1st (macro) tier.
        ## lam2p/mb : the intensity of the parent/daughter points of the PCP in the 2nd (small) tier.
        ## sigma/rd : the clustering variance/the cluster radius for the TCP/MCP in the 2nd (small) tier.
        ## P1/P2    : the transmission power in the 1st/2nd (macro/small) tier.
        ## beta     : the common path-loss exponent in the 1st/2nd (macro/small) tier.
    ## PCP_id  : the indicator of the specific PCP [TCP/MCP, TCP -> Thomas Cluster Process, MCP -> Matern Cluster Process].



#### main function ####

def NoticeMessages():
    print('<< Notice >>')
    print('11 arguments are required to be written in the input csv.')
    print('the 11 args -> 1 : skipping time for ppp')
    print('               2 : skipping time for pcp')
    print('               3 : moving velocity of a user')
    print('               4 : intensity for ppp')
    print('               5 : parent intensity for pcp')
    print('               6 : daughter intensity for pcp')
    print('               7 : daighter variance for pcp')
    print('               8 : transmitting power for ppp')
    print('               9 : transmitting power for pcp')
    print('               10: path-loss exponent')
    print('               11: curve type <<must be either of "straight"/"circle"/"spiral">>  ')



if __name__ == "__main__":
    argvs = sys.argv

    if not len(argvs[1:]) == 2:
        print('InputError: 2 arguments need to be input.')
        print('                1st: csv filename')
        print('                2nd: PCP indicator')
        print('                     (TCP(tcp) <- Thomas Cluster Process / (MCP) <- Matern Cluster Process)')
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

    PCP_id = argvs[2].upper()
    if PCP_id not in ['TCP', 'MCP']:
        print('InputError: Please set a valid PCP_id; TCP(tcp)/MCP(mcp).')
        print('            (TCP(tcp) <- Thomas Cluster Process / (MCP) <- Matern Cluster Process)')
        sys.exit()


    ## Make the result csv from the input csv.
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

    ## get the data base matrix (datacsv_DR-per-period) from the directory ./database/
    def get_DataBaseMatrix():
#        ## test start ##
#        DataCsvFilename = 'old__result_calc-DR-approx-{}_'.format(PCP_id) + InputCsvFilename
#        ## test end ##

        DataCsvFilename = 'datacsv_DR-per-period-{}_'.format(PCP_id) + InputCsvFilename
        f = open('database/'+DataCsvFilename, 'r')
        csvreader = csv.reader(f)
        header = next(csvreader)
        matrix = [v for v in csvreader][1:]    ## the first row (t1, t2 = 0.0) is being cut.
        f.close()
        return matrix
    def get_DataBaseMatrixRowIndex(s_ppp, s_pcp):
        index_dic = {float(i): int(i)-1 for i in range(1, 101)}   ## index_dic returns the corresponding index in the db_matrix for the given s_ppp/s_pcp (float type)
        return index_dic[s_ppp], index_dic[s_pcp]




    ## Make matrix data of the result csv.
    db_matrix = get_DataBaseMatrix()
    cid_ppp, cid_pcp = 4, 5    ## column id of the db_matrix (for getting period_DR_ppp/period_DR_pcp)

    result_rows_L = []
    for i in range(len(matrix)):
        ParamSets, CurveType, ErrorFlag = make_parameters_set(matrix[i])

        s_ppp, s_pcp = ParamSets[:2]

        if not s_ppp or not s_pcp:
            pass
            #result_row2 = [None, None]
        else:
            ## get period_DR_ppp/period_DR_pcp from the db_matrix ##
            rid_ppp, rid_pcp = get_DataBaseMatrixRowIndex(s_ppp, s_pcp)    ## row id of the db martix (for getting period_DR_ppp/period_DR_pcp)
            period_DR_ppp = float(db_matrix[rid_ppp][cid_ppp])
            period_DR_pcp = float(db_matrix[rid_pcp][cid_pcp])

            ## calculate mean_period ##
            mean_period = f_MPv2.MeanPeriod(ParamSets, DivNum=200)

            result_row1 = [s_ppp, s_pcp]
            result_row2 = [period_DR_ppp/mean_period, period_DR_pcp/mean_period]
            #result_row2 = (np.sum(np.array(result_rows_L)[:, 2:4], axis=0) / mean_period).tolist()

            result_rows_L.append(result_row1 + result_row2)


#            ## test start ##
#            print('s_ppp={}, s_pcp={}, rid_ppp={}, rid_pcp={}, period_DR_ppp={}, period_DR_pcp={}, mean_period={}'.format(s_ppp, s_pcp, rid_ppp, rid_pcp, period_DR_ppp, period_DR_pcp, mean_period))
#            ## test end ##

    ## make header data of the result csv.
    header1 = ['skipping time1 (macro)', 'skipping time2 (small)']
    header2 = ['mean data rate1 (macro)', 'mean data rate2 (small)']
    result_header = header1 + header2


    ## Output the result csv file.
    OutputCsvFilename = 'result_calc-DR-approx-{}_'.format(PCP_id) + InputCsvFilename
    f = open(OutputCsvFilename, 'w')
    csvwriter = csv.writer(f)
    csvwriter.writerow(result_header)
    csvwriter.writerows(result_rows_L)
    f.close()




