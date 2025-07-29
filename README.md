■ main_cDR.py:
- The function that calculates the expected data rate, shown in eq.(11) of the submitted paper.
- Available for the 2 types of the PPCP classes; Thomas point process (->TCP) and Mattern cluster process (->MCP).
- 1st argument: the system parameter set csv; see "param_set_example1/2/3.csv".
- 2nd argument: the PPCP class indicator; TCP(tcp)/MCP(mcp).
- Run this script AFTER RUNNNING "main_cDR-per-period.py" (otherwise this script does not work).
- Running command example: python main_cDR.py param_set_example1.csv MCP

Required modules:
sys, csv, time, NumPy, SciPy, sobol_seq, func_MeanPeriod (see "func_MeanPeriod.py")

Operation confirmed in Python 3.6.13.




■ main_cDR-per-period.py:
- The function that calculates the expected received data in a cycle, shown in eq.(9) of the submitted paper.
- Available for the 2 types of the PPCP classes; Thomas point process (->TCP) and Mattern cluster process (->MCP).
- 1st argument: the system parameter set csv; see "param_set_example1/2/3.csv".
- 2nd argument: the PPCP class indicator; TCP(tcp)/MCP(mcp).
- Running command example: python main_cDR-per-period.py param_set_example1.csv MCP

Required modules:
sys, csv, time, NumPy, SciPy, sobol_seq, func_MeanPeriod (see "func_MeanPeriod.py")

Operation confirmed in Python 3.6.13.




■ main_cHOR_ver2.py:
- The function that calculates the expected longterm average of the handover rate, shown in eq.(11) of the submitted paper.
- Available for only Mattern cluster process (->MCP).
- 1st argument: the system parameter set csv; see "param_set_example1/2/3.csv".
- Running command example: python main_cHOR_ver2.py param_set_example1.csv

Required modules:
sys, csv, random, time, NumPy, SciPy, sobol_seq, func_MeanPeriod (see "func_MeanPeriod.py")

Operation confirmed in Python 3.6.13.
