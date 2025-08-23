from scipy.stats import f_oneway
from scipy.stats import tukey_hsd

groups = [[1,2,3,4,5], [6,7,8,9,10], [11,12,13,14,15]]

F, p = f_oneway(groups[0], groups[1], groups[2])
print(F, p)

tukey_res = tukey_hsd(groups[0], groups[1], groups[2])
print(tukey_res)
