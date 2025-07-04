from scipy.stats import f_oneway

groups = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
alpha = 0.05

F, p = f_oneway(groups[0], groups[1], groups[2])
print(F, p)
