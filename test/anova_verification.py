from scipy.stats import f_oneway

groups = [[1, 2, 1, 0], [2, 3, 2, 1], [10, 7, 10, 8], [5, 4, 5, 6]]

alpha = 0.05

F, p = f_oneway(groups[0], groups[1], groups[2], groups[3])
print(F, p)
