import pandas as pd
import numpy as np
import random
from functions import epanechnikov_kernel, lambda_function, psihat_statistic

random.seed(1)
sigma = 2
T = 3
y_data = np.random.normal(0, sigma, T)
sigmahat = 2

#Creating g_t_set like in Section 2.1
u = np.linspace(1/T, 1, num = T)
h = np.linspace(1/T, 1, num = T)
'''
u = np.linspace(4/T, 1, num = T//4)
h = np.linspace(3/T, 1/4 + 3/T, num = T//20, endpoint = False)
'''
g_t_set_temp = [(x, y, 0) for x in u for y in h]
labels = ['u', 'h', 'values']
g_t_set_t = pd.DataFrame.from_records(g_t_set_temp, columns = labels)
g_t_set = g_t_set_t.loc[(g_t_set_t['u'] - g_t_set_t['h'] >=0) & (g_t_set_t['u'] + g_t_set_t['h'] <=1)]

g_t_set = g_t_set.assign(lambda_ = lambda_function(g_t_set.h))
g_t_set = g_t_set.reset_index(drop=True)

a = psihat_statistic(y_data, g_t_set, sigmahat, epanechnikov_kernel)

psistar_statistic_list = []

'''for i in range(1000):
    z = np.random.normal(0, 1, T)
    z_temp = sigma * z
    psistar_statistic = psihat_statistic(z_temp, g_t_set, sigma, epanechnikov_kernel)
    psistar_statistic_list.append(psistar_statistic)
'''