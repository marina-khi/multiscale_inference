import numpy as np
from math import sqrt, log

def epanechnikov_kernel(x):
    '''
    Epanechnikov kernel function,
    which is defined as f(x) = 3/4(1-x^2)
    for |x|<=1 and 0 elsewhere
    ''' 
    epanechnikov = 3/4 * (1 - x*x) if abs(x)<=1 else 0
    return epanechnikov

def lambda_function(h_list):
    '''
    Additive correction tern \lambda(h)
    that depends only on the bandwidth h
    '''
    for h in h_list:
        try:
            lambda_h = sqrt(2*log(1/(2*h)))
        except ValueError:
            print("h is exceeding h_max")
        return lambda_h

def omega_function(T, u, h, k_function):
    '''
    Kernel weight function that calculates
    the vector of (\omega_{1,T}(u,h), \ldots, \omega_{T,T}(u,h)).
    Meanwhile, the function also calculates
    the norm ||K||_{u,h,T} of the kernel function.
    '''
    result_temp = np.zeros(T)
    K_norm_temp = 0
    for i in range(T):
        x = (u - (i+1)/T)/h
        k = k_function(x)
        result_temp[i] = k
        K_norm_temp = K_norm_temp + k*k
    K_norm = sqrt(K_norm_temp)
    result = result_temp / K_norm
    return result



def psi_average(data, u_list, h_list, sigmahat, k_function):
    '''
    Kernel average function \psi_T(u,h) that takes u, h, data
    and the type of kernel function as arguments.
    The data can be y_data for \hat{\Psi}_T or independent gaussian rv
    z_temp = sigma*z for \Psi^star_T.
    The output is one value for each u and h.
    '''
    T = len(data)
    '''psi_average_normed_list = []
    for u, h in zip(u_list, h_list):
        psi_average_normed = sum(omega_function(T, u, h, k_function)*data)
        psi_average_normed_list.append(abs(psi_average_normed)/sigmahat)
    '''
    psi_average_normed_list = [abs(sum(omega_function(T, u, h, k_function) * data))/sigmahat for u, h in zip(u_list, h_list)]
    return psi_average_normed_list


def psihat_statistic(y_data, g_t_set, sigmahat, k_function = epanechnikov_kernel):
    '''
    Function that calculates the multiscale statistic \hat{\Psi}_T.
    It takes the following entities as arguments.
        y-data: the data
        g_t_set: range of different locations u and bandwidths h
        k_function: type of kernel function
        sigmahat: the estimator of the square root of the long-run error variance \sigma^2
    It produces the value of the test statistic as an output
    '''
    g_t_set.loc[:,'values'] = psi_average(y_data, g_t_set['u'], g_t_set['h'], sigmahat, k_function) - g_t_set['lambda_']
    result = g_t_set['values'].max()
    return result
'''
# Function that calculates the auxiliary statistic \Psi^star_T.
# It takes the following entities as arguments.
#   z_t: the independent standard normal random variables
#   g_t_set: range of different locations u and bandwidths h
#   k_function: type of kernel function
#   sigma: the square root of the long-run error variance \sigma^2
# It produces the value of the test statistic as an output
# which is used further to calculate the critical values

psistar_statistic <- function(z, g_t_set, k_function = epanechnikov_kernel, sigma) {
  g_t_set_card = nrow(g_t_set)
  z_temp = sigma*z
  for (i in 1:g_t_set_card) {
    g_t_set[['values']][i] <- abs(psi_average(z_temp, g_t_set[['u']][i], g_t_set[['h']][i], k_function)/sigma) - lambda(g_t_set[['h']][i])
  }
  result = max(g_t_set$values)
#  cat("Statistic with a star:", result)
  return(result)
}'''