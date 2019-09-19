# Required packages
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.stats import norm
import pdb
import warnings
import SolveLinSys
import time
import sys
sys.stdout.flush()

def centra_diff(data, dim, order, dlt, cap = None):  # compute the central difference derivatives for given input and dimensions
    res = np.zeros(data.shape)
    if order == 1:                    # first order derivatives
        
        if dim == 0:                  # to first dimension

            res[1:-1,:,:] = (1 / (2 * dlt)) * (data[2:,:,:] - data[:-2,:,:])
            res[-1,:,:] = (1 / dlt) * (data[-1,:,:] - data[-2,:,:])
            res[0,:,:] = (1 / dlt) * (data[1,:,:] - data[0,:,:])

        elif dim == 1:                # to second dimension

            res[:,1:-1,:] = (1 / (2 * dlt)) * (data[:,2:,:] - data[:,:-2,:])
            res[:,-1,:] = (1 / dlt) * (data[:,-1,:] - data[:,-2,:])
            res[:,0,:] = (1 / dlt) * (data[:,1,:] - data[:,0,:])

        elif dim == 2:                # to third dimension

            res[:,:,1:-1] = (1 / (2 * dlt)) * (data[:,:,2:] - data[:,:,:-2])
            res[:,:,-1] = (1 / dlt) * (data[:,:,-1] - data[:,:,-2])
            res[:,:,0] = (1 / dlt) * (data[:,:,1] - data[:,:,0])

        else:
            raise ValueError('wrong dim')
            
    elif order == 2:
        
        if dim == 0:                  # to first dimension

            res[1:-1,:,:] = (1 / dlt ** 2) * (data[2:,:,:] + data[:-2,:,:] - 2 * data[1:-1,:,:])
            res[-1,:,:] = (1 / dlt ** 2) * (data[-1,:,:] + data[-3,:,:] - 2 * data[-2,:,:])
            res[0,:,:] = (1 / dlt ** 2) * (data[2,:,:] + data[0,:,:] - 2 * data[1,:,:])

        elif dim == 1:                # to second dimension

            res[:,1:-1,:] = (1 / dlt ** 2) * (data[:,2:,:] + data[:,:-2,:] - 2 * data[:,1:-1,:])
            res[:,-1,:] = (1 / dlt ** 2) * (data[:,-1,:] + data[:,-3,:] - 2 * data[:,-2,:])
            res[:,0,:] = (1 / dlt ** 2) * (data[:,2,:] + data[:,0,:] - 2 * data[:,1,:])

        elif dim == 2:                # to third dimension

            res[:,:,1:-1] = (1 / dlt ** 2) * (data[:,:,2:] + data[:,:,:-2] - 2 * data[:,:,1:-1])
            res[:,:,-1] = (1 / dlt ** 2) * (data[:,:,-1] + data[:,:,-3] - 2 * data[:,:,-2])
            res[:,:,0] = (1 / dlt ** 2) * (data[:,:,2] + data[:,:,0] - 2 * data[:,:,1])

        else:
            raise ValueError('wrong dim')
        
    else:
        raise ValueError('wrong order')
        
    if cap is not None:
        res[res < cap] = cap
    return res

def quad_points_legendre(n):
    u = np.sqrt(1 / (4 - 1 / np.linspace(1,n-1,n-1)**2))  # upper diag
    [lambda0,V] = np.linalg.eig(np.diagflat(u,1) + np.diagflat(u,-1))  # V's column vectors are the main d
    i = np.argsort(lambda0)
    Vtop = V[0,:]
    Vtop = Vtop[i]
    w = 2 * Vtop ** 2
    return (lambda0[i],w)

def quad_points_hermite(n):
    i = np.linspace(1,n-1,n-1)
    a = np.sqrt(i / 2.0)
    [lambda0,V] = np.linalg.eig(np.diagflat(a,1) + np.diagflat(a,-1))
    i = np.argsort(lambda0)
    Vtop = V[0,:]
    Vtop = Vtop[i]
    w = np.sqrt(np.pi) * Vtop ** 2
    return (lambda0[i],w)


def quad_int(f,a,b,n,method):
    """
This function takes a function f to integrate from the multidimensional
interval specified by the row vectors a and b. N different points are used
in the quadrature method. Legendre and Hermite basis functions are
currently supported. In the case of Hermite methodology b is the normal
density and a is the normal mean.

Created by John Wilson (johnrwilson@uchicago.edu) & Updaed by Jiaming Wang (Jiamingwang@uchicago.edu)
    """
    if method == 'legendre':
        
        (xs,ws) = quad_points_legendre(n)
        g = lambda x: f((b-a) * 0.5  * x + (a + b) * 0.5)
        s = np.prod((b-a) * 0.5)                ######## Why using prod here?
        
    elif method == 'hermite':
        
        (xs,ws) = quad_points_hermite(n)
        g = lambda x: f(np.sqrt(2) * b * x + a)
        s = 1 / np.sqrt(np.pi)
        
    else:
        raise TypeError('Wrong polynomials specification')
    
    
    tp = type(a)
    if tp is np.float64 or tp is int or tp is np.double:
        res = 0
        for i in range(n):
#             pdb.set_trace()
            res += ws[i] * g(xs[i])
    else:
        raise ValueError('dimension is not 1')
    
    return s * res

# Data Cleansing
x = loadmat('climate_pars.mat')

for key,val in x.items():
    if key[0:2] == '__':
        pass
    elif key == 'filename' or key == 's' or key == 's0' or key =='s1' or key == 's2':
        pass
    else:
        (height, width) = val.shape
        if width == 1:
            if height == 1:    # one number
                exec(key +'=val[0,0]')
            else:              # a column array
                exec(key +'=val[:,0]')
        else:
            if height == 1:    # a row array
                exec(key +'=val[0,:]')
            else:
                exec(key +'=val')

start_time = time.time()

quadrature = 'legendre'   ###
theta = 5000
weight = 0.5

sigma_o = sigma_k
indic = 0    # redundant
beta_f = np.mean(par_lambda_McD)
var_beta_f = np.var(par_lambda_McD)
# var_beta_f = for_check['var_beta_f'][0][0]               # Worth discussion further
par_beta_f = par_lambda_McD
lambda0 = 1.0 / var_beta_f   # noted as lambda in MATLAB

xi_d = -1 * (1-alpha)
r_min = 0
r_max = 9
t_min = 0
t_max = 4000
k_min = 0
k_max = 9

tol = 1e-4

nr = 30
nt = 30
nk = 25
r = np.linspace(r_min, r_max, nr)
t = np.linspace(t_min, t_max, nt)
k = np.linspace(k_min, k_max, nk)

hr = r[1] - r[0]
hk = k[1] - k[0]
ht = t[1] - t[0]

(r_mat, t_mat, k_mat) = np.meshgrid(r,t,k, indexing = 'ij')
# r_mat = for_check['r_mat']
# t_mat = for_check['t_mat']
# k_mat = for_check['k_mat']

# time step
dt = 0.5
n = 30

# Power selection
power = 2


if power == 2:
    # Quad, max temp = 5:
    gamma_1 =  0.00017675
    gamma_2 =  2 * 0.0022
    gamma_2_plus = 2 * 0.0197
    sigma_1 = 0
    sigma_2 = 0
    rho_12 = 0
elif power == 3:
    # Cubic, max temp = 5:
    gamma_1 =  0.00017675
    gamma_2 =  2 * 0.0022
    gamma_2_plus = 3 * 0.0079 * 1.5
    sigma_1 = 0
    sigma_2 = 0
    rho_12 = 0
    
f_bar = 2
crit = 2

F0 = 1
sigma = np.matrix([[sigma_1 ** 2, rho_12], [rho_12, sigma_2 ** 2]])
Sigma = np.matrix([[var_beta_f,0,0],[0,sigma_1 ** 2, rho_12],[0, rho_12, sigma_2 ** 2]])
dee = np.matrix([gamma_1 + gamma_2 * F0 + gamma_2_plus * (F0 - f_bar) ** 2 * (F0>=2), beta_f, beta_f * F0])
sigma_d = float(np.sqrt(dee * Sigma * dee.T))

weight = 0.5  # on nordhaus
bar_gamma_2_plus = weight * 0 + (1 - weight) * gamma_2_plus

a = beta_f - 5 * np.sqrt(var_beta_f)
b = beta_f + 5 * np.sqrt(var_beta_f)

step_val1 = round(len(r) / 5)
step_val2 = round(len(t) / 5)
step_val3 = round(len(k) / 5)

v0 = alpha * r_mat + (1-alpha) * k_mat - beta_f * t_mat
v1_initial = v0 * np.ones(r_mat.shape)

v0_dr = centra_diff(v0,0,1,hr,1e-8)
v0_dt = centra_diff(v0,1,1,ht)
v0_dk = centra_diff(v0,2,1,hk)

v0_drr = centra_diff(v0,0,2,hr)
v0_dtt = centra_diff(v0,1,2,ht)
v0_dkk = centra_diff(v0,2,2,hk)

### FOC_all combinations
B1 = v0_dr - xi_d * (gamma_1 + gamma_2 * t_mat * beta_f + gamma_2_plus * (t_mat * beta_f - f_bar) ** (power - 1) * (t_mat >= (crit / beta_f))) * beta_f * np.exp(r_mat) - v0_dt * np.exp(r_mat)
C1 = - delta * alpha
e = -C1 / B1
e_hat = e
Acoeff = np.exp(r_mat - k_mat)
Bcoeff = delta * (1-alpha) / (np.exp(-r_mat + k_mat) * v0_dr * Gamma_r * 0.5) + v0_dk * Gamma / (np.exp(-r_mat + k_mat) * v0_dr * Gamma_r * 0.5)
Ccoeff = -A_O  - Theta
f = ((-Bcoeff + np.sqrt(Bcoeff ** 2 - 4 * Acoeff * Ccoeff)) / (2 * Acoeff)) ** 2
i_k = (v0_dk * Gamma / (np.exp(-r_mat + k_mat) * v0_dr * Gamma_r * 0.5)) * (f ** 0.5) - Theta

# lower damage model:
a_1 = np.zeros(r_mat.shape);
b_1 = xi_d * e_hat * np.exp(r_mat) * gamma_1
c_1 = 2 * xi_d * e_hat * np.exp(r_mat) * t_mat * gamma_2
lambda_tilde_1 = lambda0 + c_1 * theta
beta_tilde_1 = beta_f - c_1 * theta / lambda_tilde_1 * beta_f - theta / lambda_tilde_1 * b_1
I_1 = a_1 - 0.5 * np.log(lambda0) / theta + 0.5 * np.log(lambda_tilde_1) / theta + 0.5 * lambda0 * beta_f ** 2 / theta - 0.5 * lambda_tilde_1 * (beta_tilde_1) ** 2 / theta
#     R_1 = theta.*(I_1-(a_1+b_1.*beta_tilde_1+c_1./2.*(beta_tilde_1).^2+c_1./2./lambda_tilde_1));
R_1 = theta * (I_1 - (a_1 + b_1 * beta_tilde_1 + c_1 * beta_tilde_1 ** 2 + c_1 / lambda_tilde_1))
J_1_without_e = xi_d * (gamma_1 * beta_tilde_1 + gamma_2 * t_mat * (beta_tilde_1 ** 2 + 1 / lambda_tilde_1)) * np.exp(r_mat)

# J_1_without_e = xi_d.*(gamma_1.*beta_tilde_1 ...
#         +gamma_2.*t_mat.*(beta_tilde_1.^2+1./lambda_tilde_1)).*exp(r_mat) ;

pi_tilde_1 = weight * np.exp(-theta * I_1)

r_multi = np.exp(r)
r_multi = r_multi[:,np.newaxis,np.newaxis]

# high damage model
if quadrature == 'legendre':
    
    scale_2_fnc = lambda x: np.exp(-theta * xi_d * (gamma_1 * x + gamma_2 * x ** 2 * t_mat + gamma_2_plus * x * (x * t_mat - f_bar) ** (power - 1) * ((x * t_mat - f_bar) >= 0)) * r_multi * e_hat)  * norm.pdf(x,beta_f,np.sqrt(var_beta_f))
    scale_2 = quad_int(scale_2_fnc, a, b, n, 'legendre')   
    
elif quadrature == 'hermite':
    
    scale_2_fnc = lambda x: np.exp(-theta * xi_d * (gamma_1 * x + gamma_2 * x ** 2 * t_mat + gamma_2_plus * x * (x * t_mat - f_bar) ** (power - 1) * ((x * t_mat - f_bar) >= 0)) * r_multi * e_hat)  # * norm.pdf(x,beta_f,np.sqrt(var_beta_f))
    scale_2 = quad_int(scale_2_fnc, beta_f, np.sqrt(var_beta_f), n, 'hermite')
    
else: 
    
    raise TypeError('Wrong polynomials specification')

q2_tilde_fnc = lambda x: np.exp(-theta * xi_d * (gamma_1 * x + gamma_2 * x ** 2 * t_mat + gamma_2_plus * x * (x * t_mat - f_bar) ** (power - 1) * ((x * t_mat - f_bar) >= 0)) * r_multi * e_hat) / scale_2
I_2 = -1 / theta * np.log(scale_2)

#     q2_tilde_fnc = @(x) exp(-theta.*xi_d.*(gamma_1.*x ...
#         +gamma_2.*x.^2.*t_mat ...
#         +gamma_2_plus.*x.*(x.*t_mat-f_bar).^(power-1).*((x.*t_mat-f_bar)>=0)).*exp(r).*e_hat)./scale_2;

if quadrature == 'legendre':
    
    J_2_without_e_fnc = lambda x: xi_d * np.exp(r_mat) * q2_tilde_fnc(x) * (gamma_1 * x + gamma_2 * t_mat * x ** 2 + gamma_2_plus * x * (x * t_mat - f_bar) ** (power - 1) * ((x * t_mat - f_bar) >= 0)) * norm.pdf(x,beta_f,np.sqrt(var_beta_f))
    J_2_without_e = quad_int(J_2_without_e_fnc, a, b, n, 'legendre')
    
elif quadrature == 'hermite':
    
    J_2_without_e_fnc = lambda x: xi_d * np.exp(r_mat) * q2_tilde_fnc(x) * (gamma_1 * x + gamma_2 * t_mat * x ** 2 + gamma_2_plus * x * (x * t_mat - f_bar) ** (power - 1) * ((x * t_mat - f_bar) >= 0)) # * norm.pdf(x,beta_f,np.sqrt(var_beta_f))
    J_2_without_e = quad_int(J_2_without_e_fnc, beta_f, np.sqrt(var_beta_f), n, 'hermite')
    
J_2_with_e = J_2_without_e * e_hat

R_2 = theta * (I_2 - J_2_with_e)
pi_tilde_2 = (1 - weight) * np.exp(-theta * I_2)
pi_tilde_1_norm = pi_tilde_1 / (pi_tilde_1 + pi_tilde_2)
pi_tilde_2_norm = 1 - pi_tilde_1_norm

expec_e_sum = (pi_tilde_1_norm * J_1_without_e + pi_tilde_2_norm * J_2_without_e)

B1 = v0_dr - v0_dt * np.exp(r_mat) - expec_e_sum
C1 = -delta * alpha
e = -C1 / B1
e_star = e

J_1 = J_1_without_e * e_star
J_2 = J_2_without_e * e_star

I_term = -1 / theta * np.log(pi_tilde_1 + pi_tilde_2)

R_1 = theta * (I_1 - J_1)
R_2 = theta * (I_2 - J_2)
drift_distort = (pi_tilde_1_norm * J_1 + pi_tilde_2_norm * J_2)

RE = pi_tilde_1_norm * R_1 + pi_tilde_2_norm * R_2 + pi_tilde_1_norm * np.log(pi_tilde_1_norm / weight) + pi_tilde_2_norm * np.log(pi_tilde_2_norm / (1 - weight))
RE_total = 1 / theta * RE

### PDE Inputs

A = -delta * np.ones(r_mat.shape)
B_r = -e_star + Gamma_r * (f ** Theta_r) - 0.5 * (sigma_r ** 2)
B_k = Alpha + Gamma * np.log(1 + i_k / Theta) - 0.5 * (sigma_k ** 2)
B_t = e_star * np.exp(r_mat)
C_rr = 0.5 * sigma_r ** 2 * np.ones(r_mat.shape)
C_kk = 0.5 * sigma_k ** 2 * np.ones(r_mat.shape)
C_tt = np.zeros(r_mat.shape)

D = delta * alpha * np.log(e_star) + delta * alpha * r_mat + delta * (1 - alpha) * (np.log(A_O - i_k - f * np.exp(r_mat - k_mat)) + k_mat) + I_term #  + drift_distort + RE_total

stateSpace = np.hstack([r_mat.reshape(-1,1,order = 'F'),t_mat.reshape(-1,1,order = 'F'),k_mat.reshape(-1,1,order = 'F')])
A = A.reshape(-1,1,order = 'F')
B = np.hstack([B_r.reshape(-1,1,order = 'F'),B_t.reshape(-1,1,order = 'F'),B_k.reshape(-1,1,order = 'F')])
C = np.hstack([C_rr.reshape(-1,1,order = 'F'), C_tt.reshape(-1,1,order = 'F'), C_kk.reshape(-1,1,order = 'F')])
D = D.reshape(-1,1,order = 'F')
v01 = v0.reshape(-1,1,order = 'F')
dt = dt

out = SolveLinSys.solvels(stateSpace, A,B,C,D,v01,dt)

out_comp = out[2].reshape(v0.shape,order = "F")
episode = 0

print("Episode %d: PDE Error: %f, Conjugate Gradient Error %f after %d iterations" % (episode, np.max((out_comp - v0) / dt), out[1], out[0]))
v0 = out_comp
vold = v1_initial

while np.max((out_comp - vold) / dt) > tol:
# for iter in range(1):
    vold = v0
    v0_dr = centra_diff(v0,0,1,hr,1e-8)
    v0_dt = centra_diff(v0,1,1,ht)
    v0_dk = centra_diff(v0,2,1,hk)

    v0_drr = centra_diff(v0,0,2,hr)
    v0_dtt = centra_diff(v0,1,2,ht)
    v0_dkk = centra_diff(v0,2,2,hk)

    e_hat = e_star 

    f = ((A_O + Theta) * np.exp(-r_mat + k_mat) * (v0_dr * Gamma_r * Theta_r) / ((v0_dr * Gamma_r * Theta_r) * f ** (Theta_r) + (delta * (1-alpha) + v0_dk * Gamma))) ** (1 / (1 - Theta_r))
    f = f * (v0_dr > 1e-8)
    i_k = ((v0_dk * Gamma / (np.exp(-r_mat + k_mat) * v0_dr * Gamma_r * Theta_r)) * (f ** (1 - Theta_r)) - Theta) * (v0_dr > 1e-8) + (v0_dr <= 1e-8) * (v0_dk * Gamma * A_O - delta * (1-alpha) * Theta) / (delta * (1-alpha) + v0_dk * Gamma)

    # lower damage model:
    a_1 = np.zeros(r_mat.shape);
    b_1 = xi_d * e_hat * np.exp(r_mat) * gamma_1
    c_1 = 2 * xi_d * e_hat * np.exp(r_mat) * t_mat * gamma_2
    lambda_tilde_1 = lambda0 + c_1 * theta
    beta_tilde_1 = beta_f - c_1 * theta / lambda_tilde_1 * beta_f - theta / lambda_tilde_1 * b_1

    I_1 = a_1 - 0.5 * np.log(lambda0) / theta + 0.5 * np.log(lambda_tilde_1) / theta + 0.5 * lambda0 * beta_f ** 2 / theta - 0.5 * lambda_tilde_1 * (beta_tilde_1) ** 2 / theta

    R_1 = theta * (I_1 - (a_1 + b_1 * beta_tilde_1 + c_1 * beta_tilde_1 ** 2 + c_1 / lambda_tilde_1))
    J_1_without_e = xi_d * (gamma_1 * beta_tilde_1 + gamma_2 * t_mat * (beta_tilde_1 ** 2 + 1 / lambda_tilde_1)) * np.exp(r_mat)
    pi_tilde_1 = weight * np.exp(-theta * I_1)  

    #     r_multi = np.exp(r)
    #     r_multi = r_multi[:,np.newaxis,np.newaxis]

    ## high damage model
    if quadrature == 'legendre':

        scale_2_fnc = lambda x: np.exp(-theta * xi_d * (gamma_1 * x + gamma_2 * x ** 2 * t_mat + gamma_2_plus * x * (x * t_mat - f_bar) ** (power - 1) * ((x * t_mat - f_bar) >= 0)) * r_multi * e_hat)  * norm.pdf(x,beta_f,np.sqrt(var_beta_f))
        scale_2 = quad_int(scale_2_fnc, a, b, n, 'legendre')   

    elif quadrature == 'hermite':

        scale_2_fnc = lambda x: np.exp(-theta * xi_d * (gamma_1 * x + gamma_2 * x ** 2 * t_mat + gamma_2_plus * x * (x * t_mat - f_bar) ** (power - 1) * ((x * t_mat - f_bar) >= 0)) * r_multi * e_hat)  * norm.pdf(x,beta_f,np.sqrt(var_beta_f))
        scale_2 = quad_int(scale_2_fnc, beta_f, np.sqrt(var_beta_f), n, 'hermite')

    else: 

        raise TypeError('Wrong polynomials specification')

    q2_tilde_fnc = lambda x: np.exp(-theta * xi_d * (gamma_1 * x + gamma_2 * x ** 2 * t_mat + gamma_2_plus * x * (x * t_mat - f_bar) ** (power - 1) * ((x * t_mat - f_bar) >= 0)) * r_multi * e_hat) / scale_2
    I_2 = -1 / theta * np.log(scale_2)

    if quadrature == 'legendre':

        J_2_without_e_fnc = lambda x: xi_d * np.exp(r_mat) * q2_tilde_fnc(x) * (gamma_1 * x + gamma_2 * t_mat * x ** 2 + gamma_2_plus * x * (x * t_mat - f_bar) ** (power - 1) * ((x * t_mat - f_bar) >= 0)) * norm.pdf(x,beta_f,np.sqrt(var_beta_f))
        J_2_without_e = quad_int(J_2_without_e_fnc, a, b, n, 'legendre')

    elif quadrature == 'hermite':

        J_2_without_e_fnc = lambda x: xi_d * np.exp(r_mat) * q2_tilde_fnc(x) * (gamma_1 * x + gamma_2 * t_mat * x ** 2 + gamma_2_plus * x * (x * t_mat - f_bar) ** (power - 1) * ((x * t_mat - f_bar) >= 0)) * norm.pdf(x,beta_f,np.sqrt(var_beta_f))
        J_2_without_e = quad_int(J_2_without_e_fnc, beta_f, np.sqrt(var_beta_f), n, 'hermite')

    J_2_with_e = J_2_without_e * e_hat
    R_2 = theta * (I_2 - J_2_with_e)
    pi_tilde_2 = (1 - weight) * np.exp(-theta * I_2)
    pi_tilde_1_norm = pi_tilde_1 / (pi_tilde_1 + pi_tilde_2)
    pi_tilde_2_norm = 1 - pi_tilde_1_norm

    expec_e_sum = (pi_tilde_1_norm * J_1_without_e + pi_tilde_2_norm * J_2_without_e)

    B1 = v0_dr - v0_dt * np.exp(r_mat) - expec_e_sum
    C1 = -delta * alpha
    e = -C1 / B1
    e_star = e

    J_1 = J_1_without_e * e_star
    J_2 = J_2_without_e * e_star

    I_term = -1 / theta * np.log(pi_tilde_1 + pi_tilde_2)

    R_1 = theta * (I_1 - J_1)
    R_2 = theta * (I_2 - J_2)
    drift_distort = (pi_tilde_1_norm * J_1 + pi_tilde_2_norm * J_2)

    RE = pi_tilde_1_norm * R_1 + pi_tilde_2_norm * R_2 + pi_tilde_1_norm * np.log(pi_tilde_1_norm / weight) + pi_tilde_2_norm * np.log(pi_tilde_2_norm / (1 - weight))
    RE_total = 1 / theta * RE 

    ### PDE Inputs

    A = -delta * np.ones(r_mat.shape)
    B_r = -e_star + Gamma_r * (f ** Theta_r) - 0.5 * (sigma_r ** 2)
    B_k = Alpha + Gamma * np.log(1 + i_k / Theta) - 0.5 * (sigma_k ** 2)
    B_t = e_star * np.exp(r_mat)
    C_rr = 0.5 * sigma_r ** 2 * np.ones(r_mat.shape)
    C_kk = 0.5 * sigma_k ** 2 * np.ones(r_mat.shape)
    C_tt = np.zeros(r_mat.shape)

    D = delta * alpha * np.log(e_star) + delta * alpha * r_mat + delta * (1 - alpha) * (np.log(A_O - i_k - f * np.exp(r_mat - k_mat)) + k_mat) +  I_term # + drift_distort + RE_total

    stateSpace = np.hstack([r_mat.reshape(-1,1,order = 'F'),t_mat.reshape(-1,1,order = 'F'),k_mat.reshape(-1,1,order = 'F')])
    A = A.reshape(-1,1,order = 'F')
    B = np.hstack([B_r.reshape(-1,1,order = 'F'),B_t.reshape(-1,1,order = 'F'),B_k.reshape(-1,1,order = 'F')])
    C = np.hstack([C_rr.reshape(-1,1,order = 'F'), C_tt.reshape(-1,1,order = 'F'), C_kk.reshape(-1,1,order = 'F')])
    D = D.reshape(-1,1,order = 'F')
    v01 = v0.reshape(-1,1,order = 'F')
    dt = dt

    out = SolveLinSys.solvels(stateSpace, A,B,C,D,v01,dt)
    out_comp = out[2].reshape(v0.shape,order = "F")
    episode += 1

    v0 = out_comp
    print("Episode %d: PDE Error: %f, Conjugate Gradient Error %f after %d iterations" % (episode, np.max((out_comp - vold) / dt), out[1], out[0]))

print("--- %s seconds ---" % (time.time() - start_time))