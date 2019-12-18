import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import os
from scipy.io import loadmat
import scipy.linalg as la

def piecewise_est(x, y1, y2, order):

    Tbar = 2

    xLo = x[x < Tbar]
    xHi = x[x >= Tbar]
    y1Hi = y1[x >= Tbar]

    X = np.array([x, x**2]).T
    Y = np.log(y2)

    beta = np.linalg.inv(X.T@X)@X.T@Y
    g1 = beta[0]
    g2 = beta[1]

    X = np.array([(xHi - Tbar)**order]).T
    Y2 = np.log(y1Hi) - g1*xHi - g2*xHi**2

    g3 = np.linalg.inv(X.T@X)@X.T@Y2
    #g4 = g3 * 1.5

    yhatLo = np.exp(g1*xLo + g2*xLo**2)
    yhat1Hi = np.exp(g1*xHi + g2*xHi**2 + g3*(xHi - Tbar)**order)
    #yhat1_15Hi = np.exp(g1*xHi + g2*xHi**2 + g4*(xHi - Tbar)**order)
    #yhat1_05Hi = np.exp(g1*xHi + g2*xHi**2 + .5*g3*(xHi - Tbar)**order)

    #yhat1_05 = np.append(yhatLo, yhat1_05Hi)
    yhat1 = np.append(yhatLo, yhat1Hi)
    #yhat1_15 = np.append(yhatLo, yhat1_15Hi)
    yhat2 = np.exp(g1*x + g2*x**2)

    coeffs = [g1, g2, g3]

    return yhat1, yhat2, Tbar, coeffs

def piecewise_est_double(x, y1, y2, order):

    Tbar = 2

    xLo = x[x < Tbar]
    xHi = x[x >= Tbar]
    y1Hi = y1[x >= Tbar]

    X = np.array([x, x**2]).T
    Y = np.log(y2)

    beta = np.linalg.inv(X.T@X)@X.T@Y
    g1 = beta[0]
    g2 = beta[1]

    X = np.array([(xHi - Tbar)**order]).T
    Y2 = np.log(y1Hi) - g1*xHi - g2*xHi**2

    g3 = np.linalg.inv(X.T@X)@X.T@Y2
    #g4 = g3 * 1.5

    yhatLo = np.exp(g1*xLo + g2*xLo**2)
    yhat1Hi = np.exp(g1*xHi + g2*xHi**2 + 2 * g3*(xHi - Tbar)**order)
    #yhat1_15Hi = np.exp(g1*xHi + g2*xHi**2 + g4*(xHi - Tbar)**order)
    #yhat1_05Hi = np.exp(g1*xHi + g2*xHi**2 + .5*g3*(xHi - Tbar)**order)

    #yhat1_05 = np.append(yhatLo, yhat1_05Hi)
    yhat1 = np.append(yhatLo, yhat1Hi)
    #yhat1_15 = np.append(yhatLo, yhat1_15Hi)
    yhat2 = np.exp(g1*x + g2*x**2)

    coeffs = [g1, g2, g3]

    return yhat1, yhat2, Tbar, coeffs

def piecewise_est_quad(x, y1, y2, order):

    Tbar = 2

    xLo = x[x < Tbar]
    xHi = x[x >= Tbar]
    y1Hi = y1[x >= Tbar]

    X = np.array([x, x**2]).T
    Y = np.log(y2)

    beta = np.linalg.inv(X.T@X)@X.T@Y
    g1 = beta[0]
    g2 = beta[1]

    X = np.array([(xHi - Tbar)**order]).T
    Y2 = np.log(y1Hi) - g1*xHi - g2*xHi**2

    g3 = np.linalg.inv(X.T@X)@X.T@Y2
    #g4 = g3 * 1.5

    yhatLo = np.exp(g1*xLo + g2*xLo**2)
    yhat1Hi = np.exp(g1*xHi + g2*xHi**2 + 4 * g3*(xHi - Tbar)**order)
    #yhat1_15Hi = np.exp(g1*xHi + g2*xHi**2 + g4*(xHi - Tbar)**order)
    #yhat1_05Hi = np.exp(g1*xHi + g2*xHi**2 + .5*g3*(xHi - Tbar)**order)

    #yhat1_05 = np.append(yhatLo, yhat1_05Hi)
    yhat1 = np.append(yhatLo, yhat1Hi)
    #yhat1_15 = np.append(yhatLo, yhat1_15Hi)
    yhat2 = np.exp(g1*x + g2*x**2)

    coeffs = [g1, g2, g3]

    return yhat1, yhat2, Tbar, coeffs

def gen_distributions(xi):

    data = np.loadtxt("./data/TCRE_MacDougallEtAl2017_update.csv", skiprows=1, delimiter=',')

    sigma = np.std(data, ddof = 1)
    mu = np.mean(data)

#     dom = np.arange(0, 4.01, .05)
#     macdougall = norm.pdf(dom, mu, sigma)

#     with open("{}/mean_distort_0yr.txt".format(xi), 'r') as f:
#         lines = f.readlines()
#     for line in lines:
#         words = line.split(" ")
#         if "mean-distortion-nordhaus" in words[0]:
#             mean_dist_n_0 = float(words[1])
#         elif "mean-distortion-weitzman" in words[0]:
#             mean_dist_w_0 = float(words[1])

#     with open("{}/mean_distort_100yr.txt".format(xi), 'r') as f:
#         lines = f.readlines()
#     for line in lines:
#         words = line.split(" ")
#         if "mean-distortion-nordhaus" in words[0]:
#             mean_dist_n_100 = float(words[1])
#         elif "mean-distortion-weitzman" in words[0]:
#             mean_dist_w_100 = float(words[1])

#     n_0 = norm.pdf(dom, mu + mean_dist_n_0 * 1000, sigma)
#     w_0 = norm.pdf(dom, mu + mean_dist_w_0 * 1000, sigma)
#     n_100 = norm.pdf(dom, mu + mean_dist_n_100 * 1000, sigma)
#     w_100 = norm.pdf(dom, mu + mean_dist_w_100 * 1000, sigma)

#     return dom, macdougall, n_0, w_0, n_100, w_100, sigma + mu, sigma, mu
    return sigma, mu

def Burke_bootstrap(x, n_sims):
    mu_1 = 1.272e-02
    mu_2 = -4.871e-04
    sigma_1 = 3.248e-03
    sigma_2 = 1.029e-04
    rho_12 = -2.859133e-07
    Tbar = 13

    R = np.random.multivariate_normal(np.array([mu_1, mu_2]), \
                                      np.array([[sigma_1**2,rho_12],
                                                [rho_12, sigma_2**2]]),\
                                      n_sims)

    t_bars = -R[:,0] / (2 * R[:,1])
    trans_dom = np.tile(x, (n_sims, 1)) + Tbar
    damg_func = R[:, 1, np.newaxis] * trans_dom**2 + R[:,0, np.newaxis] * trans_dom
    maxs = R[:,1] * t_bars**2 + R[:,0] * t_bars
#     maxs = np.amax(damg_func, axis = 1)
    damg_func = damg_func - maxs[:,np.newaxis]

    dec2 = np.percentile(damg_func, 20, axis = 0)
    dec4 = np.percentile(damg_func, 40, axis = 0)
    dec6 = np.percentile(damg_func, 60, axis = 0)
    dec8 = np.percentile(damg_func, 80, axis = 0)

    return dec2, dec4, dec6, dec8

def quad_int(f, a, b, n, method):
    #This function takes a function f to integrate from the multidimensional
    #interval specified by the row vectors a and b. N different points are used
    #in the quadrature method. Legendre and Hermite basis functions are
    #currently supported. In the case of Hermite methodology b is the normal
    #density and a is the normal mean.

    #Created by John Wilson (johnrwilson@uchicago.edu)

    if method == "legendre":
        xs, ws = np.polynomial.legendre.leggauss(n)
        g = lambda x: f((b - a) / 2 * x + (a+b)/2);
        s = np.prod((b - a) / 2);

    elif method == "hermite":
        xs, ws = np.polynomial.hermite.hermgauss(n)
        g = lambda x: f(np.sqrt(2) * b * x + a)
        s = 1 / np.sqrt(np.pi)

    else:
        raise ValueError("Invalid 'method' parameter used in quadrature.")

    sum = 0
    if hasattr(a, "__len__"):
        dim = len(a)
    else:
        dim = 1

    if dim == 3:
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    sum = sum + (ws[i] * ws[j] * ws[k]) * g([xs[i], xs[j], xs[k]])

    if dim == 2:
        for i in range(n):
            for j in range(n):
                sum = sum + (ws[i] * ws[j]) * g([xs[i], xs[j]])

    elif dim == 1:
        for i in range(n):
            sum = sum + ws[i] * g(xs[i])

    result = s * sum;
    return result

def get_emissions(xi):
    files = os.listdir(str(xi)+"/")
    file_name = [f for f in files if "emission" in f][0]
    data = loadmat("{}/{}".format(str(xi), file_name))
    data_key = [k for k in data.keys() if 'e_value' in k][0]
    data = data[data_key]
    return np.array(data)[:,0]

def get_SCC(xi):
    files = os.listdir(str(xi)+"/")
    file_name = [f for f in files if "SCC" in f][0]
    data = loadmat("{}/{}".format(str(xi), file_name))
    total_SCC = np.array(data['SCC'])[:,0]
    external_SCC = np.array(data['SCC2'])[:,0]
    uncertainty_SCC = np.array(data['SCC3'])[:,0]
    private_SCC = np.array(data['SCC1'])[:,0]
    return total_SCC, external_SCC, uncertainty_SCC, private_SCC

def get_low_dmg_SCC(xi):
    files = os.listdir(str(xi)+"/")
    file_name = [f for f in files if "Low_dmg_SCC" in f][0]
    data = loadmat("{}/{}".format(str(xi), file_name))
    total_SCC = np.array(data['SCC'])[:,0]
    external_SCC = np.array(data['SCC2'])[:,0]
    uncertainty_SCC = np.array(data['SCC3'])[:,0]
    private_SCC = np.array(data['SCC1'])[:,0]
    return total_SCC, external_SCC, uncertainty_SCC, private_SCC

if __name__ == "__main__":

    est_max = 5
    plot_max = 3.5
    order = 2
    x = np.arange(0, est_max + 0.01, 0.01)
    y_w = (1 / (1 + (x / 20.46) **2 + (x / 6.081) ** 6.754))
    y_n = (1 / (1 + 0.00227 * x ** 2))

    yhat_w, yhat_n, Tbar, coeffs = piecewise_est(x, y_w, y_n, order)
    # plt.plot(x, yhat_w)
    # plt.plot(x, yhat_n)
    # plt.xlim([0, plot_max])
    # plt.ylim([np.min(yhat_w[x <= plot_max]) * .98, 1.01])
    # plt.show()

    dom, macdougall, n_0, w_0, n_100, w_100, std1, _, _ = gen_distributions(0.0001)
    # plt.plot(dom, macdougall)
    # plt.plot(dom, n_0)
    # plt.vlines(std1, 0, 1)
    # plt.ylim([0,1])
    # plt.xlim([0,4])
    # plt.show()

    # plt.plot(dom, macdougall)
    # plt.plot(dom, n_100)
    # plt.plot(dom, w_100)
    # plt.ylim([0,1])
    # plt.xlim([0,4])
    # plt.show()

    x = np.arange(0, 3.51, .1)
    dec2, dec4, dec6, dec8 = Burke_bootstrap(x, 1000)
    # plt.plot(x, dec2)
    # plt.plot(x, dec4)
    # plt.plot(x, dec6)
    # plt.plot(x, dec8)
    # plt.show()

    e = get_emissions(0.0002)
    # print(e)

    e = get_SCC(0.0002)
    # print(e)

    h = lambda x: np.array([x[0]**2 * x[1], x[1], x[0] - x[1]])
    res = quad_int(h, np.array([1,0]), np.array([2,-1]), 3, 'legendre')
    print(res)
