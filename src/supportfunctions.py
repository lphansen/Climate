import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.stats import norm
import pdb
import warnings
import time
from scipy.interpolate import RegularGridInterpolator, CubicSpline
from plotly.subplots import make_subplots
import plotly.graph_objs as go
from plotly.offline import init_notebook_mode, iplot
import pickle
import SolveLinSys

def finiteDiff(data, dim, order, dlt, cap = None):  
    # compute the central difference derivatives for given input and dimensions
    res = np.zeros(data.shape)
    l = len(data.shape)
    if l == 3:
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
    elif l == 2:
        if order == 1:                    # first order derivatives
            
            if dim == 0:                  # to first dimension

                res[1:-1,:] = (1 / (2 * dlt)) * (data[2:,:] - data[:-2,:])
                res[-1,:] = (1 / dlt) * (data[-1,:] - data[-2,:])
                res[0,:] = (1 / dlt) * (data[1,:] - data[0,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1] = (1 / (2 * dlt)) * (data[:,2:] - data[:,:-2])
                res[:,-1] = (1 / dlt) * (data[:,-1] - data[:,-2])
                res[:,0] = (1 / dlt) * (data[:,1] - data[:,0])

            else:
                raise ValueError('wrong dim')
                
        elif order == 2:
            
            if dim == 0:                  # to first dimension

                res[1:-1,:] = (1 / dlt ** 2) * (data[2:,:] + data[:-2,:] - 2 * data[1:-1,:])
                res[-1,:] = (1 / dlt ** 2) * (data[-1,:] + data[-3,:] - 2 * data[-2,:])
                res[0,:] = (1 / dlt ** 2) * (data[2,:] + data[0,:] - 2 * data[1,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1] = (1 / dlt ** 2) * (data[:,2:] + data[:,:-2] - 2 * data[:,1:-1])
                res[:,-1] = (1 / dlt ** 2) * (data[:,-1] + data[:,-3] - 2 * data[:,-2])
                res[:,0] = (1 / dlt ** 2) * (data[:,2] + data[:,0] - 2 * data[:,1])

            else:
                raise ValueError('wrong dim')
            
        else:
            raise ValueError('wrong order')
    else:
        raise ValueError("Dimension NOT supported")
        
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


def cap(x, lb, ub):
    if x <= ub or x >= lb:
        return x
    else:
        if x > ub:
            return ub
        else:
            return lb


def PDESolver(stateSpace, A, B_r, B_f, B_k, C_rr, C_ff, C_kk, D, v0, ε = 1, tol = -10, smartguess = False, solverType = 'False Transient'):

    if solverType == 'False Transient':
        A = A.reshape(-1,1,order = 'F')
        B = np.hstack([B_r.reshape(-1,1,order = 'F'),B_f.reshape(-1,1,order = 'F'),B_k.reshape(-1,1,order = 'F')])
        C = np.hstack([C_rr.reshape(-1,1,order = 'F'), C_ff.reshape(-1,1,order = 'F'), C_kk.reshape(-1,1,order = 'F')])
        D = D.reshape(-1,1,order = 'F')
        v0 = v0.reshape(-1,1,order = 'F')
        out = SolveLinSys.solveFT(stateSpace, A, B, C, D, v0, ε, tol)

        return out

    elif solverType == 'Feyman Kac':
        
        if smartguess:
            iters = 1
        else:
            iters = 400000
            
        A = A.reshape(-1, 1, order='F')
        B = np.hstack([B_r.reshape(-1, 1, order='F'), B_f.reshape(-1, 1, order='F'), B_k.reshape(-1, 1, order='F')])
        C = np.hstack([C_rr.reshape(-1, 1, order='F'), C_ff.reshape(-1, 1, order='F'), C_kk.reshape(-1, 1, order='F')])
        D = D.reshape(-1, 1, order='F')
        v0 = v0.reshape(-1, 1, order='F')
        out = SolveLinSys.solveFK(stateSpace, A, B, C, D, v0, iters)

        return out

def densityPlot(beta_f_space, Dists, key = 'Weighted'):
    years = [50, 75, 100]

    titles = ["Year {}".format(year) for year in years]

    fig = make_subplots(1, len(years), print_grid = False, subplot_titles = titles)

    dom =beta_f_space
    inds = ((dom>=0) & (dom<=5e-3))

    for i, year in enumerate(years):
        # data = loadmat("{}/50-50 weight/Dist_{}yr.mat".format(quad_rule, year))
        data = Dists
        if key == 'Weighted': 
            if i == 0:
                fig.add_scatter(x = dom[inds] * 1000, y = data['Original'][inds], row = 1, col = i + 1,
                    name = 'Original Distribution', line = dict(color = '#1f77b4', width = 3), showlegend = True, legendgroup = 'Original Distribution')
                fig.add_scatter(x = dom[inds] * 1000, y = data['Nordhaus_year' + str(year)][inds], row = 1, col = i + 1,
                    name = 'Low Damage Function', line = dict(color = 'red', dash='dashdot', width = 3), showlegend = True, legendgroup = 'Low Damage Function')
                fig.add_scatter(x = dom[inds] * 1000, y = data['Weitzman_year' + str(year)][inds], row = 1, col = i + 1,
                    name = 'High Damage Function', line = dict(color = 'green', dash='dash', width = 3), showlegend = True, legendgroup = 'High Damage Function')
            else:
                fig.add_scatter(x = dom[inds] * 1000, y = data['Original'][inds], row = 1, col = i + 1,
                    name = 'Original Distribution', line = dict(color = '#1f77b4', width = 3), showlegend = False, legendgroup = 'Original Distribution')
                fig.add_scatter(x = dom[inds] * 1000, y = data['Nordhaus_year' + str(year)][inds], row = 1, col = i + 1,
                    name = 'Low Damage Function', line = dict(color = 'red', dash='dashdot', width = 3), showlegend = False, legendgroup = 'Low Damage Function')
                fig.add_scatter(x = dom[inds] * 1000, y = data['Weitzman_year' + str(year)][inds], row = 1, col = i + 1,
                    name = 'High Damage Function', line = dict(color = 'green', dash='dash', width = 3), showlegend = False, legendgroup = 'High Damage Function')

        elif key == 'High':
            if i == 0:
                fig.add_scatter(x = dom[inds] * 1000, y = data['Original'][inds], row = 1, col = i + 1,
                    name = 'Original Distribution', line = dict(color = '#1f77b4', width = 3), showlegend = True, legendgroup = 'Original Distribution')
                fig.add_scatter(x = dom[inds] * 1000, y = data['Weitzman_year' + str(year)][inds], row = 1, col = i + 1,
                    name = 'High Damage Function', line = dict(color = 'green', dash='dash', width = 3), showlegend = True, legendgroup = 'High Damage Function')
            else:
                fig.add_scatter(x = dom[inds] * 1000, y = data['Original'][inds], row = 1, col = i + 1,
                    name = 'Original Distribution', line = dict(color = '#1f77b4', width = 3), showlegend = False, legendgroup = 'Original Distribution')
                fig.add_scatter(x = dom[inds] * 1000, y = data['Weitzman_year' + str(year)][inds], row = 1, col = i + 1,
                    name = 'High Damage Function', line = dict(color = 'green', dash='dash', width = 3), showlegend = False, legendgroup = 'High Damage Function')


        elif key == 'Low':
            if i == 0:
                fig.add_scatter(x = dom[inds] * 1000, y = data['Original'][inds], row = 1, col = i + 1,
                    name = 'Original Distribution', line = dict(color = '#1f77b4', width = 3), showlegend = True, legendgroup = 'Original Distribution')
                fig.add_scatter(x = dom[inds] * 1000, y = data['Nordhaus_year' + str(year)][inds], row = 1, col = i + 1,
                    name = 'Low Damage Function', line = dict(color = 'red', dash='dashdot', width = 3), showlegend = True, legendgroup = 'Low Damage Function')
            else:
                fig.add_scatter(x = dom[inds] * 1000, y = data['Original'][inds], row = 1, col = i + 1,
                    name = 'Original Distribution', line = dict(color = '#1f77b4', width = 3), showlegend = False, legendgroup = 'Original Distribution')
                fig.add_scatter(x = dom[inds] * 1000, y = data['Nordhaus_year' + str(year)][inds], row = 1, col = i + 1,
                    name = 'Low Damage Function', line = dict(color = 'red', dash='dashdot', width = 3), showlegend = False, legendgroup = 'Low Damage Function')

    fig['layout'].update(title = key + " Damage Specification", showlegend = True, titlefont = dict(size = 20), height = 450)

    for i in range(len(years)):

        fig['layout']['yaxis{}'.format(i+1)].update(showgrid = False)
        fig['layout']['xaxis{}'.format(i+1)].update(showgrid = False)

    fig['layout']['yaxis1'].update(title=go.layout.yaxis.Title(
                                    text="Probability Density", font=dict(size=16)))
    fig['layout']['xaxis2'].update(title=go.layout.xaxis.Title(
                                    text="Climate Sensitivity", font=dict(size=16)), showgrid = False)

    fig = go.FigureWidget(fig)
    iplot(fig)

def SCCDecomposePlot(SCCs, key = 'Weighted'):

    if key == 'Low':

        data = SCCs
        x1, y1, x2, y2, x3, y3 = 60, 195, 93, 330, 96, 100

    elif key == 'Weighted':

        data = SCCs
        x1, y1, x2, y2, x3, y3 = 60, 320, 80, 315, 90, 350

    elif key == 'High':

        data = SCCs
        x1, y1, x2, y2, x3, y3 = 60, 340, 93, 495, 96, 430


    total_SCC = np.array(data['SCC'])
    external_SCC = np.array(data['SCC2'])
    uncertainty_SCC = np.array(data['SCC3'])
    private_SCC = np.array(data['SCC1'])
    x = np.linspace(0,100,400)

    total = go.Scatter(x = x, y = total_SCC,
                   name = 'Total', line = dict(color = '#1f77b4', dash = 'solid', width = 3),\
                       showlegend = False)
    external = go.Scatter(x = x, y = external_SCC,
                   name = 'Uncertainty', line = dict(color = 'red', dash = 'dot', width = 3),\
                          showlegend = False)
    uncertainty = go.Scatter(x = x, y = uncertainty_SCC,
                   name = 'External', line = dict(color = 'green', dash = 'dashdot', width = 3),\
                             showlegend = False)
    private = go.Scatter(x = x, y = private_SCC,
                   name = 'Private', line = dict(color = 'black', width = 3),\
                         showlegend = False)

    annotations=[dict(x=x1, y=y1, text="Total", textangle=0, ax=-100,
                ay=-75, font=dict(color="black", size=12), arrowcolor="black",
                arrowsize=3, arrowwidth=1, arrowhead=1),

                dict(x=x2, y=y2, text="Uncertainty", textangle=0, ax=-100,
                ay=0, font=dict(color="black", size=12), arrowcolor="black",
                arrowsize=3, arrowwidth=1, arrowhead=1)]

                # dict(x=x3, y=y3, text="External", textangle=0, ax=-80,
                # ay=80, font=dict(color="black", size=12), arrowcolor="black",
                # arrowsize=3, arrowwidth=1, arrowhead=1)]

    layout = dict(title = 'Social Cost of Carbon, {} Damage Specification'.format(key),
                  titlefont = dict(size = 20),
                  xaxis = go.layout.XAxis(title=go.layout.xaxis.Title(
                                    text='Years', font=dict(size=16)),
                                         tickfont=dict(size=12), showgrid = False),
                  yaxis = go.layout.YAxis(title=go.layout.yaxis.Title(
                                    text='Dollars per Ton of Carbon', font=dict(size=16)),
                                         tickfont=dict(size=12), showgrid = False), 
                  annotations=annotations
                  )

    fig = dict(data = [total,  uncertainty], layout = layout)
    iplot(fig)
    
def emissionPlot(damageSpec, ξ, e_hists):

	colors = {'High': 'red', 'Low': 'green', 'Weighted': '#1f77b4'}
	lines = {'Averse': 'solid', "Neutral": 'dashdot'}

	# damageSpecs = ['High', 'Low', 'Weighted']
	# aversionSpecs = ['Averse', 'Neutral']
	# colors = ['green', '#1f77b4', 'red']
	# lines = ['solid', 'dashdot'] 

	x = np.linspace(0, 100, 400)
	data = []

	data.append(go.Scatter(x = x, y = e_hists[:,0], name = damageSpec +  ' Damage w/ ξ= {}'.format(ξ),
	    line = dict(width = 2, dash = 'solid', color = colors[damageSpec]), showlegend = True))

	layout = dict(title = 'Emissions Plot with {} Damage Setting, ξ = {}'.format(damageSpec, ξ),
	  titlefont = dict(size = 20),
	  xaxis = go.layout.XAxis(title=go.layout.xaxis.Title(
	                    text='Years', font=dict(size=16)),
	                         tickfont=dict(size=12), showgrid = False, showline = True),
	  yaxis = go.layout.YAxis(title=go.layout.yaxis.Title(
	                    text='Gigatons of Carbon', font=dict(size=16)),
	                         tickfont=dict(size=12), showgrid = False),
	  legend = dict(orientation = 'h', y = 1.15)
	  )

	fig = dict(data = data, layout = layout)
	iplot(fig)

def growthdensityPlot(beta_f_space, Dists):
    years = [50, 75, 100]

    titles = ["Year {}".format(year) for year in years]
    tilt_colors = ['#FFFF00', '#FFE600', '#FFCC00', '#FFB300', '#FF9900', '#FF8000', '#FF6600', '#FF4D00', '#FF3300', '#FF1A00', '#FF0000']

    fig = make_subplots(1, len(years), print_grid = False, subplot_titles = titles)

    dom = beta_f_space
    inds = ((dom>=0) & (dom<=5e-3))

    for i, year in enumerate(years):
        data = Dists
        if i == 0:
            fig.add_scatter(x = dom[inds] * 1000, y = data['Original'][inds], row = 1, col = i + 1,
                name = 'Original Distribution', line = dict(color = '#1f77b4', width = 3), showlegend = True, legendgroup = 'Original Distribution')
            for j, tilt in enumerate(Dists['Year{}'.format(year)]['tilt_dist']):
                fig.add_scatter(x = dom[inds] * 1000, y = tilt[inds], row = 1, col = i + 1,
                    name = 'Tilted {}'.format(j+1), line = dict(color = tilt_colors[j], dash='dash', width = 2), showlegend = True, legendgroup = 'Tilted Densities {}'.format(j))
        else:
            fig.add_scatter(x = dom[inds] * 1000, y = data['Original'][inds], row = 1, col = i + 1,
                name = 'Original Distribution', line = dict(color = '#1f77b4', width = 3), showlegend = False, legendgroup = 'Original Distribution')
            for j, tilt in enumerate(Dists['Year{}'.format(year)]['tilt_dist']):
                fig.add_scatter(x = dom[inds] * 1000, y = tilt[inds], row = 1, col = i + 1,
                    name = 'Tilted {}'.format(j+1), line = dict(color = tilt_colors[j], dash='dash', width = 2), showlegend = False, legendgroup = 'Tilted Densities {}'.format(j))


    fig['layout'].update(title = "Worst Case Probabilities, Growth Damage Specification", showlegend = True, titlefont = dict(size = 20), height = 450)

    for i in range(len(years)):

        fig['layout']['yaxis{}'.format(i+1)].update(showgrid = False)
        fig['layout']['xaxis{}'.format(i+1)].update(showgrid = False)

    fig['layout']['yaxis1'].update(title=go.layout.yaxis.Title(
                                    text="Probability Density", font=dict(size=16)))
    fig['layout']['xaxis2'].update(title=go.layout.xaxis.Title(
                                    text="Climate Sensitivity", font=dict(size=16)), showgrid = False)

    fig = go.FigureWidget(fig)
    iplot(fig)

def growthemissionPlot(ξ, e_hists):

    colors = {'High': 'red', 'Low': 'green', 'Weighted': '#1f77b4'}
    lines = {'Averse': 'solid', "Neutral": 'dashdot'}

    # damageSpecs = ['High', 'Low', 'Weighted']
    # aversionSpecs = ['Averse', 'Neutral']
    # colors = ['green', '#1f77b4', 'red']
    # lines = ['solid', 'dashdot'] 
    if ξ < 1:  # Averse
        key = 'Growth Averse'

    else:      # Neutral
        key = 'Growth Neutral'

    x = np.linspace(0, 100, 400)
    data = []

    data.append(go.Scatter(x = x, y = e_hists[:,0], name = key +  ' Damage w/ ξ= {:4f}'.format(ξ),
        line = dict(width = 2, dash = 'solid', color = '#1f77b4'), showlegend = True))

    layout = dict(title = 'Emissions Plot with {} Setting, ξ = {:4f}'.format(key, ξ),
      titlefont = dict(size = 20),
      xaxis = go.layout.XAxis(title=go.layout.xaxis.Title(
                        text='Years', font=dict(size=16)),
                             tickfont=dict(size=12), showgrid = False, showline = True),
      yaxis = go.layout.YAxis(title=go.layout.yaxis.Title(
                        text='Gigatons of Carbon', font=dict(size=16)),
                             tickfont=dict(size=12), showgrid = False),
      legend = dict(orientation = 'h', y = 1.15)
      )

    fig = dict(data = data, layout = layout)
    iplot(fig)

def growthSCCDecomposePlot(SCCs, ξ):

    if ξ < 1:  # Averse
        data = SCCs
        x1, y1, x2, y2 = 60, 1500, 80, 1120
        lgd = False
        key = 'Growth Averse'
        
        annotations=[dict(x=x1, y=y1, text="Total", textangle=0, ax=-100,
            ay=-75, font=dict(color="black", size=12), arrowcolor="black",
            arrowsize=3, arrowwidth=1, arrowhead=1),

            dict(x=x2, y=y2, text="Uncertainty", textangle=0, ax=-100,
            ay=0, font=dict(color="black", size=12), arrowcolor="black",
            arrowsize=3, arrowwidth=1, arrowhead=1)]

    else:      # Neutral
        data = SCCs
        x1, y1, x2, y2 = 60, 900, 80, 0
        lgd = True
        key = 'Growth Neutral'
        annotations=[dict(x=x1, y=y1, text="Total", textangle=0, ax=-100,
            ay=-75, font=dict(color="black", size=12), arrowcolor="black",
            arrowsize=3, arrowwidth=1, arrowhead=1),

            dict(x=x2, y=y2, text="Uncertainty", textangle=0, ax=-100,
            ay=-50, font=dict(color="black", size=12), arrowcolor="black",
            arrowsize=3, arrowwidth=1, arrowhead=1)]


    total_SCC = np.array(data['SCC'])
    external_SCC = np.array(data['SCC2'])
    uncertainty_SCC = np.array(data['SCC3'])
    private_SCC = np.array(data['SCC1'])
    x = np.linspace(0,100,400)

    total = go.Scatter(x = x, y = total_SCC,
           name = 'Total', line = dict(color = '#1f77b4', dash = 'solid', width = 3),\
               showlegend = lgd)
    external = go.Scatter(x = x, y = external_SCC,
           name = 'Ambiguity', line = dict(color = 'green', dash = 'dashdot', width = 3),\
                  showlegend = lgd)
    uncertainty = go.Scatter(x = x, y = uncertainty_SCC,
           name = 'Uncertainty', line = dict(color = 'red', dash = 'dot', width = 3),\
                     showlegend = lgd)
    private = go.Scatter(x = x, y = private_SCC,
           name = 'Private', line = dict(color = 'black', width = 3),\
                 showlegend = False)

    

        # dict(x=x3, y=y3, text="External", textangle=0, ax=-80,
        #     ay=80, font=dict(color="black", size=12), arrowcolor="black",
        #     arrowsize=3, arrowwidth=1, arrowhead=1)]

    layout = dict(title = 'Social Cost of Carbon, {} Specification'.format(key), 
              titlefont = dict(size = 20),
              xaxis = go.layout.XAxis(title=go.layout.xaxis.Title(
                                text='Years', font=dict(size=16)),
                                     tickfont=dict(size=12), showgrid = False),
              yaxis = go.layout.YAxis(title=go.layout.yaxis.Title(
                                text='Dollars per Ton of Carbon', font=dict(size=16)),
                                     tickfont=dict(size=12), showgrid = False), 
              annotations=annotations
              )

    fig = dict(data = [total, uncertainty], layout = layout)
    iplot(fig)