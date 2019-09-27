import numpy as np

def finiteDiff(data, dim, order, dlt, cap = None):  
# compute the central difference derivatives for given input and dimensions
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


def cap(x, lb, ub):
    if x <= ub or x >= lb:
        return x
    else:
        if x > ub:
            return ub
        else:
            return lb