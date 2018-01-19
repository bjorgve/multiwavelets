import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.legendre as leg

def legendre_polynomial(order, x):
    coef = []
    if type(x) is float:

        if x<= 0 or x>1:
            return 0
        else:
            for n in range(0, order+1):
                if(n == order):
                    coef.append(1)
                else:
                    coef.append(0)
            return (2.*order +1.)**.5*leg.Legendre(coef, domain=[-1,1])(2*x-1)
    else:
        y = np.zeros(x.size)

        for i in range(x.size):

            if x[i] <= 0 or x[i] > 1:
                y[i] = 0
            else:
                for n in range(0, order+1):
                    if(n == order):
                        coef.append(1)
                    else:
                        coef.append(0)

                y[i]=(2.*order +1.)**.5*leg.Legendre(coef, domain=[-1,1])(2*x[i]-1)
                coef.clear()
    return y

def scaling(order,scale, translation,x):
    return 2**(scale/2)*legendre_polynomial(order,2**scale*x-translation)



def roots(order):
    coef=[]
    for n in range(0,order+1):
        if n == order:
            coef.append(1)
        else:
            coef.append(0)
    return .5*(leg.legroots(coef)+1)

#return a np.array
def weights(order):

    x = roots(order)
    w = np.empty(order)
    coef = []
    for n in range(0, order+1):
        if(n == order):
            coef.append(1)
        else:
            coef.append(0)

    derivative_coefficients = leg.legder(coef)
    del coef[-1]
    del coef[-1]
    coef.append(1)
    for i in range(0,order):
        w[i] = 1/(order*leg.Legendre(derivative_coefficients, domain=[-1,1])(2*x[i] -1)* \
            leg.Legendre(coef, domain = [-1,1])(2*x[i]-1))
            #possibly some bug here, get a factor 2 different using mathematica
    return w

def np_scaling(order,f):
    x = np.arange(0,1,.01)

    return .5*np.trapz(legendre_polynomial(order,x)*f(2**(-2)*x),dx=0.01)

def approximated_scaling_coef(order,scale,translation,f):
    x = np.arange(0,1,.01)
    return 2**(-.5*scale)*np.trapz(legendre_polynomial(order,x)*f(2**(-scale)*(x + translation)),dx=0.01)

def scaling_coef(order,order_max,scale,translation,f):

    tmp = 0;
    w = weights(order_max)
    x = roots(order_max)
    for q in range(0,w.size):
        tmp=tmp + w[q]*f(2**(-scale)*(x[q]+translation))*legendre_polynomial(order,float(x[q]))
    return 2**(-.5*scale)*tmp

def func_proj(scale,max_order,f,x):
    tmp = 0
    for j in range(0,max_order):
        for l in range(0,2**scale):
            tmp = tmp+scaling_coef(j,max_order,scale,l,f)*scaling(j,scale,l,x)
    return tmp
