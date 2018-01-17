import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.legendre as leg

def legendre_polynomial(order, x):
    coef = []
    if type(x) is float:
       # print("here")
        if x<0 or x>1:
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
        print(y.shape)
        for i in range(x.size):
           # print(i)
            if x[i] < 0 or x[i] > 1:
                y[i] = 0
            else:
                for n in range(0, order+1):
                    if(n == order):
                        coef.append(1)
                    else:
                        coef.append(0)
               # print(coef)
                #print(x[i])
                #print(y[i])
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
