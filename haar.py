
import numpy as np
import matplotlib.pyplot as plt

def gauss(x):
    return np.exp(-100*(x-.5)**2)

def chi(x):
    if type(x) is int or float == True:
        if x < 0 or x > 1:
            return 0
        else:
            return 1
    else:
        y = np.empty(x.shape)
        for i in range(x.size):
            if x[i] < 0 or x[i] > 1:
                y[i]=0
            else:
                y[i]=1
        return y


def haar(x):
    return chi(2*x) - chi(2*x-1)

def chijk(j,k,x):
    return 2**(j/2)*chi(2**j * x -k)

def haarjk(j,k,x):
    return 2**(j/2)*haar(2**j * x -k)


def snk_start(f,n,k):
    return 2**(n/2)* np.trapz(f(  np.arange(2**(-n)*k,2**(-n)*(k+1),0.00001) ),dx= 0.00001)

def snk(f,n,k):
    return 1/np.sqrt(2)*(snk_start(f,n+1,2*k)+snk_start(f,n+1,2*k+1))

def dnk(f,n,k):
    return 1/np.sqrt(2)*(snk_start(f,n+1,2*k)-snk_start(f,n+1,2*k+1))

def f_proj(f,x,n):
    tmp = chijk(n,0,x)*snk_start(f,n,0)

    for k in range(1,2**n):
        tmp = tmp + chijk(n,k,x)*snk(f,n,k)
    return tmp


x=np.arange(0.0, 1, 0.003)


plt.plot(x,f_proj(gauss,x,10),x,gauss(x))
plt.show()
