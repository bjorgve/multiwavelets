import numpy, sys, struct
from math import sin, pi
import matplotlib.pyplot as plt

# luca.frediani@uit.no

nPolynomials = 4
max_scale = 5

#
# Given an index j and a position x, returns the value of the j-th scaling function in x
# Tips: use the legendre polynomial defined by numpy
#
def scaling(j, x):
    coef = []
    for n in range(0, j+1):
        if(n == j):
            coef.append(1)
        else:
            coef.append(0)
    factor = numpy.sqrt(2*j+1)
    if(x >= 0 and x < 1):
        polynom = numpy.polynomial.legendre.Legendre(coef, domain=[0,1])
        value = polynom(x)
    else:
        value = 0.0
    return factor * value

#
# Given the index j, the number of polynomials np, the position x,
# returns the value of the j-th wavelet function of order np in x
#
def wavelet(j,np,x):
    Sum = 0
    G0 = getfilter("G", 0, np)
    G1 = getfilter("G", 1, np)
#    print G0
#    print G1
    for i in range(np):
        wavelet_formel = (G0[j, i]*scaling_nl(i, 1, 0, x))/numpy.sqrt(2) - (G1[j, i]*scaling_nl(i, 1, 1, x))/numpy.sqrt(2)
        Sum += wavelet_formel
    return Sum

#
# Given the number of polynomials np and the position x,
# returns the value of all wavelet functions of order np in x
#
def wavelets(np,x):
    print("FIXME")

#
# Given the number of polynomials np, the scale n, the translation l and the position x,
# returns the value of all wavelet functions of order np in x at scale n and translation l
#
def wavelets_nl(n,l,np,x):
    print("FIXME")

#
# Given the index j, the scale n, the translation l and the position x,
# returns the value of the scaling function j at scale n and translation l
#
def scaling_nl(j, n, l, x):
    factor = 2**n
    return numpy.sqrt(factor) * scaling(j,factor*x-l)
#
# Given the index j, the number of polynomials np, the scale n, the translation l and the position x,
# returns the value of the j-th wavelet function of order np in x at scale n and translation l
#
def wavelet_nl(j, n, l, np, x):
#    return 2**(n/2)*wavelet(j,np,x)*((2**n)*x-l)
        return 2**(n/2)*wavelet(j,np,(2**n)*x-l)

#
# Load the requested filter matrix (H/G and 0/1) for a given order np
#
def getfilter(type="H",flag=0,np=4):
    order = np - 1
    name = "mwfilters/L_" + type + "0_" + str(order)
    filterFile = open(name, 'rb')
    filterData = filterFile.read()
    filterList = []
    for i in range((np)*(np)):
        rstart = i * 8
        rend = (i + 1) * 8
        filterList.append(struct.unpack_from("d",filterData[rstart:rend]))
    filterArray = numpy.array(filterList)
    filterMatrix = numpy.ndarray((np, np),buffer=filterArray)
    if (flag == 1):
        if (type == "G" and np%2 != 0 ):
            sign_np = -1
        else:
            sign_np = 1
        for i in range(np):
            for j in range(np):
                sign_ij = 1
                if (((i+j)%2) == 0):
                    sign_ij = -1
                sign = sign_ij * sign_np
                filterMatrix[i,j] *= sign
    return filterMatrix



if __name__ == '__main__':
    x = numpy.arange(-0.4, 1.4001, 0.001);
    y0 = numpy.array([wavelet_nl(0,2,0,4,x[i]) for i in range(x.size)])
    y1 = numpy.array([wavelet_nl(1,2,1,4,x[i]) for i in range(x.size)])
    y2 = numpy.array([wavelet_nl(2,2,2,4,x[i]) for i in range(x.size)])
    y3 = numpy.array([wavelet_nl(3,2,3,4,x[i]) for i in range(x.size)])
#    w0 = numpy.array([wavelet_nl(0,3,0,8,x[i]) for i in range(x.size)])
#    w1 = numpy.array([wavelet_nl(1,3,1,8,x[i]) for i in range(x.size)])
#    w2 = numpy.array([wavelet_nl(2,3,2,8,x[i]) for i in range(x.size)])
#    w3 = numpy.array([wavelet_nl(3,3,3,8,x[i]) for i in range(x.size)])
#    w4 = numpy.array([wavelet_nl(4,3,4,8,x[i]) for i in range(x.size)])
#    w5 = numpy.array([wavelet_nl(5,3,5,8,x[i]) for i in range(x.size)])
#    w6 = numpy.array([wavelet_nl(6,3,6,8,x[i]) for i in range(x.size)])
#    w7 = numpy.array([wavelet_nl(7,3,7,8,x[i]) for i in range(x.size)])
#    plt.plot(x, w0)
#    plt.plot(x, w1)
#    plt.plot(x, w2)
#    plt.plot(x, w3)
#    plt.plot(x, w4)
#    plt.plot(x, w5)
#    plt.plot(x, w6)
    plt.plot(x, y0)
    plt.plot(x, y1)
    plt.plot(x, y2)
    plt.plot(x, y3)
    plt.show()
