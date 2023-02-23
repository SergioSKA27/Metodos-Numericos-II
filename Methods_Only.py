#!/usr/bin/python
import PySimpleGUI as sg
import numpy as np
import math
import sys
from sympy import *
from matplotlib import pyplot as plt
from sympy.plotting.plot import MatplotlibBackend, Plot

def jacobian(ff,symb):
    """
    It takes a vector of functions and a vector of symbols and returns the Jacobian matrix of the functions with respect to
    the symbols
    :param ff: the function
    :param symb: the symbols that are used in the function
    :return: A matrix of the partial derivatives of the function with respect to the variables.
    """
    m = []

    for i in range(0,len(ff)):
        aux  = []
        for j in range(0,len(symb)):
            aux.append(diff(ff[i],symb[j]))
        m.append(aux)

    return np.array(m)

def hessian(ff,symb):
    """
    It takes a vector of functions and a vector of symbols and returns the Hessian matrix of the functions with respect to
    the symbols
    :param ff: a list of functions of the form f(x,y,z)
    :param symb: the symbols that are used in the function
    :return: A matrix of the second derivatives of the function.
    """

    m = []

    for i in range(0,len(ff)):
        aux  = []
        for j in range(0,len(symb)):
            aux.append(diff(ff[i],symb[j],2))
        m.append(aux)
    return np.array(m)

def eval_matrix(matrix , v,symb):
    """
    It takes a matrix, a list of symbols and a list of values, and returns the matrix with the symbols substituted by the
    values

    :param matrix: the matrix of the system of equations
    :param v: the vector of values for the variables
    :param symb: the symbols that will be used in the matrix
    :return: the matrix with the values of the variables substituted by the values of the vector v.
    """
    e = 0
    mm = []
    for i in range(0,len(matrix)):
        aux = []
        ev = []
        for k in range(0,len(symb)):
            ev.append((symb[k],v[k]))
        for j in range(len(matrix[i])):
            aux.append(matrix[i][j].subs(ev).evalf())
        mm.append(aux)
    return np.array(mm)

def evalVector(ff, x0,symb):
    """
    > Given a list of symbolic expressions, a list of values for the symbols, and a list of the symbols, evaluate the
    symbolic expressions at the given values

    :param ff: the vector of functions
    :param x0: initial guess
    :param symb: the symbols that are used in the symbolic expression
    :return: the value of the function at the point x0.
    """
    v = []
    for i in range(0,len(ff)):
        ev = []

        for k in range(0,len(x0)):
            ev.append((symb[k],x0[k]))

        v.append(ff[i].subs(ev).evalf())
    return np.array(v)

def NewtonMethod( ff, x0,symb ):
    """
    The function takes in a vector of functions, a vector of initial guesses, and a vector of symbols. It then calculates
    the Jacobian matrix, the Jacobian matrix evaluated at the initial guess, the inverse of the Jacobian matrix evaluated at
    the initial guess, the vector of functions evaluated at the initial guess, and then the Newton step.

    The function returns the Newton step.

    :param ff: the function we want to find the root of
    :param x0: initial guess
    :param symb: the symbols used in the function
    :return: The return value is the x_np1 value.
    """
    j = jacobian(ff,symb)
    #print("Jacobian Matrix")
    #pprint(Matrix(j))
    jev = Matrix( eval_matrix(j,x0,symb))
    #print("J(",x0,")")
    #pprint(jev)

    jinv = jev.inv()
    #print("F(",x0,")")
    ffev = Matrix(evalVector(np.transpose(ff),x0,symb))
    #print("J^-1(",x0,")*","F(",x0,")")
    mm = Matrix(jinv)*ffev
    #pprint(mm)
    x_np1 = Matrix(np.transpose(np.array(x0)))
    #pprint(x_np1-mm)
    return list(x_np1-mm)

def norm_inf(x_0,x_1):
    """
    > The function `norm_inf` takes two vectors `x_0` and `x_1` and returns the maximum absolute difference between the two
    vectors

    :param x_0: the initial guess
    :param x_1: the vector of the current iteration
    :return: The maximum difference between the two vectors.
    """
    a = [abs(x_1[i]-x_0[i]) for i in range(len(x_0))]
    return max(a)

def newton_method(ff,x_0,symbs):
    """
    Given a function (x,y)$, a starting point $, and a list of symbols,
    the function will return the next point $ in the Newton's method sequence

    :param ff: the function to be minimized
    :param x_0: initial guess
    :param symbs: the symbols that we're using in the function
    :return: the final value of x_0, the list of x values, and the list of y values.
    """
    pprint(Matrix(x_0))
    xs = []
    ys = []

    while True:
        x_1 = NewtonMethod(ff,x_0,symbs)
        print(x_1)
        ninf = norm_inf(x_0,x_1)
        print(ninf)

        x_0 = list(x_1)
        xs.append(x_0[0])
        ys.append(x_0[1])
        if ninf < 1e-6:
            break

    print(x_0)
    return x_0,xs,ys

def get_sympy_subplots(plot:Plot):
    """
    It takes a plot object and returns a matplotlib figure object

    :param plot: The plot object to be rendered
    :type plot: Plot
    :return: A matplotlib figure object.
    """
    backend = MatplotlibBackend(plot)

    backend.process_series()
    backend.fig.tight_layout()
    return backend.plt

def li(v, i):
    """
    The function takes a list of numbers and an index, and returns the Lagrange interpolating polynomial for the list of
    numbers with the index'th number removed

    :param v: the list of x values
    :param i: the index of the x value you want to interpolate
    :return: the Lagrange interpolating polynomial for the given data points.
    """
    x = symbols('x')

    s = 1
    st = ''
    for k in range(0,len(v)):
        if k != i:
            st = st + '((' + str(x) + '-'+ str(v[k])+')/('+str(v[i])+'-'+str(v[k])+'))'
            s = s*((x-v[k])/(v[i]-v[k]))

    return s

def Lagrange(v,fx):
    """
    It takes in a list of x values and a list of y values, and returns the Lagrange polynomial that interpolates those
    points

    :param v: list of x values
    :param fx: The function you want to interpolate
    :return: the Lagrange polynomial.
    """
    print(v)
    print(fx)
    lis = []
    for i in range(0,len(v)):
        lis.append(li(v,i))

    sums = 0

    for k in range(0,len(v)):
        sums = sums+(fx[k]*lis[k])

    print(sums)

    simplify(sums)

    pprint(sums)

    p1 = plot(sums,(sums,0,math.pi),(math.pi/2,0),show=False)
    p2 = get_sympy_subplots(p1)
    p2.plot(v,fx,"o")
    p2.show()
    return sums

def diff_div(v,fx,order):
    """
    > The function takes in a list of values, a list of function values, and an order, and returns a list of divided
    differences

    :param v: the list of x values
    :param fx: the function you want to differentiate
    :param order: the order of the derivative you want to take
    :return: the difference quotient of the function f(x)
    """

    m = []

    for i in range(0,len(fx)):
        #print(fx[i])
        if i + 1 < len(fx) and i +order < len(v):
            #print(v[i+1],v[i],"/",fx[i+order]," ",fx[i])
            m.append((fx[i+1]-fx[i])/(v[i+order]-v[i]))
    return m

def divided_diff(fx,v):
    """
    The function takes in a list of x values and a list of f(x) values, and returns a list of lists of divided differences

    :param fx: the function to be interpolated
    :param v: list of x values
    :return: The divided difference table is being returned.
    """
    x = v
    nfx = fx
    m = []
    for i in range(0,len(v)-1):
        nx = diff_div(v,nfx,i+1)
        #print(nx)
        m.append(nx)
        nfx = nx

    #print(m)
    return m

def Newton_interpolation(fx,v):
    """
    It takes in a list of x values and a list of f(x) values, and returns a polynomial that interpolates the points

    :param fx: a list of the function values
    :param v: list of x values
    :return: The function is being returned.
    """
    diff = divided_diff(fx,v)
    x = symbols('x')

    expr = v[0]

    for i in range(0,len(diff)):
        s = diff[i][0]
        p = 1
        for k in range(0,len(v)):

            p = p*(x-v[k])
            #print(p, "p",k)
            if k == i:
                break
        s = s * p
        expr = expr + s

    pprint(expr)

    p = plot(expr,(x,-10,10),show=False)
    p2 = get_sympy_subplots(p)
    p2.plot(v,fx,"o")
    p2.show()

    return expr

def spline_natural(fx,v):
    """
    It takes a list of x values and a list of y values, and returns a list of sympy expressions that represent the cubic
    spline interpolation of the data

    :param fx: list of f(x) values
    :param v: list of x values
    """

    inter = []
    fxinter = []

    hi =[]
    for i in range(0,len(v)-1):
        inter.append((v[i],v[i+1]))
        fxinter.append((fx[i],fx[i+1]))

    print(inter)
    for i in range(0,len(inter)):
        hi.append(inter[i][1]-inter[i][0])

    m = np.zeros(len(v)**2).reshape(len(fx),len(fx))
    print(hi)
    print(m)
    for i in range(0,len(v)):
        for j in range(0,len(v)):
            if (i == j and i == 0 and j == 0) or (j == i and i == len(v)-1 and j == len(v)-1):
                m[i][j] = 1
                continue
            else:
                if (i == j):
                    m[i][j] = 2*(hi[i-1]+hi[i])
                    m[i][j-1] = hi[i-1]
                    m[i][j+1] = hi[i]

    b = np.zeros(len(v))

    for i in range(1,len(v)-1):
        b[i] = ((1/hi[i])*(fx[i+1]-fx[i]))-((1/hi[i-1])*(fx[i]-fx[i-1]))

    print(m)
    pprint(Matrix(b.transpose()))

    c = (Matrix(m).inv())*Matrix(b.transpose())
    pprint(c)
    b = []

    for i in range(0,len(hi)):
        b.append(((fx[i+1]-fx[i])/hi[i])-((((2*c[i])+c[i+1])*hi[i])/3))

    pprint(Matrix(b))

    d = []

    for i in range(0,len(hi)):
        d.append((c[i+1]-c[i])/(3*hi[i]))

    pprint(Matrix(d))


    x = symbols('x')
    spl = []
    for i in range(0,len(inter)):
        spl.append(fx[i]+ (b[i]*(x-v[i]))+(c[i]*((x-v[i])**2)) + (d[i]*((x-v[i])**3)))

    pprint(Matrix(spl))



    p = plot(spl[0], (x,inter[0][0],inter[0][1]),show=False)

    for i in range(1, len(spl)):
        paux = plot(spl[i],(x,inter[i][0],inter[i][1]),show=False)
        p.append(paux[0])


    p2 = get_sympy_subplots(p)
    p2.plot(v,fx,"o")
    p2.show()

def spline_sujeto(fx,v,fpx0,fpx1 ):
    """
    It takes a list of x values, a list of y values, and the first and second derivatives of the first and last points, and
    returns a plot of the cubic spline interpolation

    :param fx: the function values
    :param v: the x values of the points
    :param fpx0: the first derivative of the function at the first point
    :param fpx1: the derivative of the function at the last point
    """

    inter = []
    fxinter = []

    hi =[]
    for i in range(0,len(v)-1):
        inter.append((v[i],v[i+1]))
        fxinter.append((fx[i],fx[i+1]))

    print(inter)
    for i in range(0,len(inter)):
        hi.append(inter[i][1]-inter[i][0])

    m = np.zeros(len(v)**2).reshape(len(fx),len(fx))
    print(hi)
    print(m)
    for i in range(0,len(v)):
        for j in range(0,len(v)):
            if (i == j and i == 0 and j == 0) :
                m[i][j] = 2*hi[i]
                m[i][j+1] = hi[i]
                continue
            elif (j == i and i == len(v)-1 and j == len(v)-1):
                m[i][j] = 2*hi[-1]
                m[i][j-1] = hi[-1]
                continue
            else:
                if (i == j):
                    m[i][j] = 2*(hi[i-1]+hi[i])
                    m[i][j-1] = hi[i-1]
                    m[i][j+1] = hi[i]

    b = np.zeros(len(v))
    b[0] = ((3/hi[0])*(fx[1]-fx[0]))- (3*fpx0)
    b[-1] = (3*fpx1)-((3/hi[-1])*(fx[-1]-fx[len(fx)-2]))

    for i in range(1,len(v)-1):
        b[i] = ((3/hi[i])*(fx[i+1]-fx[i]))-((3/hi[i-1])*(fx[i]-fx[i-1]))

    print(m)
    pprint(Matrix(b.transpose()))

    c = (Matrix(m).inv())*Matrix(b.transpose())
    pprint(c)
    b = []

    for i in range(0,len(hi)):
        b.append(((fx[i+1]-fx[i])/hi[i])-((((2*c[i])+c[i+1])*hi[i])/3))

    pprint(Matrix(b))

    d = []

    for i in range(0,len(hi)):
        d.append((c[i+1]-c[i])/(3*hi[i]))

    pprint(Matrix(d))


    x = symbols('x')
    spl = []
    for i in range(0,len(inter)):
        spl.append(fx[i]+ (b[i]*(x-v[i]))+(c[i]*((x-v[i])**2)) + (d[i]*((x-v[i])**3)))

    pprint(Matrix(spl))



    p = plot(spl[0], (x,inter[0][0],inter[0][1]),show=False)

    for i in range(1, len(spl)):
        paux = plot(spl[i],(x,inter[i][0],inter[i][1]),show=False)
        p.append(paux[0])


    p2 = get_sympy_subplots(p)
    p2.plot(v,fx,"o")
    p2.show()






#x,y = symbols('x,y')
#
#s = np.array([x,y])
#ff = np.array([x**3+3*y**2-21,x**2+2*y+2])
#
#p = plot_implicit(ff[0],(x,-10,10),(y,-10,10),show=False)
#p.append((plot_implicit(ff[1],(x,-10,10),(y,-10,10),show=False))[0])
#p2 = get_sympy_subplots(p)
#
#x = newton_method(ff,[1,-1],s)
#print(x)
#p2.plot(x[1],x[2],"o")
#p2.show()
#
#
#vv = [0,math.pi/2,math.pi]
#fxi = [0,1,0]


#Lagrange(vv,fxi)


xx = [2,4,6,8]
fxx = [4,8,14,16]
#Newton_interpolation(fxx,xx)

spline_natural([0,0.5,2,1.5],[0,1,2,3])

spline_sujeto([0,0.5,2,1.5],[0,1,2,3],1/5.0,-1)
