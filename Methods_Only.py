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



x,y = symbols('x,y')

s = np.array([x,y])
ff = np.array([x**3+3*y**2-21,x**2+2*y+2])

p = plot_implicit(ff[0],(x,-10,10),(y,-10,10),show=False)
p.append((plot_implicit(ff[1],(x,-10,10),(y,-10,10),show=False))[0])
p2 = get_sympy_subplots(p)

x = newton_method(ff,[1,-1],s)
print(x)
p2.plot(x[1],x[2],"o")
p2.show()


vv = [0,math.pi/2,math.pi]
fxi = [0,1,0]


Lagrange(vv,fxi)
