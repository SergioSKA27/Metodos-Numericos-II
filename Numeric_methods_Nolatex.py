#!/usr/bin/python
import PySimpleGUI as sg
import numpy as np
import math
import sys
from sympy import *
from matplotlib import pyplot as plt
from sympy.plotting.plot import MatplotlibBackend, Plot
from numpy.linalg import *
import copy
from sympy.plotting import plot3d,plot3d_parametric_line
from tkinter import END, INSERT
from tkinter import *
from os import system
import platform
from sympy.plotting import textplot

if platform.system() == 'Linux':
    CLEARW = 'clear'
elif platform.system() ==  'Windows':
    CLEARW = 'cls'

#This version don't needs latex to work properly :)

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

def NewtonMethod( ff, x0,symb):
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

    #printing = []
    j = jacobian(ff,symb)
    print("Jacobian Matrix")
    pprint(Matrix(j),use_unicode=False)
    jev = Matrix( eval_matrix(j,x0,symb))
    print("J(",x0,")")
    pprint(jev,use_unicode=False)


    jinv = jev.inv()

    print("F(",x0,")")

    ffev = Matrix(evalVector(np.transpose(ff),x0,symb))
    print("J^-1(",x0,")*","F(",x0,")")

    mm = Matrix(jinv)*ffev
    pprint(mm,use_unicode=False)

    x_np1 = Matrix(np.transpose(np.array(x0)))
    print("X_n+1 = ")
    pprint(Matrix(x_np1-mm),use_unicode=False)

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
    system(CLEARW)

    print("METODO DE NEWTON SISTEMAS DE ECUACIONES NO LINEALES")

    print("Aproximación Inicial:")
    pprint(Matrix(x_0),use_unicode=False)
    xs = []
    ys = []
    itercount = 0

    while True and itercount < 1000:

        print("ITERACION ", itercount)
        x_1 = NewtonMethod(ff,x_0,symbs)
        #print(x_1)
        ninf = norm_inf(x_0,x_1)
        print("||X_n+1 - X_n || =",ninf)
        #print(ninf)


        x_0 = list(x_1)
        xs.append(x_0[0])
        ys.append(x_0[1])
        if ninf < 1e-6:
            print("||X_n+1 - X_n|| < error (True)")
            break
        else:
            print("||X_n+1 - X_n|| < error (False)")
        itercount = itercount + 1

    print("La Solucion es: ")
    pprint(Matrix(x_0),use_unicode=False)

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
    system(CLEARW)

    print("INTERPOLACION POLINOMICA DE LAGRANGE ")

    print("x = ")
    pprint(Matrix(v),use_unicode=False)
    print("f(x) = ")
    pprint(Matrix(fx),use_unicode=False)
    lis = []
    for i in range(0,len(v)):
        lis.append(li(v,i))

    print("Polinomios de Lagrange (Li)")
    for i in range(0,len(lis)):
        print("L_"+str(i)+" = ")
        pprint(lis[i],use_unicode=False)

    sums = 0

    for k in range(0,len(v)):
        sums = sums+(fx[k]*lis[k])

    print("El polinomios de lagrange esta dado por: ")

    simplify(sums)

    pprint(sums,use_unicode=False)

    p1 = textplot(sums,min(v),max(v))
    #p2 = get_sympy_subplots(p1)
    #p2.plot(v,fx,"o")
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

    system(CLEARW)

    print("INTERPOLACION POLINOMICA DE NEWTON")
    print("x = ")
    pprint(Matrix(v),use_unicode=False)
    print("f(x) = ")
    pprint(Matrix(fx),use_unicode=False)
    diff = divided_diff(fx,v)
    x = symbols('x')

    expr = v[0]

    print("Diferencias Divididas: ")
    for i in diff:
        pprint(Matrix(i),use_unicode=False)

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

    print("El Polinomio de interpolacion esta dado por: ")
    pprint(expr,use_unicode=False)

    p = textplot(expr,min(v),max(v))
    #p2 = get_sympy_subplots(p)
    #p2.plot(v,fx,"o")
    #p2.show()

    return expr

def spline_natural(fx,v):
    """
    It takes a list of x values and a list of y values, and returns a list of sympy expressions that represent the cubic
    spline interpolation of the data

    :param fx: list of f(x) values
    :param v: list of x values
    """
    system(CLEARW)

    print("Interpolacion por medio de Splines Cubicos(Natural)")

    print("x = ")
    pprint(Matrix(v),use_unicode=False)
    print("f(x) = ")
    pprint(Matrix(fx),use_unicode=False)

    inter = []
    fxinter = []

    hi =[]
    for i in range(0,len(v)-1):
        inter.append((v[i],v[i+1]))
        fxinter.append((fx[i],fx[i+1]))

    #print(inter)
    for i in range(0,len(inter)):
        hi.append(inter[i][1]-inter[i][0])

    m = np.zeros(len(v)**2).reshape(len(fx),len(fx))
    #print(hi)
    #print(m)
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

    print("A = ")
    pprint(Matrix(m),use_unicode=False)

    b = np.zeros(len(v))

    for i in range(1,len(v)-1):
        b[i] = ((1/hi[i])*(fx[i+1]-fx[i]))-((1/hi[i-1])*(fx[i]-fx[i-1]))

    print("b0 = ")
    pprint(Matrix(b.transpose()),use_unicode=False)

    print("A*b0 = c ---> c = ")
    c = (Matrix(m).inv())*Matrix(b.transpose())
    pprint(Matrix(c),use_unicode=False)
    b = []

    for i in range(0,len(hi)):
        b.append(((fx[i+1]-fx[i])/hi[i])-((((2*c[i])+c[i+1])*hi[i])/3))

    print("b = ")

    pprint(Matrix(b),use_unicode=False)

    d = []

    for i in range(0,len(hi)):
        d.append((c[i+1]-c[i])/(3*hi[i]))

    print("d = ")
    pprint(Matrix(d),use_unicode=False)


    x = symbols('x')
    spl = []
    for i in range(0,len(inter)):
        spl.append(fx[i]+ (b[i]*(x-v[i]))+(c[i]*((x-v[i])**2)) + (d[i]*((x-v[i])**3)))

    print("Los Splines son: ")
    pprint(Matrix(spl),use_unicode=False)


    print("SPLINE_0 = ")
    pprint((spl[0]),use_unicode=False)
    p = textplot(spl[0], inter[0][0],inter[0][1])

    for i in range(1, len(spl)):
        print("SPLINE_"+str(i)+" = ")
        pprint(spl[i],use_unicode=False)
        paux = textplot(spl[i],inter[i][0],inter[i][1])



    #p2 = get_sympy_subplots(p)
    #p2.plot(v,fx,"o")
    #p2.show()
    return spl

def spline_sujeto(fx,v,fpx0,fpx1 ):
    """
    It takes a list of x values, a list of y values, and the first and second derivatives of the first and last points, and
    returns a plot of the cubic spline interpolation

    :param fx: the function values
    :param v: the x values of the points
    :param fpx0: the first derivative of the function at the first point
    :param fpx1: the derivative of the function at the last point
    """
    system(CLEARW)

    print("Interpolacion por medio de Splines Cubicos(Sujeto)")

    print("x = ")
    pprint(Matrix(v),use_unicode=False)
    print("f(x) = ")
    pprint(Matrix(fx),use_unicode=False)

    print("f'(a) =", fpx0)
    print("f'(b) =", fpx1)


    inter = []
    fxinter = []

    hi =[]
    for i in range(0,len(v)-1):
        inter.append((v[i],v[i+1]))
        fxinter.append((fx[i],fx[i+1]))

    #print(inter)
    for i in range(0,len(inter)):
        hi.append(inter[i][1]-inter[i][0])

    m = np.zeros(len(v)**2).reshape(len(fx),len(fx))
    #print(hi)
    #print(m)
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

    print("A =")
    pprint(Matrix(m),use_unicode=False)
    print("b0 = ")
    pprint(Matrix(b.transpose()),use_unicode=False)

    print("A*b0 = c ---> c = ")

    c = (Matrix(m).inv())*Matrix(b.transpose())
    pprint(Matrix(c),use_unicode=False)
    b = []

    for i in range(0,len(hi)):
        b.append(((fx[i+1]-fx[i])/hi[i])-((((2*c[i])+c[i+1])*hi[i])/3))

    print("b = ")
    pprint(Matrix(b),use_unicode=False)

    d = []

    for i in range(0,len(hi)):
        d.append((c[i+1]-c[i])/(3*hi[i]))

    print("d = ")
    pprint(Matrix(d),use_unicode=False)


    x = symbols('x')
    spl = []
    for i in range(0,len(inter)):
        spl.append(fx[i]+ (b[i]*(x-v[i]))+(c[i]*((x-v[i])**2)) + (d[i]*((x-v[i])**3)))

    print("Los Splines son: ")
    pprint(Matrix(spl),use_unicode=False)



    print("SPLINE_0 = ")
    pprint((spl[0]),use_unicode=False)
    p = textplot(spl[0], inter[0][0],inter[0][1])

    for i in range(1, len(spl)):
        print("SPLINE_"+str(i)+" = ")
        pprint(spl[i],use_unicode=False)
        paux = textplot(spl[i],inter[i][0],inter[i][1])


    #p2 = get_sympy_subplots(p)
    #p2.plot(v,fx,"o")
    #p2.show()

    return spl

def discrete_minimun_quads_aprox(xs,y,functionss,symbs):
    """
    Given a set of points $(x_i,y_i)$ and a set of functions (x)$, it finds the linear combination of the functions that
    best fits the points

    :param xs: list of x values
    :param y: the y values of the data points
    :param functionss: a list of functions that will be used to approximate the data
    :param symbs: the symbol that you want to use for the function
    :return: The expression of the function that best fits the data.
    """

    system(CLEARW)
    print("APROXIMACION DISCRETA POR MINIMOS CUADRADOS")

    m = []


    for i in range(0,len(xs)):
        aux = []

        for j in range(0,len(functionss)):

            aux.append(functionss[j])
        m.append(aux)

    print("M = ")
    #pprint(Matrix(m))

    mev = []
    for i in range(0,len(m)):
        aux = []

        for j in range(0,len(m[0])):
            if len(m[i][j].free_symbols) > 0:
                aux.append(m[i][j].subs(symbs,xs[i]))
            else:
                aux.append(m[i][j])
        mev.append(aux)

    pprint(Matrix(mev),use_unicode=False)

    mevT = Matrix(mev).transpose()
    print("M^T = ")
    pprint(mevT,use_unicode=False)

    a = mevT*Matrix(mev)

    print("A = M^T * M")
    pprint(a,use_unicode=False)

    b = mevT*Matrix(y)


    print("b = M^T * y")
    pprint(b,use_unicode=False)

    ainv = a.inv()
    print("Ax = b ---> x = ")
    xsol = ainv*b

    pprint(xsol,use_unicode=False)

    print("El la recta esta dada esta dado por: ")
    expr = xsol[0]+xsol[1]*symbs
    pprint(expr,use_unicode=False)


    p = textplot(expr,min(xs),max(xs))

    return expr

def continium_minimun_quads_aprox(fx,interval,symb,degree):
    """
    Given a function, an interval, a symbol and a degree, it returns the polynomial that best approximates the function in
    the given interval

    :param fx: The function to be approximated
    :param interval: the interval in which the function is defined
    :param symb: the symbol that will be used to represent the variable in the function
    :param degree: The degree of the polynomial
    :return: The function that is the best aproximation of the given function in the given interval.
    """
    system(CLEARW)
    print("APROXIMACION DISCRETA POR MINIMOS CUADRADOS")
    m = []


    for i in range(0,degree+1):
        aux = []
        for j in range(0,degree+1):
            aux.append(integrate((symb**i)*(symb**j),(symb,interval[0],interval[1])))
        m.append(aux)

    print("M =  ")
    pprint(Matrix(m),use_unicode=False)


    b = []

    for i in range(0,degree+1):
        b.append(integrate((symb**i)*fx,(symb,interval[0],interval[1])))

    print("b = ")
    pprint(Matrix(b))

    print("M*x = b ---> x = ")

    sol = Matrix(m).inv() * Matrix(b)
    pprint(Matrix(sol),use_unicode=False)
    expr = 0

    for i in range(0,degree+1):
        expr = expr + (sol[i]*symb**i)

    print("La solucion esta dada por \n")
    pprint(expr,use_unicode=False)

    pprint(fx,use_unicode=False)
    p = textplot(fx,interval[0],interval[1])
    pprint(expr,use_unicode=False)
    p2 = textplot(expr,interval[0],interval[1])



    return expr



def is_numeric(sn):
    """
    If the input is a string that is a single digit, return True, otherwise return False

    :param sn: the string to be checked
    :return: True or False
    """
    if sn=='0' or sn=='1' or sn=='2' or sn=='3' or sn=='4' or sn=='5' or sn=='6' or sn=='7' or sn == '8' or sn== '9':
        return True
    else:
        return False

def is_operator(sn):
    """
    It checks if the input is an operator.

    :param sn: the string to be checked
    :return: True or False
    """
    if sn=='+' or sn=='-' or sn=='*' or sn=='/' or sn=='=':
        return True
    else:
        return False

def fix_expression(s):
    """
    It takes a string, removes all spaces, replaces all '^' with '**', and adds a negative sign to the right side of the
    equation

    :param s: the string to be fixed
    :return: the expression in the form of a string.
    """
    sn = ''

    for i in range(0,len(s)):
        if s[i] != ' ':
            sn = sn + str(s[i])

    #print(sn)
    sfix = ''
    eq = ''
    f = True
    for i in range(0,len(sn)):
        if sn[i] == '=':
            f = False
            continue
        if f == False:
            eq = eq + sn[i]
        if sn[i] == '^' and f:
            sfix = sfix + "**"
        elif f:
            sfix = sfix + sn[i]

    #print(sfix)
    #print(eq)

    return sfix+str(-1*float(eq))

def parse_system(st):
    """
    It takes a string of equations and returns a list of strings, each string being an equation

    :param s: the string to be parsed
    :return: A list of strings
    """
    eq = []
    s = ''

    for i in range(0,len(st)):
        if st[i] != '\n':
            s = s + st[i]
        else:
            eq.append(fix_expression(copy.copy(s)))
            s = ''
    if len(s) > 0:
        eq.append(fix_expression(s))

    return eq


def parse_interpolation(sinter):
    """
    It takes a string of the form "a = b\na = b\n" and returns a list of tuples of the form [(a,b),(a,b)]

    :param sinter: the string to be parsed
    :return: A list of tuples.
    """
    params = []
    s1 = ''
    s2 = ''
    flag = False
    for i in range(0,len(sinter)):
        if sinter[i] == ' ':
            continue
        if sinter[i] == '=':
            flag = True
            continue
        if sinter[i] == '\n':
            params.append((s1,s2))
            s1 = ''
            s2 = ''
            flag = False
            continue

        if flag == False:
            s1 = s1+sinter[i]
        else:
            s2 = s2 +sinter[i]

    if len(s1) > 0 and len(s2) > 0:
            params.append((s1,s2))

    return params


def parse_array(sarray):
    """
    It takes a string of numbers separated by commas and returns a list of floats.

    :param sarray: the string array that contains the data
    :return: A list of floats
    """
    sn = ''

    for i in range(0,len(sarray)):
        if sarray[i] != ' ' and (sarray[i] != '[' and sarray[i] != ']') and (sarray[i] != '{' and sarray[i] != '}'):
            sn = sn + sarray[i]

    return list(map(float,sn.split(',')))


lisimageL10_1 = []
countL10_1 = 0

def add_monkeys(window, counter,outputim,filename):
    tktext = window[outputim].widget
    global my_image
    my_image = PhotoImage(file=filename)
    lisimageL10_1.append(my_image)
    position = tktext.index(INSERT)
    tktext.image_create(position, image=lisimageL10_1[counter-1])
    #print(str(my_image))

lisimageL11_1 = []
countL11_1 = 0

def add_monkeys2(window, counter,outputim,filename):
    tktext = window[outputim].widget
    global my_image2
    my_image2 = PhotoImage(file=filename)
    lisimageL11_1.append(my_image2)
    position = tktext.index(INSERT)
    tktext.image_create(position, image=lisimageL11_1[counter-1])
    #print(str(my_image))


lisimageL11_2 = []
countL11_2 = 0

def add_monkeys3(window, counter,outputim,filename):
    tktext = window[outputim].widget
    global my_image3
    my_image3 = PhotoImage(file=filename)
    lisimageL11_2.append(my_image3)
    position = tktext.index(INSERT)
    tktext.image_create(position, image=lisimageL11_2[counter-1])
    #print(str(my_image))



lisimageL11_4 = []
countL11_4 = 0

def add_monkeys4(window, counter,outputim,filename):
    tktext = window[outputim].widget
    global my_image4
    my_image4 = PhotoImage(file=filename)
    lisimageL11_4.append(my_image4)
    position = tktext.index(INSERT)
    tktext.image_create(position, image=lisimageL11_4[counter-1])
    #print(str(my_image))



layout0 = [
[sg.Text('MENU PRINCIPAL')],

[sg.Button('Solución de  Sistemas Ecuaciones Lineales',size=(10,5),key='-LinearEqButton-'),
sg.Text('          '),
sg.Button('Solución de  Sistemas Ecuaciones No Lineales',size=(10,5),key='-NoLinearEqButton-'),
sg.Text('          '),
sg.Button('Métodos de Interpolación',size=(10,5),key='-InterpolationButton-')],

[sg.Text('')],

[sg.Button('Métodos de Aproximación',size=(10,5),key='-AproximationsButton-'),
sg.Text('          '),
sg.Button('Derivación Numérica',size=(10,5),key='-NumDerivationButton-'),
sg.Text('          '),
sg.Button('Integración Numérica',size=(10,5),key='-NumIntegrationButton-')]
]



layout10 = [
[sg.Button('<--',size=(2,2),key='-returnL0-L10-'),sg.Text('MENU SOLUCIÓN DE SISTEMAS DE ECUACIONES NO LINEALES')],

[sg.Button('Metodo de Newton Para Sistemas de Ecuaciones No Lineales',size=(10,5),key='-NoLineaSystemNewton-')],
]



layout10_1 =  [

[sg.Button('<--',size=(2,2),key='-returnL10-L10_1-'),
sg.Text('Método de Newton Para Sistemas de Ecuaciones No Lineales')],

[sg.Text("Ingrese el sistema de Ecuaciones a Resolver\t\t"),sg.Text('Registro del proceso'),
 sg.Button('SOLVE',key='-SolveL10_1-',size=(5,5)),sg.Button('Graficar',key='-PlotL10_1-',size=(5,5))],
[sg.Text("Valor inicial"),sg.Input("X_0",key='-x0L10_1-',size=(10,10) ),
 sg.Text("Error"),sg.Input("X_0",key='-errorL10_1-',size=(10,10) ),
 sg.Text("plot range"),sg.Input("X_range",key='-x_rangeplotL10_1-',size=(10,10) ),
 sg.Input("Y_range",key='-y_rangeplotL10_1-',size=(10,10) )],

[sg.Multiline("x**3+3*y**2-21 = 0\nx**2+2*y+2 = 0",key='-NoLinearSystemL10_1-',size=(30,30),horizontal_scroll=True),
 sg.Multiline("",key='-LogNewtonNoLinear-',size=(80,20),horizontal_scroll=True,auto_refresh=True)],

]


layout11 = [
[sg.Button('<--',size=(2,2),key='-returnL0-L11-'),sg.Text('MENU INTERPOLACIÓN')],

[sg.Button('Interpolación Polinomica de Lagrange',size=(10,5),key='-LagrangeInterpol-')],
[sg.Button('Interpolación Polinomica de Newton',size=(10,5),key='-NewtonInterpol-')],
[sg.Button('Interpolación Polinomica de Hermite',size=(10,5),key='-HermiteInterpol-')],
[sg.Button('Interpolación Mediante Splines Cubicos',size=(10,5),key='-SplineInterpol-')],
]


layout11_1 =  [

[sg.Button('<--',size=(2,2),key='-returnL11-L11_1-'),
sg.Text('Interpolación de Lagrange')],

[sg.Text("Ingrese los puntos x y f(x) \t\t"),sg.Text('Registro del proceso'),
 sg.Button('SOLVE',key='-SolveL11_1-',size=(5,5))],


[sg.Multiline("x = [0,1.570796327,3.141592654] \n f(x) = [0,1,0]",key='-DataInterpolationL11_1-',size=(30,30),horizontal_scroll=True),
 sg.Multiline("",key='-LogLagrangeInter-',size=(80,20),horizontal_scroll=True,auto_refresh=True)],

]



layout11_2 =  [

[sg.Button('<--',size=(2,2),key='-returnL11-L11_2-'),
sg.Text('Interpolación de Newton')],

[sg.Text("Ingrese los puntos x y f(x) \t\t"),sg.Text('Registro del proceso'),
 sg.Button('SOLVE',key='-SolveL11_2-',size=(5,5))],


[sg.Multiline("x = [2,4,6,8] \n f(x) = [4,8,14,16]",key='-DataInterpolationL11_2-',size=(30,30),horizontal_scroll=True),
 sg.Multiline("",key='-LogNewtonInter-',size=(80,20),horizontal_scroll=True,auto_refresh=True)],

]

layout11_3 =  [

[sg.Button('<--',size=(2,2),key='-returnL11-L11_3-'),
sg.Text('Interpolación de Hermite')],

[sg.Text("Ingrese los puntos x y f(x) \t\t"),sg.Text('Registro del proceso'),
 sg.Button('SOLVE',key='-SolveL11_3-',size=(5,5))],


[sg.Multiline("x = [1,2,3,4] \n f(x) = [1,2,4,5,6]",key='-DataInterpolationL11_3-',size=(30,30),horizontal_scroll=True),
 sg.Multiline("",key='-LogHermiteInter-',size=(80,20),horizontal_scroll=True,auto_refresh=True)],

]


layout11_4 =  [

[sg.Button('<--',size=(2,2),key='-returnL11-L11_4-'),
sg.Text('Interpolación Mediante Splines Cubicos')],

[sg.Text("Ingrese los puntos x y f(x) \t\t"),sg.Text('Registro del proceso'),sg.Text("Tipo de Spline:"),
 sg.Combo(['Natural','Sujeto'],default_value='Natural',key='-SplineType-',enable_events=True),
 sg.Button('SOLVE',key='-SolveL11_4-',size=(5,5))],


[sg.Multiline("x = [0,1,2,3] \n f(x) = [0,0.5,2,1.5]",key='-DataInterpolationL11_4-',size=(30,30),horizontal_scroll=True),
 sg.Multiline("",key='-LogSplineInter-',size=(80,20),horizontal_scroll=True,auto_refresh=True)],

]



#Layout menu de Metodos de Aproximacion
layout12 = [
[sg.Button('<--',size=(2,2),key='-returnL0-L12-'),sg.Text('MENU METODOS DE APROXIMACION')],

[sg.Button('Aproximacion discreta de minimos cuadrados',size=(10,5),key='-DisminQuadAprox-')],
[sg.Button('Aproximacion continua de minimos cuadrados',size=(10,5),key='-ContminQuadAprox-')],
]


layout12_1 =  [

[sg.Button('<--',size=(2,2),key='-returnL12-L12_1-'),
sg.Text('Aproximacion discreta de minimos cuadrados')],

[sg.Text("Ingrese los puntos x e y las funciones f_i(x) \t\t"),sg.Text('Registro del proceso'),
 sg.Button('SOLVE',key='-SolveL12_1-',size=(5,5))],


[sg.Multiline("x = [-1,1,3] \n y = [6,1,11]\n f_1(x) = 1 \n f_2(x) = x",key='-DataAproximationL12_1-',size=(30,30),horizontal_scroll=True),
 sg.Multiline("",key='-LogDisMinQuad-',size=(80,20),horizontal_scroll=True,auto_refresh=True)],

]


layout12_2 =  [

[sg.Button('<--',size=(2,2),key='-returnL12-L12_2-'),
sg.Text('Aproximacion continua de minimos cuadrados')],

[sg.Text("Ingrese la funcion y el intervalo \t\t"),sg.Text('Registro del proceso'),
 sg.Button('SOLVE',key='-SolveL12_2-',size=(5,5))],


[sg.Multiline("f(x) = sin(3.141592*x) \n I = [0,1]\n D = 2",key='-DataAproximationL12_2-',size=(30,30),horizontal_scroll=True),
 sg.Multiline("",key='-LogContMinQuad-',size=(80,20),horizontal_scroll=True,auto_refresh=True)],

]



layout = [
    [sg.Column(layout=layout0,key='-COL{0}-',visible=True),

     sg.Column(layout=layout10,key='-COL{10}-',visible=False),
     sg.Column(layout=layout10_1,key='-COL{101}-',visible=False),

     sg.Column(layout=layout11,key='-COL{11}-',visible=False),
     sg.Column(layout=layout11_1,key='-COL{111}-',visible=False),
     sg.Column(layout=layout11_2,key='-COL{112}-',visible=False),
     sg.Column(layout=layout11_3,key='-COL{113}-',visible=False),
     sg.Column(layout=layout11_4,key='-COL{114}-',visible=False),

     sg.Column(layout=layout12,key='-COL{12}-',visible=False),
     sg.Column(layout=layout12_1,key='-COL{121}-',visible=False),
     sg.Column(layout=layout12_2,key='-COL{122}-',visible=False),

    ]
]






#fix_expression("2x^2 + 3yx - y^2 = 0")

window = sg.Window('Métodos Numéricos ', layout,size=(720,480),resizable=True)

matop = '+'
#Event Loop to process "eSujetovents" and get the "values" of the inputs
while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'Cancel': # if user closes window or clicks cancel
        break


    if event == '-NoLinearEqButton-':
        window['-COL{0}-'].update(visible=False)
        window['-COL{10}-'].update(visible=True)

#--------------------------------------------------

    if event == '-NoLineaSystemNewton-':
        window['-COL{10}-'].update(visible=False)
        window['-COL{101}-'].update(visible=True)

    if event == '-PlotL10_1-':
        if len(values['-NoLinearSystemL10_1-']) > 0:
            try:
                eq = parse_system(values['-NoLinearSystemL10_1-'])
                print(eq)
                expr = []

                for i in range(0,len(eq)):
                    expr.append(parse_expr(eq[i]))

                print(expr)

                symbs = list(expr[0].free_symbols)

                if len(symbs) == 2:



                    try:
                        xs = list(map(float,values['-x_rangeplotL10_1-'].split(',')))
                        ys = list(map(float,values['-y_rangeplotL10_1-'].split(',')))
                        xr_1 = xs[0]
                        xr_2 = xs[1]
                        yr_1 = ys[0]
                        yr_2 = ys[1]
                        zr_1 = -10
                        zr_2 = 10
                    except:
                        xr_1 = -10
                        xr_2 = 10
                        yr_1 = -10
                        yr_2 = 10
                        zr_1 = -10
                        zr_2 = 10


                    try:

                        pl = plot_implicit(expr[0],(symbs[0],xr_1,xr_2),(symbs[1],yr_1,yr_2),show=False)
                        pl.append(plot_implicit(expr[1],(symbs[0],xr_1,xr_2),(symbs[1],yr_1,yr_2),show=False)[0])
                        pt3d = plot3d(expr[0],(symbs[0],xr_1,xr_2),(symbs[1],yr_1,yr_2),show=False)
                        pt3d.append(plot3d(expr[1],(symbs[0],xr_1,xr_2),(symbs[1],yr_1,yr_2),show=False)[0])
                        pt3d.show()
                        pl.show()
                    except:
                        sg.popup_ok("El sistema no se puede graficar :(")



            except:
                sg.popup_ok("Algo Salio Mal intente otra vez :(")

    if event == '-SolveL10_1-':
        try:
            window.Element('-LogNewtonNoLinear-').update(value="")
            eq = parse_system(values['-NoLinearSystemL10_1-'])
            print(eq)
            expr = []

            for i in range(0,len(eq)):
                expr.append(parse_expr(eq[i]))

            print(expr)

            symbs = list(expr[0].free_symbols)

            if str(symbs[0]) == 'x':
                symbx = symbs[0]
                symby = symbs[1]
            elif str(symbs[0]) == 'y':
                symbx = symbs[1]
                symby = symbs[0]
            else:
                symbx = symbs[0]
                symby = symbs[1]

            if len(expr) == 2:
                jac = jacobian(expr,[symbx,symby])
            else:
                jac = jacobian(expr,symbs)

            slog10_1 = ''

            slog10_1 = slog10_1 + "Sistema de " +str(len(expr))+"x"+str(len(expr))+'\n'





            for i in range(0,len(expr)):
                window['-LogNewtonNoLinear-'].update('f(x,y) ='+str(expr[i])+'\n', append=True)


            window['-LogNewtonNoLinear-'].update("Matriz Jacobiana"+'\n\n', append=True)

            window['-LogNewtonNoLinear-'].update(str(jac)+'\n\n', append=True)


            x_0L10_1 = list(map(float,values['-x0L10_1-'].split(',')))


            window['-LogNewtonNoLinear-'].update("Inicial Aprox."+'\n\n', append=True)

            window['-LogNewtonNoLinear-'].update(str(x_0L10_1)+'\n\n', append=True)

            solutionL10_1 = newton_method(expr,x_0L10_1,symbs)


            window['-LogNewtonNoLinear-'].update("La solucion es: "+'\n\n', append=True)

            window['-LogNewtonNoLinear-'].update(str(solutionL10_1[0])+'\n\n', append=True)



            if len(expr) == 2:
                p = plot_implicit(expr[0],(symbx,-10,10),(symby,-10,10),show=False)
                p.append((plot_implicit(expr[1],(symbx,-10,10),(symby,-10,10),show=False))[0])
                p2 = get_sympy_subplots(p)
                p2.plot(solutionL10_1[1],solutionL10_1[2],"o")
                p2.show()

            #window.Element('-LogNewtonNoLinear-').update(value=slog10_1)

        except Exception as e:
            print(e)
            sg.popup_ok("Algo Salio Mal intente otra vez :(")

#--------------------------------------------------
    if event == '-InterpolationButton-':
        window['-COL{0}-'].update(visible=False)
        window['-COL{11}-'].update(visible=True)

#--------------------------------------------------

    if event == '-LagrangeInterpol-':
        window['-COL{11}-'].update(visible=False)
        window['-COL{111}-'].update(visible=True)

    if event == '-SolveL11_1-':

        try:
            window.Element('-LogLagrangeInter-').update(value="")
            parseL11_1 = parse_interpolation(values['-DataInterpolationL11_1-'])
            print(parseL11_1)
            v11_1 = []
            fx11_1 = []

            for i in range(0,len(parseL11_1)):
                if parseL11_1[i][0] == 'x':
                    v11_1 = parse_array(parseL11_1[i][1])
                if parseL11_1[i][0] == 'f(x)':
                    fx11_1 = parse_array(parseL11_1[i][1])
            #print(v11_1)
            #print(fx11_1)
            countL11_1 = 0


            solutionL11_1 = Lagrange(v11_1, fx11_1)






        except Exception as e:
            print(e)
            sg.popup_ok("Algo Salio Mal intente otra vez :(")




#--------------------------------------------------
    if event == '-NewtonInterpol-':
        window['-COL{11}-'].update(visible=False)
        window['-COL{112}-'].update(visible=True)

    if event == '-SolveL11_2-':

        try:
            window.Element('-LogNewtonInter-').update(value="")
            parseL11_2 = parse_interpolation(values['-DataInterpolationL11_2-'])
            print(parseL11_2)
            v11_2 = []
            fx11_2 = []

            for i in range(0,len(parseL11_2)):
                if parseL11_2[i][0] == 'x':
                    v11_2 = parse_array(parseL11_2[i][1])
                if parseL11_2[i][0] == 'f(x)':
                    fx11_2 = parse_array(parseL11_2[i][1])
            print(v11_2)
            print(fx11_2)
            countL11_2 = 0


            solutionL11_2 = Newton_interpolation(fx11_2,v11_2 )

            window['-LogNewtonInter-'].update("El polinomio de interpolacion es: "+'\n\n', append=True)

            window['-LogNewtonInter-'].update(str(solutionL11_2)+'\n\n', append=True)





        except Exception as e:
            print(e)
            sg.popup_ok("Algo Salio Mal intente otra vez :(")


#--------------------------------------------------
    if event == '-HermiteInterpol-':
        window['-COL{11}-'].update(visible=False)
        window['-COL{113}-'].update(visible=True)

#--------------------------------------------------
    if event == '-SplineInterpol-':
        window['-COL{11}-'].update(visible=False)
        window['-COL{114}-'].update(visible=True)


    if event == '-SplineType-':
        if values['-SplineType-'] == 'Sujeto':
            window.Element('-DataInterpolationL11_4-').update(value="x = [0,1,2,3] \n f(x) = [0,0.5,2,1.5]\n f'(a) = 0.2\n f'(b) =-1")

        if values['-SplineType-'] == 'Natural':
            window.Element('-DataInterpolationL11_4-').update(value="x = [0,1,2,3] \n f(x) = [0,0.5,2,1.5]")

    if event == '-SolveL11_4-':

        window.Element('-LogSplineInter-').update(value="")
        if values['-SplineType-'] == 'Natural':
            try:
                parseL11_4 = parse_interpolation(values['-DataInterpolationL11_4-'])
                print(parseL11_4)
                v11_4 = []
                fx11_4 = []

                for i in range(0,len(parseL11_4)):
                    if parseL11_4[i][0] == 'x':
                        v11_4 = parse_array(parseL11_4[i][1])
                    if parseL11_4[i][0] == 'f(x)':
                        fx11_4 = parse_array(parseL11_4[i][1])
                print(v11_4)
                print(fx11_4)
                countL11_4 = 0


                solutionL11_4 = spline_natural(fx11_4,v11_4 )


                window['-LogSplineInter-'].update("Los splines son: "+'\n\n', append=True)

                window['-LogSplineInter-'].update("S(x): "+'\n\n', append=True)


                for i in range(0,len(solutionL11_4)):
                    window['-LogSplineInter-'].update(str(solutionL11_4[i])+'\n\n', append=True)

            except Exception as e:
                print(e)
                sg.popup_ok("Algo Salio Mal intente otra vez :(")
        else:

            try:
                parseL11_4 = parse_interpolation(values['-DataInterpolationL11_4-'])
                print(parseL11_4)
                v11_4 = []
                fx11_4 = []
                fpa = 0
                fpb = 0

                for i in range(0,len(parseL11_4)):
                    if parseL11_4[i][0] == 'x':
                        v11_4 = parse_array(parseL11_4[i][1])
                    if parseL11_4[i][0] == 'f(x)':
                        fx11_4 = parse_array(parseL11_4[i][1])
                    if parseL11_4[i][0] == "f'(a)":
                        fpa = float(parseL11_4[i][1])
                    if parseL11_4[i][0] == "f'(b)":
                        fpb = float(parseL11_4[i][1])

                print(v11_4)
                print(fx11_4)
                print(fpa)
                print(fpb)
                countL11_4 = 0


                solutionL11_4 = spline_sujeto(fx11_4,v11_4,fpa,fpb )


                window['-LogSplineInter-'].update("Los splines son: "+'\n\n', append=True)


                window['-LogSplineInter-'].update("S(x): "+'\n\n', append=True)


                for i in range(0,len(solutionL11_4)):
                    window['-LogSplineInter-'].update(str(solutionL11_4[i])+'\n\n', append=True)


            except Exception as e:
                print(e)
                sg.popup_ok("Algo Salio Mal intente otra vez :(")

#--------------------------------------------------

    if event == '-AproximationsButton-':
        window['-COL{0}-'].update(visible=False)
        window['-COL{12}-'].update(visible=True)

    if event == '-DisminQuadAprox-':
        window['-COL{12}-'].update(visible=False)
        window['-COL{121}-'].update(visible=True)

    if event == '-SolveL12_1-':
        window.Element('-LogDisMinQuad-').update(value="")
        lisimageL12_1 = []
        try:
            parseL12_1 = parse_interpolation(values['-DataAproximationL12_1-'])
            print(parseL12_1)

            v12_1 = []
            fx12_1 = []
            funcs = []

            for i in range(0,len(parseL12_1)):
                if parseL12_1[i][0] == 'x':
                    v12_1 = parse_array(parseL12_1[i][1])
                elif parseL12_1[i][0] == 'y':
                    fx12_1 = parse_array(parseL12_1[i][1])
                else:
                    funcs.append(parse_expr(parseL12_1[i][1]))

            print(v12_1)
            print(fx12_1)
            pprint(Matrix(funcs))

            countL12_1 = 0


            window['-LogDisMinQuad-'].update("X = "+'\n\n', append=True)


            window['-LogDisMinQuad-'].update(str(v12_1)+'\n\n', append=True)



            window['-LogDisMinQuad-'].update("Y = "+'\n\n', append=True)


            window['-LogDisMinQuad-'].update(str(fx12_1)+'\n\n', append=True)


            for i in range(0,len(funcs)):
                sr = "f"+str(i)+"(x) ="
                window['-LogDisMinQuad-'].update(sr+'\n\n', append=True)
                window['-LogDisMinQuad-'].update(str(funcs[i])+'\n\n', append=True)

            symb = ''

            for i in range(0,len(funcs)):
                if  len (funcs[i].free_symbols) > 0:
                    symb = list(funcs[i].free_symbols)
                    break

            #print(symb)
            solutionL12_1 = discrete_minimun_quads_aprox(v12_1,fx12_1,funcs,symb[0])


            window['-LogDisMinQuad-'].update("La solución esta dada por: "+'\n\n', append=True)

            window['-LogDisMinQuad-'].update(str(solutionL12_1)+'\n\n', append=True)

















        except Exception as e:

            print(e)
            sg.popup_ok("Algo Salio Mal intente otra vez :(")

    if event == '-ContminQuadAprox-':
        window['-COL{12}-'].update(visible=False)
        window['-COL{122}-'].update(visible=True)

    if event == '-SolveL12_2-':
        window.Element('-LogContMinQuad-').update(value="")
        lisimageL12_2 = []
        try:
            parseL12_2 = parse_interpolation(values['-DataAproximationL12_2-'])
            print(parseL12_2)

            v12_2 = []
            fx12_2 = []
            degr= 0

            for i in range(0,len(parseL12_2)):
                if parseL12_2[i][0] == 'f(x)':
                    fx12_2 = parse_expr(parseL12_2[i][1])
                elif parseL12_2[i][0] == 'I':
                    v12_2 = parse_array(parseL12_2[i][1])
                elif parseL12_2[i][0] == 'D':
                    degr = (float(parseL12_2[i][1]))

            pprint(Matrix(v12_2))
            pprint(fx12_2)
            print(degr)

            countL12_2 = 0


            window['-LogContMinQuad-'].update("f(x) = "+'\n\n', append=True)

            window['-LogContMinQuad-'].update(str(fx12_2)+'\n\n', append=True)



            window['-LogContMinQuad-'].update("[a,b] = "+'\n\n', append=True)

            window['-LogContMinQuad-'].update(str(v12_2)+'\n\n', append=True)




            symb = list(fx12_2.free_symbols)



            #print(symb)
            solutionL12_2 = continium_minimun_quads_aprox(fx12_2,v12_2,symb[0],int(degr))



            window['-LogContMinQuad-'].update("El polinomio de grado "+str(degr)+" esta dado por :"+'\n\n', append=True)

            window['-LogContMinQuad-'].update(str(solutionL12_2)+'\n\n', append=True)

















        except Exception as e:

            print(e)
            sg.popup_ok("Algo Salio Mal intente otra vez :(")

#--------------------------------------------------

    if event == '-returnL0-L10-':
        window['-COL{0}-'].update(visible=True)
        window['-COL{10}-'].update(visible=False)

    if event == '-returnL10-L10_1-':
        window['-COL{10}-'].update(visible=True)
        window['-COL{101}-'].update(visible=False)

#--------------------------------------------------

    if event == '-returnL0-L11-':
        window['-COL{0}-'].update(visible=True)
        window['-COL{11}-'].update(visible=False)


    if event == '-returnL11-L11_1-':
        window['-COL{11}-'].update(visible=True)
        window['-COL{111}-'].update(visible=False)

    if event == '-returnL11-L11_2-':
        window['-COL{11}-'].update(visible=True)
        window['-COL{112}-'].update(visible=False)

    if event == '-returnL11-L11_3-':
        window['-COL{11}-'].update(visible=True)
        window['-COL{113}-'].update(visible=False)

    if event == '-returnL11-L11_4-':
        window['-COL{11}-'].update(visible=True)
        window['-COL{114}-'].update(visible=False)

#--------------------------------------------------

    if event == '-returnL0-L12-':
        window['-COL{0}-'].update(visible=True)
        window['-COL{12}-'].update(visible=False)

    if event == '-returnL12-L12_1-':
        window['-COL{12}-'].update(visible=True)
        window['-COL{121}-'].update(visible=False)

    if event == '-returnL12-L12_2-':
        window['-COL{12}-'].update(visible=True)
        window['-COL{122}-'].update(visible=False)



    #print('You entered ', event)
    #print(values,"Valuessss")

window.close()

