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

    return list(x_np1-mm),

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

    prints = [["Aproximación Inicial",Matrix(x_0)]]
    xs = []
    ys = []

    while True:
        x_1 = NewtonMethod(ff,x_0,symbs)
        print(x_1[0])
        ninf = norm_inf(x_0,x_1[0])
        #print(ninf)
        prints.append(x_1[1])

        x_0 = list(x_1[0])
        xs.append(x_0[0])
        ys.append(x_0[1])
        if ninf < 1e-6:
            break

    print(prints)
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


'NewtonM3.png'
layout10_1 =  [

[sg.Button('<--',size=(2,2),key='-returnL10-L10_1-'),
sg.Text('Método de Newton Para Sistemas de Ecuaciones No Lineales')],

[sg.Text("Ingrese el sistema de Ecuaciones a Resolver\t\t"),sg.Text('Registro del proceso'),
 sg.Button('SOLVE',key='-SolveL10_1-',size=(5,5)),sg.Button('Graficar',key='-PlotL10_1-',size=(5,5))],
[sg.Text("Valor inicial"),sg.Input("X_0",key='-x0L10_1-',size=(10,10) ),
 sg.Text("plot range"),sg.Input("X_range",key='-x_rangeplotL10_1-',size=(10,10) ),
 sg.Input("Y_range",key='-y_rangeplotL10_1-',size=(10,10) ),sg.Input("Z_range",key='-z_rangeplotL10_1-' ,size=(10,10))],

[sg.Multiline("x**3+3*y**2-21 = 0\nx**2+2*y+2 = 0",key='-NoLinearSystemL10_1-',size=(30,30),horizontal_scroll=True),
 sg.Multiline("",key='-LogNewtonNoLinear-',size=(80,20),horizontal_scroll=True,auto_refresh=True)],

]






layout = [
    [sg.Column(layout=layout0,key='-COL{0}-',visible=True),

     sg.Column(layout=layout10,key='-COL{10}-',visible=False),
     sg.Column(layout=layout10_1,key='-COL{101}-',visible=False),
    ]
]






fix_expression("2x^2 + 3yx - y^2 = 0")

window = sg.Window('Métodos Numéricos ', layout,size=(720,480),resizable=True)

matop = '+'
#Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'Cancel': # if user closes window or clicks cancel
        break


    if event == '-NoLinearEqButton-':
        window['-COL{0}-'].update(visible=False)
        window['-COL{10}-'].update(visible=True)

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
                    xr_1 = -10
                    xr_2 = 10
                    yr_1 = -10
                    yr_2 = 10
                    zr_1 = -10
                    zr_2 = 10

                    if len(values['-y_rangeplotL10_1-']) > 0 and values['-y_rangeplotL10_1-'] != "Y_range" :
                        a1 = map(float,values['-y_rangeplotL10_1-'].split(','))
                        print(list(a1))
                        yr_1 = float(list(a1)[0])
                        yr_2 = float(list(a1)[1])
                        print("HI,",yr_1,yr_2)

                    if len(values['-x_rangeplotL10_1-']) > 0 and values['-x_rangeplotL10_1-'] != "X_range" :
                        a2 = map(float,values['-x_rangeplotL10_1-'].split(','))
                        print(list(a2))
                        xr_1 = min(list(a2))
                        xr_2 = max(list(a2))

                    if len(values['-z_rangeplotL10_1-']) > 0 and values['-z_rangeplotL10_1-'] != "Z_range" :
                        a3 = map(float,values['-z_rangeplotL10_1-'].split(','))
                        print(list(a3))
                        zr_1 = min(list(a3))
                        zr_2 = max(list(a3))

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

            jac = jacobian(expr,[symbx,symby])

            slog10_1 = ''

            slog10_1 = slog10_1 + "Sistema de " +str(len(expr))+"x"+str(len(expr))+'\n'




            preview(('f(x,y) =',expr[0],'\n'),output='png',viewer='file',filename='NewtonM1.png')
            preview(('f(x,y) = ',expr[1]),output='png',viewer='file',filename='NewtonM2.png')


            countL10_1 = countL10_1 + 1
            add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM1.png')
            window['-LogNewtonNoLinear-'].update('\n\n', append=True)

            countL10_1 = countL10_1+1
            add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM2.png')
            window['-LogNewtonNoLinear-'].update('\n\n', append=True)

            preview("Matriz Jacobiana",output='png',viewer='file',filename='NewtonM3.png')
            countL10_1 = countL10_1+1
            add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM3.png')
            window['-LogNewtonNoLinear-'].update('\n\n', append=True)

            preview(Matrix(jac),output='png',viewer='file',filename='NewtonM4.png')
            countL10_1 = countL10_1+1
            add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM4.png')
            window['-LogNewtonNoLinear-'].update('\n\n', append=True)


            x_0L10_1 = [-1,1]

            iterL10_1 = 0

            while True:
                preview("Iteración " +str(iterL10_1),output='png',viewer='file',filename='NewtonM1'+str(iterL10_1)+'.png')
                countL10_1 = countL10_1+1
                add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM1'+str(iterL10_1)+'.png')
                window['-LogNewtonNoLinear-'].update('\n\n', append=True)

                jev = Matrix(eval_matrix(jac,x_0L10_1,[symbx,symby]))

                preview("J(" + str(x_0L10_1) + ")",output='png',viewer='file',filename='NewtonM2'+str(iterL10_1)+'.png')
                countL10_1 = countL10_1+1
                add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM2'+str(iterL10_1)+'.png')
                window['-LogNewtonNoLinear-'].update('\n\n', append=True)

                preview(jev,output='png',viewer='file',filename='NewtonM3'+str(iterL10_1)+'.png')
                countL10_1 = countL10_1+1
                add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM3'+str(iterL10_1)+'.png')
                window['-LogNewtonNoLinear-'].update('\n\n', append=True)


                jinv = jev.inv()

                preview("Inv(J(" + str(x_0L10_1) + "))",output='png',viewer='file',filename='NewtonM4'+str(iterL10_1)+'.png')
                countL10_1 = countL10_1+1
                add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM4'+str(iterL10_1)+'.png')
                window['-LogNewtonNoLinear-'].update('\n\n', append=True)


                preview(jinv,output='png',viewer='file',filename='NewtonM5'+str(iterL10_1)+'.png')
                countL10_1 = countL10_1+1
                add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM5'+str(iterL10_1)+'.png')
                window['-LogNewtonNoLinear-'].update('\n\n', append=True)


                ffev = Matrix(evalVector(np.transpose(expr),x_0L10_1,[symbx,symby]))


                preview("F(" + str(x_0L10_1) + ")",output='png',viewer='file',filename='NewtonM6'+str(iterL10_1)+'.png')
                countL10_1 = countL10_1+1
                add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM6'+str(iterL10_1)+'.png')
                window['-LogNewtonNoLinear-'].update('\n\n', append=True)


                preview(ffev,output='png',viewer='file',filename='NewtonM7'+str(iterL10_1)+'.png')
                countL10_1 = countL10_1+1
                add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM7'+str(iterL10_1)+'.png')
                window['-LogNewtonNoLinear-'].update('\n\n', append=True)


                mm = jinv*ffev

                preview((jinv,'*',ffev),output='png',viewer='file',filename='NewtonM8'+str(iterL10_1)+'.png')
                countL10_1 = countL10_1+1
                add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM8'+str(iterL10_1)+'.png')
                window['-LogNewtonNoLinear-'].update('\n\n', append=True)


                preview(mm,output='png',viewer='file',filename='NewtonM9'+str(iterL10_1)+'.png')
                countL10_1 = countL10_1+1
                add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM9'+str(iterL10_1)+'.png')
                window['-LogNewtonNoLinear-'].update('\n\n', append=True)

                x_np1 = Matrix(np.array(x_0L10_1))


                preview(("X_n+1 = ",x_np1,'-',mm),output='png',viewer='file',filename='NewtonM8'+str(iterL10_1)+'.png')
                countL10_1 = countL10_1+1
                add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM8'+str(iterL10_1)+'.png')
                window['-LogNewtonNoLinear-'].update('\n\n', append=True)


                preview(x_np1-mm,output='png',viewer='file',filename='NewtonM9'+str(iterL10_1)+'.png')
                countL10_1 = countL10_1+1
                add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM9'+str(iterL10_1)+'.png')
                window['-LogNewtonNoLinear-'].update('\n\n', append=True)

                ninsf  = norm_inf(x_0L10_1,x_np1-mm)



                preview(ninsf,output='png',viewer='file',filename='NewtonM10'+str(iterL10_1)+'.png')
                countL10_1 = countL10_1+1
                add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM10'+str(iterL10_1)+'.png')
                window['-LogNewtonNoLinear-'].update('\n\n', append=True)


                x_0L10_1 = x_np1-mm




                if ninsf < 1e-6:


                    break

                iterL10_1 = iterL10_1+1

            solutionL10_1 = newton_method(expr,x_0L10_1,symbs)

            preview("La solucion es: " + str(solutionL10_1),output='png',viewer='file',filename='NewtonM12'+str(iterL10_1)+'.png')
            countL10_1 = countL10_1+1
            add_monkeys(window,countL10_1,'-LogNewtonNoLinear-','NewtonM12'+str(iterL10_1)+'.png')
            window['-LogNewtonNoLinear-'].update('\n\n', append=True)

            #window.Element('-LogNewtonNoLinear-').update(value=slog10_1)

        except Exception as e:
            print(e)
            sg.popup_ok("Algo Salio Mal intente otra vez :(")



    if event == '-returnL0-L10-':
        window['-COL{0}-'].update(visible=True)
        window['-COL{10}-'].update(visible=False)

    if event == '-returnL10-L10_1-':
        window['-COL{10}-'].update(visible=True)
        window['-COL{101}-'].update(visible=False)



    #print('You entered ', event)
    #print(values,"Valuessss")

window.close()
