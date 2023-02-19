#!/usr/bin/python
import PySimpleGUI as sg
import sys
from sympy import *



sg.theme('DarkGreen')

#sg.theme_previewer()



# It's a class that represents a matrix and has functions that allow you to add, subtract, multiply, and find the inverse
# of a matrix
class Matrix:

    def __init__(self, matrix):
        """
        The function takes a matrix as an argument and assigns the number of rows and columns to the variables self.rows and
        self.columns
        :param matrix: the matrix we're going to be working with
        """
        self.matrix = matrix
        self.rows = len(matrix)
        self.columns = len(matrix[0])

    def __repr__(self):
        """
        The function __repr__() returns a string representation of the object
        :return: The string "Matrix()"
        """
        return "Matrix()"

    def __str__(self):
        """
        This function returns a string representation of the matrix
        :return: The matrix is being returned.
        """
        mat = ''
        for i in range(0,self.rows):
            for j in range(0,self.columns):
                mat = mat + str(self.matrix[i][j])+'\t '
            mat = mat + "\n"
        return mat

    def __add__(self,other):
        """
        It takes two matrices, checks if they have the same size, and if they do, it adds them together and returns the
        result
        :param other: The other matrix to add to this one
        :return: A new matrix with the sum of the two matrices
        """
        if self.rows != other.rows and  self.columns != other.columns:
            raise Exception("Matrices doesn't have the same size :(")

        m = []

        for i in range(0,self.rows):
            aux = []
            for j in range(self.columns):
                aux.append(self.matrix[i][j] + other.matrix[i][j])
            m.append(aux)

        return Matrix(m)

    def __sub__(self,other):
        """
        It subtracts two matrices.
        :param other: The other matrix to be subtracted from the current matrix
        :return: A new matrix with the subtraction of the two matrices.
        """
        if self.rows != other.rows and  self.columns != other.columns:
            raise Exception("Matrices doesn't have the same size :(")

        m = []

        for i in range(0,self.rows):
            aux = []
            for j in range(self.columns):
                aux.append(self.matrix[i][j] - other.matrix[i][j])
            m.append(aux)

        return Matrix(m)

    def __mul__(self,other):
        """
        For each row in the first matrix, multiply each element in that row by the corresponding element in the column of
        the second matrix, and sum the results
        :param other: The other matrix to multiply with
        :return: A matrix
        """
        if self.columns != other.rows:
            raise Exception("Matrices cannot be multiplied :(")
        m = []

        for i in range(0,self.rows):
            aux = []
            for j in range(other.columns):
                sum = 0
                for k in range(0,self.columns):
                    sum += self.matrix[i][k]*other[k][j]
                aux.append(sum)
            m.append(aux)

        return Matrix(m)

    def __getitem__(self, row):
        """
        The function takes in a row and column number and returns the value at that position in the matrix
        :param row: The row of the matrix you want to access
        :param column: The column number to return
        :return: The value of the matrix at the given row and column.
        """
        return self.matrix[row]

    def printeq_sys(self,v):
        """
        It takes a matrix and a vector and returns a string that is the matrix and vector printed in a way that is easy to
        read
        :param v: the vector of values
        :return: The matrix is being returned.
        """
        mat = ''
        for i in range(0,self.rows):
            for j in range(0,self.columns):
                mat = mat + str(self.matrix[i][j])+' '
            mat = mat + "|"+str(v[i])
            mat = mat + "\n"
        return mat

    def transpose(self):
        """
        It transposes the matrix.
        :return: A matrix object
        """
        mat = []
        for i in range(0,self.columns):
            aux = []
            for j in range(0,self.rows):
                aux.append(self.matrix[j][i])
            mat.append(aux)
        return Matrix(mat)

    def extract(self,row,column):
        """
        It returns a matrix with the row and column removed.
        :param row: The row to be removed
        :param column: The column to be removed
        :return: A matrix with the row and column removed.
        """
        mat = []

        for i in range(0,self.rows):
            if i == row:
                continue
            aux = []
            for j in range(0,self.columns):

                if j != column:
                    aux.append(self.matrix[i][j])

            mat.append(aux)

        return Matrix(mat)

    def det(self):
        """
        It takes the determinant of a matrix by recursively calling itself on the submatrices of the matrix until it reaches
        a 2x2 matrix, at which point it returns the determinant of the 2x2 matrix
        :return: The determinant of the matrix.
        """
        if self.rows == 2:
            return (self.matrix[0][0]*self.matrix[1][1])-(self.matrix[0][1]*self.matrix[1][0])
        dx = 0
        for i in range(0,self.rows):
            exmat = self.extract(0,i)
            #print(exmat)
            d = exmat.det()
            #print(d,'det')
            if i % 2 == 0:
                dx = dx + ((d*self.matrix[0][i]))
            else:
                dx = dx + ((d*self.matrix[0][i])*(-1))

            #print(i,"i,j,k,...",self.matrix[0][i], "-> " , dx)
        return dx

    def determinat(self):
        """
        It returns the determinant of a matrix.
        :return: The determinant of the matrix.
        """
        if self.rows != self.columns:
            raise Exception("the matrix not is square")
        return self.det()

    def adjoint(self):
        """
        The adjoint of a matrix is the transpose of the cofactor matrix
        :return: The adjoint of the matrix
        """
        if self.rows != self.columns:
            raise Exception("the matrix not is square")
        mat = []
        for i in range(0,self.rows):
            aux = []
            for j in range(0,self.columns):
                mc = self.extract(i,j)
                mcd = mc.determinat()
                if ((i+1) + (j+1)) % 2 != 0:
                    aux.append((-1)*mcd)
                else:
                    aux.append(mcd)

            mat.append(aux)

        return Matrix(mat)

    def inverse(self):
        """
        The function takes a matrix and returns its inverse
        :return: The inverse of the matrix
        """
        if self.rows != self.columns:
            raise Exception("the matrix not is square")


        adj = self.adjoint()
        #print(adj)
        adjt = adj.transpose()
        #print(adjt)
        det = self.det()
        #print(det)
        mat = []

        if det == 0:
            raise Exception("the matrix determinant is zero")
        det = 1/det
        #print(det)

        for i in range(0,self.rows):
            aux = []
            for j in range(0,self.columns):
                 aux.append(adjt[i][j]*det)
            mat.append(aux)

        return Matrix(mat)

    def move_row(self, from_row,to_row):
        """
        It takes the row from_row and moves it to row to_row
        :param from_row: the row you want to move
        :param to_row: The row you want to move to
        """
        a = self.matrix[from_row]
        b = self.matrix[to_row]
        self.matrix[from_row] = b
        self.matrix[to_row] = a



def gaussj_method(pivot, matrix,solutionv):
    """
    It takes a matrix, a pivot, and a solution vector, and it performs the Gauss-Jordan elimination method on the matrix,
    using the pivot as the pivot, and the solution vector as the solution vector

    :param pivot: the row we're currently working on
    :param matrix: the matrix of the system
    :param solutionv: the solution vector
    """

    maxm = matrix[pivot][pivot]
    row = pivot
    sreg = ""

    for i in range(pivot,matrix.rows):
        if (abs(maxm) < abs(matrix[i][pivot]) ):
            maxm = matrix[i][pivot]
            row = i

    matrix.move_row(row,pivot)
    aux = solutionv[row]
    solutionv[row] = solutionv[pivot]
    solutionv[pivot] = aux

    if row != pivot:
        sreg = sreg + ("R"+str(pivot+1)+"<----> R"+str(row+1)) +'\n'


    for i in range(0,matrix.rows):
        for j in range(0,matrix.columns):
            sreg = sreg + str(matrix[i][j]) + ' / '
        sreg = sreg + ' || '+str(solutionv[i])+'\n'
    #print(matrix.printeq_sys(solutionv))
    sreg = sreg + '\n'


    if matrix[pivot][pivot] != 1:
        a = matrix[pivot][pivot]
        for i in range(pivot,matrix.columns):
            #print(matrix[pivot][i])
            matrix[pivot][i] = matrix[pivot][i]/a
        solutionv[pivot] = solutionv[pivot]/a

        sreg = sreg + ("(1/"+ str(a) + ")R"+str(pivot+1)+"-----> R"+str(pivot+1))+'\n'

    for i in range(0,matrix.rows):
        for j in range(0,matrix.columns):
            sreg = sreg + str(matrix[i][j]) + ' / '
        sreg = sreg + ' || '+str(solutionv[i])+'\n'

    #print(matrix.printeq_sys(solutionv))

    for i in range(pivot+1,matrix.rows):
        pt = -1*matrix[i][pivot]
        sreg = sreg + (str(pt)+"R"+str(pivot+1)+"+ R"+str(i+1)+"------> R"+str(i+1))+'\n'
        for j in range(pivot,matrix.columns):
            #print(pt)
            matrix[i][j] = matrix[i][j]+((pt)*matrix[pivot][j])

        solutionv[i] = solutionv[i]+((pt)*solutionv[pivot])

    for i in range(0,matrix.rows):
        for j in range(0,matrix.columns):
            sreg = sreg + str(matrix[i][j]) + ' / '
        sreg = sreg + ' || '+str(solutionv[i])+'\n'

    #print(matrix.printeq_sys(solutionv))
    sreg = sreg +'\n\n'
    return sreg

def gauss_jordan(matrix,solutionv):
    """
    It takes a matrix and a solution vector and returns the solution vector

    :param matrix: The matrix that you want to solve
    :param solutionv: The solution vector
    :return: the solution vector x.
    """

    s = ''
    for i in range(0,matrix.rows):
        s = s + '||Pivot ' + str(i+1)+'||\n\n'
        s = s + gaussj_method(i,matrix,solutionv)

    x = []


    s = s + "||Sustitución regresiva||\n\n"
    for r in range(matrix.rows-1,-1,-1):
        if matrix[r][r] == 1:
            if r == matrix.rows-1:
                x.append(solutionv[r])
                s = s + "x_" + str(r+1) + " = " + str(solutionv[r]) + '\n'
            else:
                sum = 0
                s = s + "x_" + str(r+1) + " = (" + str(solutionv[r]) + "-( "
                for j in x:
                    for k in range(r+1,matrix.columns):
                        sum = sum+(j*matrix[r][k])
                        s = s + "( " + str(matrix[r][k])+ "*"+ str(j) + ") "
                        if k != matrix.columns-1:
                                s = s + '+ '



                s = s + ")"
                sum = sum*-1
                x.append((sum+solutionv[r])/matrix[r][r])
                s = s + "/ " + str(matrix[r][r]) + '\n'

    x.reverse()
    s = s + "||Las Soluciones del sistema son: ||\n\n"
    for i in range(0,matrix.rows):
        s = s + ("x_" + str(i+1)+" = "+  str(x[i])) + '\n'

    return x,s

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

    return Matrix(m)

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
    return Matrix(m)

def eval_matrix(matrix , v):
    e = 0
    mm = []
    for i in range(0,3):
        aux = []
        for j in range(0,3):
            aux.append(matrix[i][j].subs([(x,v[0]),(y,v[1]),(z,[2])]).evalf())
        mm.append(aux)
    return Matrix(mm)












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

#layout1 main menu for linear equations solving methods
layout1 = [
[sg.Button('<--',size=(2,2),key='returnL0-L1'),
sg.Text('METODOS PARA LA SOLUCION DE SISTEMAS DE ECUACIONES LINEALES')],

[sg.Button('Métodos Directos',size=(10,5),key='-LinearDirectButton-'),
sg.Text('          '),
sg.Button('Métodos Iterativos',size=(10,5),key='-LinearIterativeButton-'),
sg.Text('          '),
sg.Button('Métodos de Factorización',size=(10,5),key='-LinearFactorButton-')],

[sg.Text('')],

[sg.Button('Métodos Generales para Matrices',size=(10,5),key='-MatrixAppButton-')]
]
#---------------------------------------------------------------------
#Metods directos Menu
layout2 = [

[sg.Button('<--',size=(2,2),key='returnL1-L2'),
sg.Text('METODOS DIRECTOS PARA LA SOLUCION DE SISTEMAS DE ECUACIONES LINEALES')],

[sg.Button('Metodo de Gauss-Jordan',size=(10,5),key='-GaussMethod-'),
sg.Text('          '),
sg.Button('Método de Gauss-Jordan particionado ',size=(10,5),key='-GaussPartMethod-'),
sg.Text('          '),
sg.Button('Metodo de Matriz Inversa Particionado ',size=(10,5),key='-InversePartMethod-')
]
]

layout2_1 = [

    [sg.Button('<--',size=(2,2),key='returnL2-L2_1'),sg.Text('Metodo de Gauss-Jordan')],

[sg.Input('',key='-M00-',size=(4,4)), sg.Input('',key='-M01-',size=(4,4)),sg.Input('',key='-M02-',size=(4,4)),
sg.Input('',key='-M03-',size=(4,4)),sg.Input('',key='-M04-',size=(4,4)),sg.Input('',key='-M05-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B0-',size=(4,4))],

[sg.Input('',key='-M10-',size=(4,4)), sg.Input('',key='-M11-',size=(4,4)),sg.Input('',key='-M12-',size=(4,4)),
sg.Input('',key='-M13-',size=(4,4)),sg.Input('',key='-M14-',size=(4,4)),sg.Input('',key='-M15-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B1-',size=(4,4))],

[sg.Input('',key='-M20-',size=(4,4)), sg.Input('',key='-M21-',size=(4,4)),sg.Input('',key='-M22-',size=(4,4)),
sg.Input('',key='-M23-',size=(4,4)),sg.Input('',key='-M24-',size=(4,4)),sg.Input('',key='-M25-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B2-',size=(4,4))],

[sg.Input('',key='-M30-',size=(4,4)), sg.Input('',key='-M31-',size=(4,4)),sg.Input('',key='-M32-',size=(4,4)),
sg.Input('',key='-M33-',size=(4,4)),sg.Input('',key='-M34-',size=(4,4)),sg.Input('',key='-M35-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B3-',size=(4,4))],

[sg.Input('',key='-M40-',size=(4,4)), sg.Input('',key='-M41-',size=(4,4)),sg.Input('',key='-M42-',size=(4,4)),
sg.Input('',key='-M43-',size=(4,4)),sg.Input('',key='-M44-',size=(4,4)),sg.Input('',key='-M45-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B4-',size=(4,4))],

[sg.Input('',key='-M50-',size=(4,4)), sg.Input('',key='-M51-',size=(4,4)),sg.Input('',key='-M52-',size=(4,4)),
sg.Input('',key='-M53-',size=(4,4)),sg.Input('',key='-M54-',size=(4,4)),sg.Input('',key='-M55-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B5-',size=(4,4))],


[sg.Button('SOLVE',key='-SolveGaussJM-',size=(5,5))],
[sg.Text('Registro del proceso')],
[sg.Multiline("",key='-LogGaussM-',size=(100,100),horizontal_scroll=True)],


]

layout2_2 = [

    [sg.Button('<--',size=(2,2),key='returnL2-L2_2'),sg.Text('Metodo de Gauss-Jordan Particionado')],

[sg.Input('',key='-M00-',size=(4,4)), sg.Input('',key='-M01-',size=(4,4)),sg.Input('',key='-M02-',size=(4,4)),
sg.Input('',key='-M03-',size=(4,4)),sg.Input('',key='-M04-',size=(4,4)),sg.Input('',key='-M05-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B0-',size=(4,4))],

[sg.Input('',key='-M10-',size=(4,4)), sg.Input('',key='-M11-',size=(4,4)),sg.Input('',key='-M12-',size=(4,4)),
sg.Input('',key='-M13-',size=(4,4)),sg.Input('',key='-M14-',size=(4,4)),sg.Input('',key='-M15-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B1-',size=(4,4))],

[sg.Input('',key='-M20-',size=(4,4)), sg.Input('',key='-M21-',size=(4,4)),sg.Input('',key='-M22-',size=(4,4)),
sg.Input('',key='-M23-',size=(4,4)),sg.Input('',key='-M24-',size=(4,4)),sg.Input('',key='-M25-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B2-',size=(4,4))],

[sg.Input('',key='-M30-',size=(4,4)), sg.Input('',key='-M31-',size=(4,4)),sg.Input('',key='-M32-',size=(4,4)),
sg.Input('',key='-M33-',size=(4,4)),sg.Input('',key='-M34-',size=(4,4)),sg.Input('',key='-M35-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B3-',size=(4,4))],

[sg.Input('',key='-M40-',size=(4,4)), sg.Input('',key='-M41-',size=(4,4)),sg.Input('',key='-M42-',size=(4,4)),
sg.Input('',key='-M43-',size=(4,4)),sg.Input('',key='-M44-',size=(4,4)),sg.Input('',key='-M45-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B4-',size=(4,4))],

[sg.Input('',key='-M50-',size=(4,4)), sg.Input('',key='-M51-',size=(4,4)),sg.Input('',key='-M52-',size=(4,4)),
sg.Input('',key='-M53-',size=(4,4)),sg.Input('',key='-M54-',size=(4,4)),sg.Input('',key='-M55-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B5-',size=(4,4))],


[sg.Button('SOLVE',key='-SolveGaussJPartM-',size=(5,5))],
[sg.Text('Registro del proceso')],
[sg.Multiline("",key='-LogGaussPartM-',size=(100,100),horizontal_scroll=True)],


]

layout2_3 = [

    [sg.Button('<--',size=(2,2),key='returnL2-L2_3'),sg.Text('Metodo de Matriz Inversa Particionado')],

[sg.Input('',key='-M00-',size=(4,4)), sg.Input('',key='-M01-',size=(4,4)),sg.Input('',key='-M02-',size=(4,4)),
sg.Input('',key='-M03-',size=(4,4)),sg.Input('',key='-M04-',size=(4,4)),sg.Input('',key='-M05-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B0-',size=(4,4))],

[sg.Input('',key='-M10-',size=(4,4)), sg.Input('',key='-M11-',size=(4,4)),sg.Input('',key='-M12-',size=(4,4)),
sg.Input('',key='-M13-',size=(4,4)),sg.Input('',key='-M14-',size=(4,4)),sg.Input('',key='-M15-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B1-',size=(4,4))],

[sg.Input('',key='-M20-',size=(4,4)), sg.Input('',key='-M21-',size=(4,4)),sg.Input('',key='-M22-',size=(4,4)),
sg.Input('',key='-M23-',size=(4,4)),sg.Input('',key='-M24-',size=(4,4)),sg.Input('',key='-M25-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B2-',size=(4,4))],

[sg.Input('',key='-M30-',size=(4,4)), sg.Input('',key='-M31-',size=(4,4)),sg.Input('',key='-M32-',size=(4,4)),
sg.Input('',key='-M33-',size=(4,4)),sg.Input('',key='-M34-',size=(4,4)),sg.Input('',key='-M35-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B3-',size=(4,4))],

[sg.Input('',key='-M40-',size=(4,4)), sg.Input('',key='-M41-',size=(4,4)),sg.Input('',key='-M42-',size=(4,4)),
sg.Input('',key='-M43-',size=(4,4)),sg.Input('',key='-M44-',size=(4,4)),sg.Input('',key='-M45-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B4-',size=(4,4))],

[sg.Input('',key='-M50-',size=(4,4)), sg.Input('',key='-M51-',size=(4,4)),sg.Input('',key='-M52-',size=(4,4)),
sg.Input('',key='-M53-',size=(4,4)),sg.Input('',key='-M54-',size=(4,4)),sg.Input('',key='-M55-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B5-',size=(4,4))],


[sg.Button('SOLVE',key='-SolveInversePartM-',size=(5,5))],
[sg.Text('Registro del proceso')],
[sg.Multiline("",key='-LogInversePartM-',size=(100,100),horizontal_scroll=True)],


]

#---------------------------------------------------------------------
#Metodos Iterativos Menu
layout3 = [

[sg.Button('<--',size=(2,2),key='returnL1-L3'),
sg.Text('METODOS ITERATIVOS PARA LA SOLUCION DE SISTEMAS DE ECUACIONES LINEALES')],

[sg.Button('Metodo de Jacobi',size=(10,5),key='-JacobiMethod-'),
sg.Text('          '),
sg.Button('Metodo de Gauss-Seidel ',size=(10,5),key='-GaussSeidelMethod-'),
sg.Text('          ')
]
]

layout3_1 = [

    [sg.Button('<--',size=(2,2),key='returnL3-L3_1'),sg.Text('Metodo de Jacobi')],

[sg.Input('',key='-M00-',size=(4,4)), sg.Input('',key='-M01-',size=(4,4)),sg.Input('',key='-M02-',size=(4,4)),
sg.Input('',key='-M03-',size=(4,4)),sg.Input('',key='-M04-',size=(4,4)),sg.Input('',key='-M05-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B0-',size=(4,4))],

[sg.Input('',key='-M10-',size=(4,4)), sg.Input('',key='-M11-',size=(4,4)),sg.Input('',key='-M12-',size=(4,4)),
sg.Input('',key='-M13-',size=(4,4)),sg.Input('',key='-M14-',size=(4,4)),sg.Input('',key='-M15-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B1-',size=(4,4))],

[sg.Input('',key='-M20-',size=(4,4)), sg.Input('',key='-M21-',size=(4,4)),sg.Input('',key='-M22-',size=(4,4)),
sg.Input('',key='-M23-',size=(4,4)),sg.Input('',key='-M24-',size=(4,4)),sg.Input('',key='-M25-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B2-',size=(4,4))],

[sg.Input('',key='-M30-',size=(4,4)), sg.Input('',key='-M31-',size=(4,4)),sg.Input('',key='-M32-',size=(4,4)),
sg.Input('',key='-M33-',size=(4,4)),sg.Input('',key='-M34-',size=(4,4)),sg.Input('',key='-M35-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B3-',size=(4,4))],

[sg.Input('',key='-M40-',size=(4,4)), sg.Input('',key='-M41-',size=(4,4)),sg.Input('',key='-M42-',size=(4,4)),
sg.Input('',key='-M43-',size=(4,4)),sg.Input('',key='-M44-',size=(4,4)),sg.Input('',key='-M45-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B4-',size=(4,4))],

[sg.Input('',key='-M50-',size=(4,4)), sg.Input('',key='-M51-',size=(4,4)),sg.Input('',key='-M52-',size=(4,4)),
sg.Input('',key='-M53-',size=(4,4)),sg.Input('',key='-M54-',size=(4,4)),sg.Input('',key='-M55-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B5-',size=(4,4))],


[sg.Button('SOLVE',key='-SolveJacobiM-',size=(5,5))],
[sg.Text('Registro del proceso')],
[sg.Multiline("",key='-LogJacobiM-',size=(100,100),horizontal_scroll=True)],


]

layout3_2 = [

    [sg.Button('<--',size=(2,2),key='returnL3-L3_2'),sg.Text('Metodo de Gauss-Seidel')],

[sg.Input('',key='-M00-',size=(4,4)), sg.Input('',key='-M01-',size=(4,4)),sg.Input('',key='-M02-',size=(4,4)),
sg.Input('',key='-M03-',size=(4,4)),sg.Input('',key='-M04-',size=(4,4)),sg.Input('',key='-M05-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B0-',size=(4,4))],

[sg.Input('',key='-M10-',size=(4,4)), sg.Input('',key='-M11-',size=(4,4)),sg.Input('',key='-M12-',size=(4,4)),
sg.Input('',key='-M13-',size=(4,4)),sg.Input('',key='-M14-',size=(4,4)),sg.Input('',key='-M15-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B1-',size=(4,4))],

[sg.Input('',key='-M20-',size=(4,4)), sg.Input('',key='-M21-',size=(4,4)),sg.Input('',key='-M22-',size=(4,4)),
sg.Input('',key='-M23-',size=(4,4)),sg.Input('',key='-M24-',size=(4,4)),sg.Input('',key='-M25-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B2-',size=(4,4))],

[sg.Input('',key='-M30-',size=(4,4)), sg.Input('',key='-M31-',size=(4,4)),sg.Input('',key='-M32-',size=(4,4)),
sg.Input('',key='-M33-',size=(4,4)),sg.Input('',key='-M34-',size=(4,4)),sg.Input('',key='-M35-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B3-',size=(4,4))],

[sg.Input('',key='-M40-',size=(4,4)), sg.Input('',key='-M41-',size=(4,4)),sg.Input('',key='-M42-',size=(4,4)),
sg.Input('',key='-M43-',size=(4,4)),sg.Input('',key='-M44-',size=(4,4)),sg.Input('',key='-M45-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B4-',size=(4,4))],

[sg.Input('',key='-M50-',size=(4,4)), sg.Input('',key='-M51-',size=(4,4)),sg.Input('',key='-M52-',size=(4,4)),
sg.Input('',key='-M53-',size=(4,4)),sg.Input('',key='-M54-',size=(4,4)),sg.Input('',key='-M55-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B5-',size=(4,4))],


[sg.Button('SOLVE',key='-SolveGaussSeidelM-',size=(5,5))],
[sg.Text('Registro del proceso')],
[sg.Multiline("",key='-LogGaussSeidelM-',size=(100,100),horizontal_scroll=True)],


]


#---------------------------------------------------------------------
#Metodos de factorizacion Menu
layout4 = [
[sg.Button('<--',size=(2,2),key='returnL1-L4'),
sg.Text('METODOS DE FACTORIZACION PARA LA SOLUCION DE SISTEMAS DE ECUACIONES LINEALES')],

[sg.Button('Método de Factorización de Doolittle',size=(10,5),key='-DoolittleMethod-'),
sg.Text('          '),
sg.Button('Método de Factorización de Crout',size=(10,5),key='-CroutMethod-'),
sg.Text('          '),
sg.Button('Método de Factorización de Cholesky',size=(10,5),key='-CholeskyMethod-')],


]

layout4_1 = [

    [sg.Button('<--',size=(2,2),key='returnL4-L4_1'),sg.Text('Método de Factorización de Doolittle')],

[sg.Input('',key='-M00-',size=(4,4)), sg.Input('',key='-M01-',size=(4,4)),sg.Input('',key='-M02-',size=(4,4)),
sg.Input('',key='-M03-',size=(4,4)),sg.Input('',key='-M04-',size=(4,4)),sg.Input('',key='-M05-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B0-',size=(4,4))],

[sg.Input('',key='-M10-',size=(4,4)), sg.Input('',key='-M11-',size=(4,4)),sg.Input('',key='-M12-',size=(4,4)),
sg.Input('',key='-M13-',size=(4,4)),sg.Input('',key='-M14-',size=(4,4)),sg.Input('',key='-M15-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B1-',size=(4,4))],

[sg.Input('',key='-M20-',size=(4,4)), sg.Input('',key='-M21-',size=(4,4)),sg.Input('',key='-M22-',size=(4,4)),
sg.Input('',key='-M23-',size=(4,4)),sg.Input('',key='-M24-',size=(4,4)),sg.Input('',key='-M25-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B2-',size=(4,4))],

[sg.Input('',key='-M30-',size=(4,4)), sg.Input('',key='-M31-',size=(4,4)),sg.Input('',key='-M32-',size=(4,4)),
sg.Input('',key='-M33-',size=(4,4)),sg.Input('',key='-M34-',size=(4,4)),sg.Input('',key='-M35-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B3-',size=(4,4))],

[sg.Input('',key='-M40-',size=(4,4)), sg.Input('',key='-M41-',size=(4,4)),sg.Input('',key='-M42-',size=(4,4)),
sg.Input('',key='-M43-',size=(4,4)),sg.Input('',key='-M44-',size=(4,4)),sg.Input('',key='-M45-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B4-',size=(4,4))],

[sg.Input('',key='-M50-',size=(4,4)), sg.Input('',key='-M51-',size=(4,4)),sg.Input('',key='-M52-',size=(4,4)),
sg.Input('',key='-M53-',size=(4,4)),sg.Input('',key='-M54-',size=(4,4)),sg.Input('',key='-M55-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B5-',size=(4,4))],


[sg.Button('SOLVE',key='-SolveDoolittleM-',size=(5,5))],
[sg.Text('Registro del proceso')],
[sg.Multiline("",key='-LogDoolittleM-',size=(100,100),horizontal_scroll=True)],


]

layout4_2 = [

    [sg.Button('<--',size=(2,2),key='returnL4-L4_2'),sg.Text('Método de Factorización de Crout')],

[sg.Input('',key='-M00-',size=(4,4)), sg.Input('',key='-M01-',size=(4,4)),sg.Input('',key='-M02-',size=(4,4)),
sg.Input('',key='-M03-',size=(4,4)),sg.Input('',key='-M04-',size=(4,4)),sg.Input('',key='-M05-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B0-',size=(4,4))],

[sg.Input('',key='-M10-',size=(4,4)), sg.Input('',key='-M11-',size=(4,4)),sg.Input('',key='-M12-',size=(4,4)),
sg.Input('',key='-M13-',size=(4,4)),sg.Input('',key='-M14-',size=(4,4)),sg.Input('',key='-M15-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B1-',size=(4,4))],

[sg.Input('',key='-M20-',size=(4,4)), sg.Input('',key='-M21-',size=(4,4)),sg.Input('',key='-M22-',size=(4,4)),
sg.Input('',key='-M23-',size=(4,4)),sg.Input('',key='-M24-',size=(4,4)),sg.Input('',key='-M25-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B2-',size=(4,4))],

[sg.Input('',key='-M30-',size=(4,4)), sg.Input('',key='-M31-',size=(4,4)),sg.Input('',key='-M32-',size=(4,4)),
sg.Input('',key='-M33-',size=(4,4)),sg.Input('',key='-M34-',size=(4,4)),sg.Input('',key='-M35-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B3-',size=(4,4))],

[sg.Input('',key='-M40-',size=(4,4)), sg.Input('',key='-M41-',size=(4,4)),sg.Input('',key='-M42-',size=(4,4)),
sg.Input('',key='-M43-',size=(4,4)),sg.Input('',key='-M44-',size=(4,4)),sg.Input('',key='-M45-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B4-',size=(4,4))],

[sg.Input('',key='-M50-',size=(4,4)), sg.Input('',key='-M51-',size=(4,4)),sg.Input('',key='-M52-',size=(4,4)),
sg.Input('',key='-M53-',size=(4,4)),sg.Input('',key='-M54-',size=(4,4)),sg.Input('',key='-M55-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B5-',size=(4,4))],


[sg.Button('SOLVE',key='-SolveCroutM-',size=(5,5))],
[sg.Text('Registro del proceso')],
[sg.Multiline("",key='-LogCroutM-',size=(100,100),horizontal_scroll=True)],


]

layout4_3 = [

    [sg.Button('<--',size=(2,2),key='returnL4-L4_3'),sg.Text('Método de Factorización de Cholesky')],

[sg.Input('',key='-M00-',size=(4,4)), sg.Input('',key='-M01-',size=(4,4)),sg.Input('',key='-M02-',size=(4,4)),
sg.Input('',key='-M03-',size=(4,4)),sg.Input('',key='-M04-',size=(4,4)),sg.Input('',key='-M05-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B0-',size=(4,4))],

[sg.Input('',key='-M10-',size=(4,4)), sg.Input('',key='-M11-',size=(4,4)),sg.Input('',key='-M12-',size=(4,4)),
sg.Input('',key='-M13-',size=(4,4)),sg.Input('',key='-M14-',size=(4,4)),sg.Input('',key='-M15-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B1-',size=(4,4))],

[sg.Input('',key='-M20-',size=(4,4)), sg.Input('',key='-M21-',size=(4,4)),sg.Input('',key='-M22-',size=(4,4)),
sg.Input('',key='-M23-',size=(4,4)),sg.Input('',key='-M24-',size=(4,4)),sg.Input('',key='-M25-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B2-',size=(4,4))],

[sg.Input('',key='-M30-',size=(4,4)), sg.Input('',key='-M31-',size=(4,4)),sg.Input('',key='-M32-',size=(4,4)),
sg.Input('',key='-M33-',size=(4,4)),sg.Input('',key='-M34-',size=(4,4)),sg.Input('',key='-M35-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B3-',size=(4,4))],

[sg.Input('',key='-M40-',size=(4,4)), sg.Input('',key='-M41-',size=(4,4)),sg.Input('',key='-M42-',size=(4,4)),
sg.Input('',key='-M43-',size=(4,4)),sg.Input('',key='-M44-',size=(4,4)),sg.Input('',key='-M45-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B4-',size=(4,4))],

[sg.Input('',key='-M50-',size=(4,4)), sg.Input('',key='-M51-',size=(4,4)),sg.Input('',key='-M52-',size=(4,4)),
sg.Input('',key='-M53-',size=(4,4)),sg.Input('',key='-M54-',size=(4,4)),sg.Input('',key='-M55-',size=(4,4)),
sg.Text('|'),sg.Input('',key='-B5-',size=(4,4))],


[sg.Button('SOLVE',key='-SolveCholeskyM-',size=(5,5))],
[sg.Text('Registro del proceso')],
[sg.Multiline("",key='-LogCholeskyM-',size=(100,100),horizontal_scroll=True)],


]


#---------------------------------------------------------------------
layout5 = [
[sg.Button('<--',size=(2,2),key='returnL1-L5'),
sg.Text('METODOS PARA MATRICES')],

[sg.Text('FILAS ') ,sg.Text('        ' ),sg.Text('COLUMNAS')],

[sg.Combo([i for i in range(2,11)],key='-RowsSize-'),sg.Text('        ' ),
    sg.Combo([i for i in range(2,11)],key='-ColsSize-')],

[sg.Button('Operaciones Basicas Para Matrices',size=(10,5),key='-BasicMatrixMethods-'),
sg.Text('          '),
sg.Button('Matriz Inversa',size=(10,5),key='-InverseMethod-'),
sg.Text('          '),
sg.Button('Matriz Adjunta',size=(10,5),key='-AdjointMethod-')],

[sg.Text('')],

[sg.Button('Matriz Traspuesta',size=(10,5),key='-TransposeMethod-'),
sg.Text('          '),
sg.Button('Determinant',size=(10,5),key='-Determinant method-')]
]


layoutM1 = [
    [sg.Input('-',key=str(i)+ str(j)+'M1' ,size=(4,4))for j in range(0,5)] for i in range(0,5)
]

layoutM2 = [
    [sg.Input('-',key=str(i)+ str(j) +'M2',size=(4,4))for j in range(0,5)] for i in range(0,5)
]
#Main layout
layout = [
    [sg.Column(layout=layout0,key='-COL{0}-',visible=True),


    sg.Column(layout=layout1,key='-COL{1}-',visible=False),


    sg.Column(layout=layout2,key='-COL{2}-',visible=False),
    sg.Column(layout=layout2_1,key='-COL{21}-',visible=False),
    sg.Column(layout=layout2_2,key='-COL{22}-',visible=False),
    sg.Column(layout=layout2_3,key='-COL{23}-',visible=False),




    sg.Column(layout=layout3,key='-COL{3}-',visible=False),
    sg.Column(layout=layout3_1,key='-COL{31}-',visible=False),
    sg.Column(layout=layout3_2,key='-COL{32}-',visible=False),



    sg.Column(layout=layout4,key='-COL{4}-',visible=False),
    sg.Column(layout=layout4_1,key='-COL{41}-',visible=False),
    sg.Column(layout=layout4_2,key='-COL{42}-',visible=False),
    sg.Column(layout=layout4_3,key='-COL{43}-',visible=False),



    sg.Column(layout=layout5,key='-COL{5}-',visible=False),
    ]
]




#m = [[1,2,3,4],[2,3,4,6],[2,6,5,6],[2,325,6,3]]
m2 = [[5,2,0],[2,1,-1],[2,3,-2]]
b = [1,2,3]
#mt = Matrix(m)
mt2 = Matrix(m2)
#mt2inv = mt2.inverse()

#print(mt2)
#gauss_jordan(mt2,[1,2,3])
#print(mt2)
x,y,z = symbols('x,y,z')
symb = [x,y,z]
expr1 = x**2 + 2*y -z
expr2 = x**2 + 2*y
expr3 = x**2*(2*y**2)
ff = [expr1,expr2,expr3]
print(ff)







j = jacobian(ff,symb)
je = eval_matrix(j,[0,1,2])
print(j.inverse())




# Create the Window
window = sg.Window('Métodos Numéricos ', layout,size=(720,480),resizable=True)
#Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'Cancel': # if user closes window or clicks cancel
        break


    if event == '-LinearEqButton-':
        window['-COL{0}-'].update(visible=False)
        window['-COL{1}-'].update(visible=True)
#----------------------------------------------------------------
    if event == '-LinearDirectButton-':
        window['-COL{1}-'].update(visible=False)
        window['-COL{2}-'].update(visible=True)
#....................
    if event == '-GaussMethod-':
        window['-COL{2}-'].update(visible=False)
        window['-COL{21}-'].update(visible=True)

    if event == '-SolveGaussJM-':

        mm = []
        bb = []

        for i in range(0,6):
            aux = []
            for j in range(0,6):
                if values['-M'+str(i)+str(j)+'-'] != '':
                    try:
                        aux.append(float(values['-M'+str(i)+str(j)+'-']))
                    except:
                        continue
                else:
                    continue
            if len(aux) != 0:
                mm.append(aux)
        for k in range(0,6):
            if values['-B'+str(k)+'-'] != '':
                try:
                    bb.append(float(values['-B'+str(k)+'-']))
                except:
                    continue

        mat = Matrix(mm)
        print(mat.printeq_sys(bb))

        try:
            mm[0]
            len(mm)
            flag = True
        except:
            flag = False


        if flag and (len(mm) == len(mm[0])) and (len(bb) == len(mm)):

            if mat.determinat() != 0:
                g = gauss_jordan(mat,bb)
                window.Element('-LogGaussM-').update(value=g[1])

                for i in range(0,mat.rows):
                    for j in range(0,mat.columns):
                        if i == j:
                            window.Element('-M'+str(i)+str(j)+'-').update(value='1')
                        else:
                            window.Element('-M'+str(i)+str(j)+'-').update(value='0')
                    window.Element('-B'+str(i)+'-').update(value=str(g[0][i]))



            else:
                sg.popup_ok("El sistema No tiene soluciones")
        else:
            sg.popup_ok("La matriz no es cuadrado")

#....................
    if event == '-GaussPartMethod-':
        window['-COL{2}-'].update(visible=False)
        window['-COL{22}-'].update(visible=True)

#....................

    if event == '-InversePartMethod-':
        window['-COL{2}-'].update(visible=False)
        window['-COL{23}-'].update(visible=True)





#----------------------------------------------------------------


    if event == '-LinearIterativeButton-':
        window['-COL{1}-'].update(visible=False)
        window['-COL{3}-'].update(visible=True)

    if event == '-LinearFactorButton-':
        window['-COL{1}-'].update(visible=False)
        window['-COL{4}-'].update(visible=True)


    if event == '-MatrixAppButton-':
        window['-COL{1}-'].update(visible=False)
        window['-COL{5}-'].update(visible=True)







    if event == 'returnL0-L1':
        window['-COL{0}-'].update(visible=True)
        window['-COL{1}-'].update(visible=False)
#----------------------------------------------------------------
    if event == 'returnL1-L2':
        window['-COL{1}-'].update(visible=True)
        window['-COL{2}-'].update(visible=False)

    if event == 'returnL2-L2_1':
        window['-COL{2}-'].update(visible=True)
        window['-COL{21}-'].update(visible=False)
        for i in range(0,6):
            for j in range(0,6):
                window.Element('-M'+str(i)+str(j)+'-').update(value='')
            window.Element('-B'+str(i)+'-').update(value='')

    if event == 'returnL2-L2_2':
        window['-COL{2}-'].update(visible=True)
        window['-COL{22}-'].update(visible=False)
        for i in range(0,6):
            for j in range(0,6):
                window.Element('-M'+str(i)+str(j)+'-').update(value='')
            window.Element('-B'+str(i)+'-').update(value='')

    if event == 'returnL2-L2_3':
        window['-COL{2}-'].update(visible=True)
        window['-COL{23}-'].update(visible=False)
        for i in range(0,6):
            for j in range(0,6):
                window.Element('-M'+str(i)+str(j)+'-').update(value='')
            window.Element('-B'+str(i)+'-').update(value='')



#----------------------------------------------------------------
    if event == 'returnL1-L3':
        window['-COL{1}-'].update(visible=True)
        window['-COL{3}-'].update(visible=False)


    if event == 'returnL1-L4':
        window['-COL{1}-'].update(visible=True)
        window['-COL{4}-'].update(visible=False)

    if event == 'returnL1-L5':
        window['-COL{1}-'].update(visible=True)
        window['-COL{5}-'].update(visible=False)



    print('You entered ', event)
    #print(values,"Valuessss")

window.close()
