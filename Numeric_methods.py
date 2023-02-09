#!/usr/bin/python
import PySimpleGUI as sg
import sys




sg.theme('DarkGreen')

#sg.theme_previewer()



# The Matrix class takes in a matrix and allows you to add, subtract, and multiply matrices
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
                mat = mat + str(self.matrix[i][j])+' '
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
            if i & 1 == 0:
                dx = dx + ((d*self.matrix[0][i])*-1)
            else:
                dx = dx + (d*self.matrix[0][i])

            #print(i,"i,j,k,...",self.matrix[0][i])
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
                if ((i+1) + (j+1)) & 1:
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
        print(adj)
        adjt = adj.transpose()
        print(adjt)
        det = adjt.determinat()
        print(det)
        mat = []

        if det == 0:
            raise Exception("the matrix determinant is zero")
        det = 1/det
        print(det)

        for i in range(0,self.rows):
            aux = []
            for j in range(0,self.columns):
                 aux.append(adjt[i][j]*det)
            mat.append(aux)

        return Matrix(mat)











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


layout1 = [
[sg.Text('METODOS PARA LA SOLUCION DE SISTEMAS DE ECUACIONES NO LINEALES')],

[sg.Button('Métodos Directos',size=(10,5),key='-LinearDirectButton-'),
sg.Text('          '),
sg.Button('Métodos Iterativos',size=(10,5),key='-LinearIterativeButton-'),
sg.Text('          '),
sg.Button('Métodos de Factorización',size=(10,5),key='-LinearFactorButton-')],

[sg.Text('')],

[sg.Button('Métodos Generales para Matrices',size=(10,5),key='-AproximationsButton-')]
]




#Main layout
layout = [
    [sg.Column(layout=layout0,key='-COL{0}-',visible=True),
    sg.Column(layout=layout1,key='-COL{1}-',visible=False),
    ]
]




m = [[1,2,3,4],[2,3,5,6],[2,6,5,6],[2,325,6,3]]
m2 = [[2,0,1],[3,0,0],[5,1,1]]
mt = Matrix(m)
mt2 = Matrix(m2)
mt2inv = mt2.inverse()
print(mt2)
print(mt2inv)
# Create the Window
#window = sg.Window('Métodos Numéricos ', layout,size=(720,480))
# Event Loop to process "events" and get the "values" of the inputs
#while True:
#    event, values = window.read()
#    if event == sg.WIN_CLOSED or event == 'Cancel': # if user closes window or clicks cancel
#        break
#
#
#    if event == '-LinearEqButton-':
#        window['-COL{0}-'].update(visible=False)
#        window['-COL{1}-'].update(visible=True)
#
#    print('You entered ', event)
#
#window.close()
#
