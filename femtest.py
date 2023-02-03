import numpy as np

def nx2n(n_Rows, n_Columns):
    Zeros = []
    for i in range(n_Rows):
        Zeros.append([])
        for j in range(n_Columns*2):
            Zeros[i].append(0)
    return Zeros

# Applying matrix coefficients
def update(inputs, n_Rows, n_Columns, Zero):
    for i in range(n_Rows):
        for j in range(n_Columns):
            Zero[i][j] = inputs[i][j]
    return Zero

# Augmenting Identity Matrix of Order n
def identity(n_Rows, n_Columns, Matrix):
    for i in range(n_Rows):
        for j in range(n_Columns):
            if i == j:
                Matrix[i][j+n_Columns] = 1
    return Matrix

# Applying & implementing the GJE algorithm
def Gussain_Jordan_Elimination(n_Rows, n_Columns, Matrix):
    for i in range(n_Rows):
        if Matrix[i][i] == 0:
            print('error cannot divide by "0"')
    
        for j in range(n_Columns):
            if i != j:
                ratio = Matrix[j][i]/Matrix[i][i]

                for k in range(2*n_Columns):
                    Matrix[j][k] = Matrix[j][k] - ratio * Matrix[i][k]
    return Matrix

# Row Operation to make Principal Diagonal Element to '1'
def row_op(n_Rows, n_Columns, Matrix):
    for i in range(n_Rows):
        divide = Matrix[i][i]
        for j in range(2*n_Columns):
            Matrix[i][j] = Matrix[i][j]/divide
    return Matrix

# Display Inversed Matix
def Inverse(Matrix):
    returnable = []
    number_Rows = int(len(Matrix))
    number_Columns = int(len(Matrix[0]))
    Inversed_Matrix = (row_op(number_Rows, number_Columns, 
        Gussain_Jordan_Elimination(number_Rows, number_Columns, 
        identity(number_Rows, number_Columns, 
        update(Matrix, number_Rows, number_Columns, 
        nx2n(number_Rows, number_Columns))))))

    for i in range(number_Rows):
        returnable.append([])
        for j in range(number_Columns, 2*number_Columns):
            returnable[i].append(Inversed_Matrix[i][j])
    return returnable


A = np.array([[ 1166.67,     0.,   -1066.67,   133.33,     0.,    -333.33,  -100.,     200.  ],
 [    0.,     666.67,   200.,    -400. ,   -333.33,     0. ,    133.33,  -266.67],
 [-1066.67,   200. ,   1166.67,  -333.33 , -100.  ,   133.33  ,   0. ,      0.  ],
 [  133.33 , -400.  ,  -333.33 ,  666.67 ,  200. ,   -266.67  ,   0. ,      0.  ],
 [    0. ,   -333.33,  -100. ,    200. ,   1166.67 ,    0. ,  -1066.67,   133.33],
 [ -333.33 ,    0. ,    133.33 , -266.67 ,    0. ,    666.67 ,  200. ,   -400.  ],
 [ -100. ,    133.33 ,    0. ,      0. ,  -1066.67,   200. ,   1166.67 , -333.33],
 [  200. ,   -266.67  ,   0. ,      0.  ,   133.33 , -400. ,   -333.33 ,  666.67]])


print(Inverse(A))