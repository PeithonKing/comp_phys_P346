from library.basic_arithmatic import Matrix

def gauss_jordan(A, B, verbose = False):
    """Perform Gauss-Jordan elimination on the augmented matrix [A|B]

    Args:
        A (Matrix): Coefficient matrix
        B (Matrix): Intercept matrix
        verbose (bool, optional): Whether to print the steps. Defaults to False.

    Returns:
        Matrix: The solution column matrix.
    """
    A.augment(B)
    for imp in range(A.shape[0]):
        if verbose:
            print(f"working with row {imp}")
        if A.mat[imp][imp] == 0:
            m = max(enumerate(A.col(imp).mat[imp:]), key=lambda x: abs(x[1][0]))[0]
            A.swap_rows(imp, imp+m, verbose)

        A.divide(imp, A.mat[imp][imp], verbose)

        for i in range(A.shape[0]):
            if imp != i:
                A.subtract(i, imp, A.mat[i][imp], verbose)
        
        if verbose: print()

    ans = A.col(-1)
    ans.name = "Sol"
    return ans

def forward_propagation(L, B):
    """forward propagation for cholesky decomposition

    Args:
        L (Matrix): Lower triangular matrix
        B (Matrix): Coefficient matrix

    Returns:
        Matrix: y matrix
    """
    p = L.precision
    L = L.mat
    B = B.mat
    y = []
    for i in range(len(L)):
        y.append([B[i][0] - sum([L[i][j]*y[j][0] for j in range(i)])])
    return Matrix(y, "y", p)

def backward_propagation(U, y):
    """backward propagation for cholesky decomposition

    Args:
        U (Matrix): Upper triangular matrix
        y (Matrix): y matrix

    Returns:
        Matrix: x matrix
    """
    p = U.precision
    U = U.mat
    y = y.mat
    x = [None for i in range(len(U))]
    for i in range(len(U)-1, -1, -1):
        x[i] = [(y[i][0] - sum([U[i][j]*x[j][0] for j in range(len(U)) if x[j]]))/U[i][i]]
    return Matrix(x, "x", p)

def L_U_Decomposition(A):
    """LU decomposition

    Args:
        A (Matrix): Coefficient matrix

    Returns:
        Matrix, Matrix: L and U matrices
    """
    p = A.precision
    mat = A.mat
    n = len(mat)
 
    lower = [[0 for x in range(n)]
             for y in range(n)]
    upper = [[0 for x in range(n)]
             for y in range(n)]
    for i in range(n):
        for k in range(i, n):
            sum = 0
            for j in range(i):
                sum += (lower[i][j] * upper[j][k])
            upper[i][k] = mat[i][k] - sum
        for k in range(i, n):
            if (i == k):
                lower[i][i] = 1
            else:
                sum = 0
                for j in range(i):
                    sum += (lower[k][j] * upper[j][i])
                lower[k][i] = (mat[k][i] - sum) / upper[i][i]
    
    return Matrix(lower, "L", p), Matrix(upper, "U", p)