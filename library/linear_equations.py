try:
    from matrix import Matrix, zeros, ones, randmat
except:
    from library.matrix import Matrix, zeros, ones, randmat


def gauss_jordan(A, B, verbose=False):
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
            m = max(enumerate(A.col(imp).mat[imp:]),
                    key=lambda x: abs(x[1][0]))[0]
            A.swap_rows(imp, imp+m, verbose)

        A.divide(imp, A.mat[imp][imp], verbose)

        for i in range(A.shape[0]):
            if imp != i:
                A.subtract(i, imp, A.mat[i][imp], verbose)

        if verbose:
            print()

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
        # if not L[i][i] == 1:
        #     print(f"L[{i}][{i}] = {L[i][i]} != 1")
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
        x[i] = [(y[i][0] - sum([U[i][j]*x[j][0]
                 for j in range(len(U)) if x[j]]))/U[i][i]]
    return Matrix(x, "x", p)


def LU_Decomposition(A):
    """LU decomposition

    Args:
        A (Matrix): Coefficient matrix

    Returns:
        Matrix, Matrix: L and U matrices
    """
    n = A.shape[0]
    LU = A

    for i in range(1, n):
        for j in range(i):
            if j == 0:
                LU[i, j] = LU[i, j] / LU[0, 0].mat[0][0]
            else:
                LU[i, j] = (LU[i, j] - LU[i, :j] @ LU[:j, j]) / LU[j, j].mat[0][0]
            # had to do index 0 separately brcause it didn't
            # have anything on it's left and I thought it would
            # have been a lot of work to generalise them to 0
        LU[i, i:] = LU[i, i:] - LU[i, :i] @ LU[:i, i:]
    
    # Was a single matrix till now, but while returning, I'm splitting it into L and U
    LU = LU.mat
    L = Matrix([[LU[i][j] if j<i else (1 if i==j else 0) for j in range(A.shape[1])] for i in range(A.shape[0])], "L", 3)
    U = Matrix([[LU[i][j] if j>=i else 0 for j in range(A.shape[1])] for i in range(A.shape[0])], "U", 3)

    return L, U


def gauss_siedel_solver(A, B, x):
    name = x.name
    precision = x.precision
    shape = A.shape
    A, B, x = A.mat, B.mat, x.mat
    x_new = []
    for i in range(len(A)):
        x_new.append([(B[i][0] - sum([A[i][j] * x[j][0] for j in range(shape[0]) if i != j])) / A[i][i]])
    return Matrix(x_new, name, precision)

def gauss_siedel(A, B, tol = 1e-6, seed = 0.1, max_iter = 100):
    x = randmat((A.shape[0], 1), seed,  "x")*5
    i = -1
    while True:
        i += 1
        # print(i, end = " ")
        x_old = x.mat
        x = gauss_siedel_solver(A, B, x)
        if sum([abs(x.mat[i][0]-x_old[i][0]) for i in range(len(x.mat))]) < tol or i > max_iter:
            break
    x.name = "x"
    return x, i

def Cholesky_Decomposition(A):
        """Returns the Cholesky decomposition of the matrix.

        Returns:
            Matrix: Lower Triangular matrix L such that A = LLᵀ. Diagonal elements of L will be changed to 1.
            Matrix: Upper Triangular matrix Lᵀ renamed as U. Here the diagonal elements are kept intact.
        """
        if A.is_symmetric():
            n = A.shape[0]
            matrix = A.mat

            l_triag = [[0 for x in range(n)]
                       for y in range(n)]

            for i in range(n):
                for j in range(i + 1):
                    sum1 = 0
                    if (j == i):
                        for k in range(j):
                            sum1 += l_triag[j][k]**2
                        l_triag[j][j] = (matrix[j][j] - sum1)**0.5
                    else:
                        for k in range(j):
                            sum1 += l_triag[i][k] * l_triag[j][k]
                        if (l_triag[j][j] > 0):
                            l_triag[i][j] = (
                                matrix[i][j] - sum1) / l_triag[j][j]
            ans2 = Matrix(l_triag, precision=3).T()
            ans2.name = "U"
            for i in range(len(l_triag)):
                l_triag[i][i] = 1

            ans1 = Matrix(l_triag, "L", 3)

            return ans1, ans2
        else:
            raise ValueError(
                "Matrix is not symmetric. Cannot perform Cholesky decomposition.")


if __name__ == "__main__":
    A = Matrix(
        [
            [1, -1, 4, 0, 2, 9],
            [0, 5, -2, 7, 8, 4],
            [1, 0, 5, 7, 3, -2],
            [6, -1, 2, 3, 0, 8],
            [-4, 2, 0, 5, -5, 3],
            [0, 7, -1, 5, 4, -2],
        ], "A", 3
    )
    L = LU_Decomposition(A)
    print(L)