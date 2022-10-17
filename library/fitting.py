try:
    from library.matrix import Matrix
    import library.matrix as m
    from library.linear_equations import Cholesky_Decomposition, forward_propagation, backward_propagation
except ImportError:
    from matrix import Matrix
    import matrix as m
    from linear_equations import Cholesky_Decomposition, forward_propagation, backward_propagation


def fit_data(x, y, n):
    A = m.ones([n, n], "A", 3)
    nem = [[len(x)]+[(x**j).sum() for j in range(1, n)]]
    A[0] = nem
    for i in range(1, n):
        nem[0].pop(0)
        nem[0].append((x**(i+n-1)).sum())
        A[i] = nem
    # coefficient matrix done
	

    Y = [[sum([y.mat[k][0]*x.mat[k][0]**i for k in range(len(x))])]
         for i in range(n)]
    Y = Matrix(Y, "Y", 3)
    # intercept matrix done
	

    L = Cholesky_Decomposition(A)
    y1 = forward_propagation(L, Y)
    a = backward_propagation(L.T(), y1)

    return a
