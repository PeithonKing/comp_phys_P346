try: from library.matrix import Matrix, randmat, truncate_p
except ImportError: from matrix import Matrix, randmat
import math


def get_eigen(A:Matrix, epsilon:float = 1e-6, seed:float = 0.1):
    precision = int(-math.log10(epsilon))
    x0 = randmat((3, 1), name="x0", precision=precision, seed=seed)
    
    z = x0
    i = 0
    lambda_old = 0
    while True:
        y = z
        z = A@y
        lambda_new = ((z.T()@x0)/(y.T()@x0)).mat[0][0]
        i += 1
        # print(f"{lambda_new = }")
        if abs(lambda_old-lambda_new)<=epsilon:
            break
        lambda_old = lambda_new
    
    y = y.normalise()
    
    y.name = "e-vector"
    return truncate_p(lambda_new, precision, str), y, i