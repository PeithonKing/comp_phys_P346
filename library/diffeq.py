try:
    from library.nonlinear_equations import solve_newton_raphson
    from library.myrandom import Random
    from library.matrix import Matrix, ones, zeros
except ImportError:
    from nonlinear_equations import solve_newton_raphson
    from myrandom import Random

import matplotlib.pyplot as plt
from tqdm import tqdm

def forward_euler(f, x, x0, y0, h=0.01):
    """Forward Euler method for solving ODEs.

    Args:
        f (function): Function to be integrated.
        x (int): x value to be evaluated.
        x0 (float): Initial x value.
        y0 (float): Initial y value.
        h (float): Step size.

    Returns:
        float: Approximate solution to the ODE.
    """
    xs, ys = [], []
    while x0 < x:
        y0 += h * f(x0, y0)
        x0 += h
        xs.append(x0)
        ys.append(y0)
    return xs, ys, y0


def backward_euler(f, x, x0, y0, h=0.01):
    """Backward Euler method for solving ODEs.

    Args:
        f (function): Function to be integrated.
        x (int): x value to be evaluated.
        x0 (float): Initial x value.
        y0 (float): Initial y value.
        h (float): Step size.

    Returns:
        float: Approximate solution to the ODE.
    """
    xs, ys = [], []
    while x0 < x:
        def fn(y): return y0 + h * f(x0 + h, y) - y
        ynp1 = solve_newton_raphson(fn)
        x0 += h
        y0 += h * f(x0, float(ynp1))
        xs.append(x0)
        ys.append(y0)
    return xs, ys, y0


def predictor_corrector(f, x, x0, y0, h=0.01):
    """Predictor-Corrector method for solving ODEs.

    Args:
        f (function): Function to be integrated.
        x (int): x value to be evaluated.
        x0 (float): Initial x value.
        y0 (float): Initial y value.
        h (float): Step size.

    Returns:
        float: Approximate solution to the ODE.
    """
    xs, ys = [], []
    while x0 < x:
        y1 = y0 + h * f(x0, y0)
        y2 = y0 + h * f(x0 + h, y1)
        y0 = y0 + h * (f(x0, y0) + f(x0 + h, y2)) / 2
        x0 += h
        xs.append(x0)
        ys.append(y0)
    return xs, ys, y0


def rk2(f, x, x0, y0, h=0.01):
    """Second-order Runge-Kutta method for solving ODEs.

    Args:
        f (function): Function to be integrated.
        x (int): x value to be evaluated.
        x0 (float): Initial x value.
        y0 (float): Initial y value.
        h (float): Step size.

    Returns:
        float: Approximate solution to the ODE.
    """
    xs, ys = [], []
    while x0 < x:
        k1 = h * f(x0, y0)
        k2 = h * f(x0 + h, y0 + k1)
        y0 += (k1 + k2) / 2
        x0 += h
        xs.append(x0)
        ys.append(y0)
    return xs, ys, y0


def rk4(f, x, x0, y0, h=0.01):
    """Fourth-order Runge-Kutta method for solving ODEs.

    Args:
        f (function): Function to be integrated.
        x (int): x value to be evaluated.
        x0 (float): Initial x value.
        y0 (float): Initial y value.
        h (float): Step size.

    Returns:
        float: Approximate solution to the ODE.
    """
    xs, ys = [], []
    while x0 < x:
        k1 = h * f(x0, y0)
        k2 = h * f(x0 + h / 2, y0 + k1 / 2)
        k3 = h * f(x0 + h / 2, y0 + k2 / 2)
        k4 = h * f(x0 + h, y0 + k3)
        y0 += (k1 + 2 * k2 + 2 * k3 + k4) / 6
        x0 += h
        xs.append(x0)
        ys.append(y0)
    return xs, ys, y0


def c_ode_2(
    f1: callable,
    f2: callable,
    x_get: float,
    x0: float,
    y0: float,
    z0: float,
    h: float = 0.01
):
    """Solves a coupled ODE system of two equations.

    Args:
        f1 (callable): First function to be integrated.
        f2 (callable): Second function to be integrated.
        x_get (float): Final x value.
        x0 (float): Initial x value.
        y0 (float): Initial y value.
        z0 (float): Initial z value.
        h (float, optional): Step size.. Defaults to 0.01.

    Returns:
        X, Y, Z: Approximate solutions to the ODE system.
    """

    xs, ys, zs = [x0], [y0], [z0]
    # i = 0
    while x0 <= x_get:
        k1y = h * f1(x0, y0, z0)
        k1z = h * f2(x0, y0, z0)

        k2y = h * f1(x0 + h/2, y0 + k1y/2, z0 + k1z/2)
        k2z = h * f2(x0 + h/2, y0 + k1y/2, z0 + k1z/2)

        k3y = h * f1(x0 + h/2, y0 + k2y/2, z0 + k2z/2)
        k3z = h * f2(x0 + h/2, y0 + k2y/2, z0 + k2z/2)

        k4y = h * f1(x0 + h, y0 + k3y, z0 + k3z)
        k4z = h * f2(x0 + h, y0 + k3y, z0 + k3z)

        x0 = x0+h
        y0 += (1/6)*(k1y+2*k2y+2*k3y+k4y)
        z0 += (1/6)*(k1z+2*k2z+2*k3z+k4z)
        xs.append(x0)
        ys.append(y0)
        zs.append(z0)

    return xs, ys, zs

def lag_interpol(zeta_h, zeta_l, yh, yl, y):
    return zeta_l + (zeta_h - zeta_l) * (y - yl)/(yh - yl)

def shooting_method(
    f1: callable,
    f2: callable,
    a: float,
    alpha: float,
    b: float,
    beta: float,
    h: float = 0.01,
    zeta_l: float = None,
    change: float = 1,
    seed:float =0.56,
    epsilon: float = 1e-6
):
    """Solves a coupled ODE system of two equations using the shooting method.

    Args:
        f1 (callable): First function to be integrated.
        f2 (callable): Second function to be integrated.
        a (float): Initial x value.
        alpha (float): Initial y value.
        b (float): Initial z value.
        beta (float): Initial z value.
        h (float, optional): Step Size. Defaults to 0.01.
        zeta_l (float, optional): Guess value of zeta. Defaults to None.
        change (float, optional): Amount to change. Defaults to 1.
        seed (float, optional): Seed if zeta_l not provided. Defaults to 0.56.
        epsilon (float, optional): Tollerance. Defaults to 1e-6.

    Returns:
        X, Y, Z: Approximate solutions to the ODE system.
    """
    if zeta_l is None:
        random = Random(seed)
        zeta_l = random.LCG()

    x, y, z = c_ode_2(f1, f2, b, a, alpha, zeta_l, h)
    yh = y[-1]
    if abs(yh - beta) < epsilon: return x, y, z
    sign0 = yh > beta  # True if yn > beta, False if yn < beta

    diff0 = abs(yh - beta)

    while True:
        zeta_h = zeta_l
        zeta_l -= change
        # print(f"{zeta_l = }")
        x, y, z = c_ode_2(f1, f2, b, a, alpha, zeta_l, h)
        yl = y[-1]
        if abs(yl - beta) < epsilon:
            return x, y, z
        if diff0 < abs(yl - beta): change = -change
        if (yl > beta) != sign0:
            break
        sign0 = yh > beta

    zeta_hope = lag_interpol(zeta_h, zeta_l, yh, yl, beta)

    x, y, z = c_ode_2(f1, f2, b, a, alpha, zeta_hope, h)
    yh = y[-1]
    # return x, y, z
    if abs(yh - beta) < epsilon:
        return x, y, z
    else: 
        return shooting_method(f1, f2, a, alpha, b, beta, h, zeta_hope, change/10, epsilon = epsilon)



def heat_eq2(temp:callable, Lx:float, Nx:int, Lt:float, Nt:int, needed:int):
    """Solves the heat equation in 1D.

    Args:
        temp (callable): Initial temperature distribution.
        Lx (float): Length of the rod.
        Nx (int): Number of grid points in x.
        Lt (float): Time to solve for.
        Nt (int): Number of time steps.
        needed (int): Upto the number of time steps to actually calculate.
    
    Returns:
        A: Approximate solution to the heat equation.
           Where each row is a time step, and each column
           is a point on the rod.
    """
    hx = Lx/Nx
    ht = Lt/Nt
    alpha = ht/(hx**2)
    print(f"{alpha=}")

    A = zeros((needed, Nx)).mat
    for i in range(Nx): A[0][i] = temp(Nx, i)
    for t in tqdm(range(1, needed)):
        for x in range(Nx):
            if x == 0:       A[t][x] = 0                 + (1 - 2*alpha)*A[t-1][x] + alpha*A[t-1][x+1]
            elif x == Nx-1:  A[t][x] = alpha*A[t-1][x-1] + (1 - 2*alpha)*A[t-1][x] + 0
            else:            A[t][x] = alpha*A[t-1][x-1] + (1 - 2*alpha)*A[t-1][x] + alpha*A[t-1][x+1]
    
    return A