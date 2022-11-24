try: from library.nonlinear_equations import solve_newton_raphson
except ImportError: from nonlinear_equations import solve_newton_raphson

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
        fn = lambda y: y0 + h * f(x0 + h, y) - y
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

def c_ode_2(f1:callable, f2:callable, x_f:float, x0:float, y0:float, z0:float, h:float=0.01):
    """Solves a coupled ODE system of two equations.

    Args:
        f1 (callable): First function to be integrated.
        f2 (callable): Second function to be integrated.
        x_f (float): Final x value.
        x0 (float): Initial x value.
        y0 (float): Initial y value.
        z0 (float): Initial z value.
        h (float, optional): Step size.. Defaults to 0.01.

    Returns:
        X, Y, Z: Approximate solutions to the ODE system.
    """

    xs, ys, zs = [x0], [y0], [z0]
    i = 0
    while x0<=x_f:
        k1y = h*f1(x0, y0, z0)
        k1z = h*f2(x0, y0, z0)

        k2y = h*f1(x0+0.5*h, y0+0.5*k1y, z0+0.5*k1z)
        k2z = h*f2(x0+0.5*h, y0+0.5*k1y, z0+0.5*k1z)

        k3y = h*f1(x0+0.5*h, y0+0.5*k2y, z0+0.5*k2z)
        k3z = h*f2(x0+0.5*h, y0+0.5*k2y, z0+0.5*k2z)

        k4y = h*f1(x0+h, y0+k3y, z0+k3z)
        k4z = h*f2(x0+h, y0+k3y, z0+k3z)
        
        x0 = x0+h
        y0 += (1/6)*(k1y+2*k2y+2*k3y+k4y)
        z0 += (1/6)*(k1z+2*k2z+2*k3z+k4z)
        xs.append(x0)
        ys.append(y0)
        zs.append(z0)
        
    return xs, ys, zs