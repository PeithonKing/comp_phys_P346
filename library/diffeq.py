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