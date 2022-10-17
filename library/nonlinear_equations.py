import warnings
try:
    from myrandom import Random
except ImportError:
    from library.myrandom import Random


def differentiate(f, x, epsilon = 1e-6):  # Numerical differentiation
	return (f(x+epsilon)-f(x))/epsilon

def P(x, coeff):  # Polynomial
    return sum([c*x**(len(coeff)-i-1) for i, c in enumerate(coeff)])

def dP(coeff):  # derivative of polynomial
    coeff = coeff[::-1]
    dcoeff = [coeff[i] * i for i in range(len(coeff))]
    return dcoeff[-1:0:-1]

def get_brackets(f, guess = None):
    if not guess:
        r = Random(0.1, [-5, 5])
        guess = [r.LCG(), r.LCG()]
    a, b = guess
    # print(f"Finding brackets...{a}, {b}")
    if f(a)*f(b) < 0:
        return [a, b]
    alpha = 0.1
    if abs(f(a)) > abs(f(b)):
        b += alpha*(b-a)
    elif abs(f(b)) > abs(f(a)):
        a += alpha*(b-a)
    else:
        r = Random(a+b, [-5, 5])  # a+b is a random number
        guess = [r.LCG(), r.LCG()]
    return get_brackets(f, [a, b])

def solve_bisection(f, guess = None, epsilon = 1e-4, delta=1e-4):
    if not guess:
        a, b = get_brackets(f, guess)
    else:
        a, b = guess
    c = (a+b)/2
    if f(a)*f(c) > 0:
        a = c
    elif f(b)*f(c) > 0:
        b = c
    else:
        raise ValueError("No root in interval")
    # print(a, b)
    return solve_bisection(f, [a, b]) if (abs(a-b) > epsilon or abs(f(a))>delta or abs(f(b))>delta) else (a+b)/2

def solve_regula_falsi(f, guess = None, epsilon = 1e-2, delta=1e-2):
    if not guess:
        a, b = get_brackets(f, guess)
    else:
        a, b = guess
    c = a + (b-a)*abs(f(a))/(abs(f(a))+abs(f(b)))
    # print("\nc =", c)
    if f(a)*f(c) > 0:
        a = c
    else:
        b = c
    return solve_regula_falsi(f, [a, b]) if (abs(f(a))>delta or abs(f(b))>delta) else [a, b][abs(f(a))>abs(f(b))]

def solve_newton_raphson(f, f_d = None, guess = None, epsilon = 1e-4, delta=1e-4):
    if f_d==None:
        warnings.warn("No derivative provided, using numerical differentiation")
        f_d = lambda x: differentiate(f, x)
    if not guess:
        r = Random(0.1, [-5, 5])
        guess = r.LCG()
    guess = guess - f(guess)/f_d(guess)
    # print(guess)
    return solve_newton_raphson(f, f_d, guess) if (abs(f(guess))>delta) else guess

def laguerre(coeff, b = None, epsilon = 1e-6):
    if b == None:
        r = Random(0.1, [-5, 5])
        b = r.LCG()
    for n in range(1000):
        if abs(P(b, coeff))<epsilon:
            return b
        G = P(b, dP(coeff))/P(b, coeff)
        H = G**2 - P(b, dP(dP(coeff)))/P(b, coeff)
        n = len(coeff)
        D = ((n-1)*(n*H-G*G))**0.5
        
        if abs(G+D) > abs(G-D):
            a = n/(G+D)
        else:
            a = n/(G-D)

        if abs(a) < epsilon:
            return b
        b = b - a

def deflation(coeff, a):  # coeff = [-1, 3, 0, -4]
    l = [coeff[0]]
    i = 1
    for i in range(len(coeff)-1):
        l.append(coeff[i+1] + a*l[-1])
    return l[:-1]

def laguerre_solve(coeff, epsilon = 1e-6):
    roots = []
    for i in range(len(coeff)-1):
        roots.append(laguerre(coeff, b = 0))
        coeff = deflation(coeff, roots[-1])
    return roots