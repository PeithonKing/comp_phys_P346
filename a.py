# from library.random import LCG
# import math

# def func1(x):
#     return x - 2*math.cos(x)

# def func2(x):
#     return math.cos(x) - x**3

# def func3(x):
#     return 3*x + math.sin(x) - math.e**x

# def func4(x):
#     return x*math.e**x - 2

# def get_brackets():
#     return [-1.5, 1.5]

# def solve(f, guess = None, epsilon = 1e-2, delta=1e-5):
#     if not guess:
#         a, b = get_brackets()
#     else:
#         a, b = guess
#     c = (a+b)/2
#     if f(a)*f(c) > 0:
#         a = c
#     elif f(b)*f(c) > 0:
#         b = c
#     else:
#         raise ValueError("No root in interval")
#     # print(a, b)
#     return solve(f, [a, b]) if (abs(a-b) > epsilon and abs(f(a))>delta and abs(f(b))>delta) else (a+b)/2

# print(solve(func2))


def laguerre_solve_equation(f, )