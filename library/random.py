def rand(seed, c = 3.5):
    """A function to generate a random number by using the formula:
        X[n+1] = c * X[n] * (1- X[n]).
        
    Args:
        seed (float): A value between 0 and 1. Seed is used if provided, else uses the last random value generated as seed.
        c (float, optional): Defaults to 3.5.
        
    Returns:
        float: a random number between 0 and 1.
    """
    return c*seed*(1-seed)

def LCG(seed, a = 1103515245, c = 12345, m = 32768):  # LCG
    """
        A function to generate a random number using the Linear Congruential Generator algorithm.
    Formula Used: X[n+1] = (a*X[n] + c) mod m.
    Value is divided by m before returning for scaling between 0 and 1.

    Args:
        seed (float, optional): A value between 0 and 1. Seed is used if provided, else uses the last random value generated as seed.
        a (int, optional): Defaults to 1103515245.
        c (int, optional): Defaults to 12345.
        m (int, optional): Defaults to 32768.

    Returns:
        float: a random number between 0 and 1.
    """
    return ((a * seed + c) % m) / m