def sum_of_AP(n, a, d):
    """Calculate the sum of the Arithmetic Progression upto N terms.

    Args:
        n (int): the number of terms in the AP
        a (float): the first term in the AP.
        d (float): the difference between the terms in the AP.

    Returns:
        float: Sum of the AP upto N terms.
    """
    return n * (2 * a + (n - 1) * d) / 2

def factorial(x):
    """Calculate the factorial of x.

    Args:
        x (int): the number whose factorial is to be calculated.

    Returns:
        int: Factorial of x.
    """
    if x == 0:
        return 1
    else:
        return x * factorial(x - 1)

def distance(a, b):
    """Calculate the distance between two points on N dimensionlist.

    Args:
        a (list): the first point.
        b (list): the second point. length of a and b must be same.

    Raises:
        ValueError: if the length of a and b are not same.

    Returns:
        float: the distance between the two points.
    """
    if len(a) != len(b):
        raise ValueError("a and b must be of same length")
    s = [(a[i]-b[i])**2 for i in range(len(a))]
    return sum(s)**0.5

def RMS(x):
    """Calculate the Root Mean Square of the list x.

    Args:
        x (list): the list whose RMS is to be calculated.

    Returns:
        float: the RMS of the list x.
    """
    return (sum([i**2 for i in x])/len(x))**0.5

class Matrix:
    """A Class for Matrix"""
    def __init__(
        self,
        mat: list,
        name: str = None,
        precision: int = 16
        ):
        """Initialize the matrix

        Args:
            mat (2D list): The matrix
            name (string, optional): The name of the matrix; used while printing.
            precision (int, optional): The number of decimal places to be printed. Defaults to 16.
        """
        self.mat = mat
        self.name = name
        self.precision = precision
        self.shape = [len(mat), len(mat[0])]

    def __str__(self):
        """Prints the matrix.

        Returns:
            string: Returns the matrix in string format.
        """
        ret = f"{self.name} = "
        t = " "*len(ret)
        ret += "|"
        for row in self.mat:
            for element in row:
                if element%1 == 0:
                    element = int(element)
                ret += f"{round(element, self.precision)}\t"
            ret += f"|\n{t}|"
            # ret += f"{t}|"
        return ret[:-len(t)-1]
    
    def row(self, i: int):
        """Returns the ith row of the matrix.

        Args:
            i (int): the row index.
        """
        return self.mat[i]
    
    def col(self, i: int):
        """Returns the ith column of the matrix.

        Args:
            i (int): the column index.
        """
        return Matrix([[round(row[i], self.precision)] for row in self.mat], f"{self.name}[:,{i}]", precision = self.precision)

    def swap_rows(self, i: int, j: int, verbose: bool = True):
        """swaps two rows of the matrix.

        Args:
            i, j (int, int): the two rows indices to be swapped.
            verbose (bool): Whether to print the message or not
        """
        if verbose: print(f"Swapping rows {i} and {j}")
        self.mat[i], self.mat[j] = self.mat[j], self.mat[i]

    def swap_cols(self, i: int, j: int, verbose: bool = True):
        """swaps two columns of the matrix.

        Args:
            i, j (int, int): the two column indices to be swapped.
            verbose (bool): Whether to print the message or not
        """
        if verbose: print(f"Swapping columns {i} and {j}")
        for k in range(self.shape[0]):
            self.mat[k][i], self.mat[k][j] = self.mat[k][j], self.mat[k][i]

    def augment(self, mat_2):
        """appends the provided matrix to the right of this matrix.

        Args:
            mat_2 (Matrix): A matrix of the same vertical dimension as this matrix.
        """
        for i in range(self.shape[0]):
            self.mat[i] += mat_2.mat[i]
        self.shape = [len(self.mat), len(self.mat[0])]
    
    def subtract(self, subtract_from: int, subtract: int, scale: float, verbose: bool = True):
        """scales up "subtract" by "scale" and subtracts it from "subtract_from".

        Args:
            subtract_from, subtract (int, int): the two row indices to be subtracted.
            scale (float): the scalar to be multiplied by row "subtract".
            verbose (bool): Whether to print the message or not
        """
        if verbose: print(f"R{subtract_from} - {round(scale, self.precision)}*R{subtract}")
        for k in range(self.shape[1]):
            self.mat[subtract_from][k] -= self.mat[subtract][k] * scale

    def divide(self, i: int, n: float, verbose: bool = True):
        """multiplies row i by n.

        Args:
            i (int): the row index to be multiplied.
            n (float): the scalar to be multiplied by row i.
            verbose (bool): Whether to print the message or not
        """
        if verbose: print(f"R{i} / {n}")
        for k in range(self.shape[1]):
            self.mat[i][k] /= n

    def matmul(self, mat_2):
        """multiplies this matrix by mat_2.

        Args:
            mat_2 (Matrix): A matrix of the same horizontal dimension as this matrix.
        
        Returns:
            matrix: The product of this matrix and mat_2.
        """
        if self.shape[1] != mat_2.shape[0]:
            raise ValueError("Matrices cannot be multiplied.")
        a = self.mat
        b = mat_2.mat
        return Matrix(
            [[sum(a[i][j]*b[j][k] for j in range(len(b))) for k in range(len(b[0]))] for i in range(len(a))],
            name = f"{self.name}×{mat_2.name}",
            precision = 4
            )
    
    def dot(self, mat_2):
        """Dot product of this matrix and mat_2.

        Args:
            mat_2 (Matrix): A matrix of the same horizontal dimension as this matrix.

        Returns:
            Matrix: The dot product of this matrix and mat_2.
        """
        if self.shape[0] != mat_2.shape[0]:
            raise ValueError("Matrices cannot be multiplied.")
        return sum([self.mat[i][0] * mat_2.mat[i][0] for i in range(len(self.mat))])
    
    # def det(self):
        

class Complex:
    """A Class for Complex Numbers"""
    def __init__(
        self,
        real: float = 0,
        imag: float = 0,
        name: str = None,
        precision: int = 16
    ):
        """A complex number.

        Args:
            real (float): The real part of the complex number. Defaults to 0.
            imag (float): The imaginary part of the complex number. Defaults to 0.
            name (str, optional): The name of the Complex number; used while printing. Defaults to "complex".
            precision (int, optional): Number of digits to be printed after the decimal point. Defaults to 16.
        """ 
        self.real = real
        self.imag = imag
        self.name = name
        self.precision = precision
    
    def add(self, complex_2):
        """A function to return the sum of this complex and complex_2.
        
        Args:
            complex_2 (Complex): The complex number to be added to this complex number.
            
        Returns:
            Complex: The sum of this complex and complex_2. """ 
        return Complex(
            self.real + complex_2.real,
            self.imag + complex_2.imag,
            name = f"{self.name}+{complex_2.name}",
            precision = self.precision
        )
    
    def multiply(self, complex_2):
        """A function to return the product of this complex and complex_2

        Args:
            complex_2 (Complex): A complex number.

        Returns:
            Complex: A complex number representing the product of this complex and complex_2.
        """
        return Complex(
            self.real * complex_2.real - self.imag * complex_2.imag,
            self.real * complex_2.imag + self.imag * complex_2.real,
            name = f"{self.name}*{complex_2.name}",
            precision = self.precision
        )
    
    def mod(self):
        """mod value of this complex number.

        Returns:
            float: mod value of this complex number.
        """
        return (self.real ** 2 + self.imag ** 2)**0.5
    
    def __str__(self, name = True):
        """A function to return a string representation of this complex number."""
        a = ""
        if name and self.name:
            a +=  f"{self.name} = "
        if self.real:
            a += f"{round(self.real, self.precision)}"
            if self.imag:
                a += " + "
        if self.imag:
            a += f"{round(self.imag, self.precision)}i"
        if a == "":
            return "0 + 0i"
        if a == f"{self.name} = ":
            return f"{self.name} = 0 + 0i"
        return a