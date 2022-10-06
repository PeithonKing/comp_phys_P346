from copy import deepcopy as copy

try:
    from myrandom import Random
except:
    from library.myrandom import Random


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
        self.name = name
        self.precision = precision
        self.shape = [0, 0] if (mat == [] or mat == [[]]) else [
            len(mat), len(mat[0])]
        self.mat = mat if len(mat) else [[]]

    def __str__(self, precision: int = None):
        """Prints the matrix.

        Returns:
            string: Returns the matrix in string format.
        """
        if not precision:
            precision = self.precision
        ret = f"{self.name} = "
        t = " "*len(ret)
        ret += "|"
        for row in self.mat:
            for element in row:
                if element % 1 == 0:
                    element = int(element)
                j = round(element, precision)
                ret += f"{j}{' '*(8-len(str(j)))}"
            ret += f"\b\b\b\b\b|\n{t}|"
            # ret += f"{t}|"
        return ret[:-len(t)-1]

    def row(self, start: int = None, end: int = None, step: int = None):
        """Returns the ith row of the matrix.

        Args:
            start (int): the row start index.
            end (int): the row end index.
        """
        # print("ravan =", self.mat[start:end:step])
        return Matrix(self.mat[start:end:step], f"{self.name}[{start}:{end}:{step},:]", precision=self.precision)

    def col(self, start: int = None, end: int = None, step: int = None):
        """Returns the ith column of the matrix.

        Args:
            i (int): the column index.
        """
        return Matrix([row[start:end:step] for row in self.mat], f"{self.name}[:,{start}:{end}:{step}]", precision=self.precision)

    def swap_rows(self, i: int, j: int, verbose: bool = True):
        """swaps two rows of the matrix.

        Args:
            i, j (int, int): the two rows indices to be swapped.
            verbose (bool): Whether to print the message or not
        """
        if verbose:
            print(f"Swapping rows {i} and {j}")
        self.mat[i], self.mat[j] = self.mat[j], self.mat[i]

    def swap_cols(self, i: int, j: int, verbose: bool = True):
        """swaps two columns of the matrix.

        Args:
            i, j (int, int): the two column indices to be swapped.
            verbose (bool): Whether to print the message or not
        """
        if verbose:
            print(f"Swapping columns {i} and {j}")
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
        if verbose:
            print(f"R{subtract_from} - {round(scale, self.precision)}*R{subtract}")
        for k in range(self.shape[1]):
            self.mat[subtract_from][k] -= self.mat[subtract][k] * scale

    def divide(self, i: int, n: float, verbose: bool = True):
        """divides row i by n.

        Args:
            i (int): the row index to be multiplied.
            n (float): the scalar to be multiplied by row i.
            verbose (bool): Whether to print the message or not
        """
        if verbose:
            print(f"R{i} / {n}")
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
            print(self.mat)
            print(mat_2.mat)
            raise ValueError(
                f"Matrices cannot be multiplied with shape {self.shape} and {mat_2.shape}.")
        a = self.mat
        b = mat_2.mat

        res = [[sum(a[i][j]*b[j][k] for j in range(len(b)))
                for k in range(len(b[0]))] for i in range(len(a))]
        return Matrix(
            [[0]] if res == [[]] else res,
            name=f"{self.name}×{mat_2.name}",
            precision=4
        )

    def is_symmetric(self):
        """Checks if the matrix is symmetric.

        Returns:
            bool: True if the matrix is symmetric, False otherwise.
        """
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                if self.mat[i][j] != self.mat[j][i]:
                    return False
        return True

    def T(self):
        """Returns the transpose of the matrix.

        Returns:
            Matrix: The transpose of this matrix.
        """
        return Matrix(
            [[self.mat[i][j] for i in range(self.shape[0])]
             for j in range(self.shape[1])],
            name=f"{self.name}ᵀ",
            precision=self.precision
        )

    def dot(self, mat_2):
        """Dot product of this matrix and mat_2.

        Args:
            mat_2 (Matrix): A matrix of the same horizontal dimension as this matrix.

        Returns:
            Matrix: The dot product of this matrix and mat_2.
        """
        if self.shape[0] != mat_2.shape[0]:
            raise ValueError(
                "Matrices don't have the same first index: dot product cannot be cannot be found.")
        return (self @ mat_2.T()).mat[0][0]

    def det(self):
        """Returns the determinant of the matrix."""
        m = self.mat
        # base case for 2x2 matrix
        if len(m) == 2:
            return m[0][0]*m[1][1]-m[0][1]*m[1][0]

        determinant = 0
        for c in range(len(m)):
            determinant += ((-1)**c)*m[0][c]*self.minor(0, c).det()
        return determinant

    def minor(self, i, j):
        return Matrix([row[:j] + row[j+1:] for row in (self.mat[:i]+self.mat[i+1:])])

    def inverse(self):
        """Returns the inverse of the matrix.

        Returns:
            Matrix: The inverse of this matrix.
        """
        determinant = self.det()
        m = self.mat
        # special cases
        if len(m) == 1:
            return Matrix([[1/m[0][0]]], f"{self.name}⁻¹", self.precision)
        if len(m) == 2:
            return Matrix(self.T(), f"{self.name}⁻¹", self.precision)

        # find matrix of cofactors
        cofactors = []
        for r in range(len(m)):
            cofactorRow = []
            for c in range(len(m)):
                cofactorRow.append(((-1)**(r+c)) * self.minor(r, c).det())
            cofactors.append(cofactorRow)
        cofactors = Matrix(cofactors).T().mat
        for r in range(len(cofactors)):
            for c in range(len(cofactors)):
                cofactors[r][c] = cofactors[r][c]/determinant
        return Matrix(cofactors, f"{self.name}⁻¹", self.precision)

    # dunder methods:

    def __abs__(self): return self.det()

    def __add__(self, b):
        if isinstance(b, Matrix):
            if self.shape != b.shape:
                raise ValueError(
                    f"Matrices cannot be subtracted with shape {self.shape} and {b.shape}.")
            return Matrix([[self.mat[i][j]+b.mat[i][j] for j in range(self.shape[1])] for i in range(self.shape[0])], self.name, self.precision)
        elif isinstance(b, (int, float)):
            return Matrix([[self.mat[i][j]+b for j in range(self.shape[1])]
                           for i in range(self.shape[0])], self.name, self.precision)
        else:
            print(f"We don't take {type(b)} as an argument.")

    def __sub__(self, b):
        return self+(-b)

    def __mul__(self, other, name=None):
        if isinstance(other, Matrix):
            if other.shape[0] == 1 and other.shape[1] == self.shape[1]:  # row matrix
                return Matrix([[self.mat[i][j]*other.mat[0][j] for j in range(self.shape[1])] for i in range(self.shape[0])], self.name, self.precision)
            if other.shape[1] == 1 and other.shape[0] == self.shape[0]:  # column matrix
                return Matrix([[self.mat[i][j]*other.mat[i][0] for j in range(self.shape[1])] for i in range(self.shape[0])], self.name, self.precision)
        elif isinstance(other, (int, float)):
            if not name:
                name = f"{self.name}.{other}"
            m = copy(self.mat)
            for i in range(len(m)):
                for j in range(len(m[i])):
                    m[i][j] *= other
            return Matrix(m, name, self.precision)

    def __matmul__(self, other): return self.matmul(other)
    def __imul__(self, b): self = self * b
    def __truediv__(self, b): return self * (1/b)
    def __itruediv__(self, b): self = self * (1/b)
    def __len__(self): return self.shape[0]
    def __neg__(self): return self * -1

    def __pow__(self, b):
        if not b:
            return Matrix.identity(self.shape[0])
        elif b == 1:
            return self
        elif b == -1:
            return self.inverse()
        elif b > 1:
            return self * self**(b-1)
        elif b < -1:
            return self.inverse() * self**(-b-1)

    def __eq__(self, __o: object):
        if not isinstance(__o, Matrix):
            print(f"Only matrices can be equal to matrices")
            return False
        if self.shape != __o.shape:
            print(f"shapes are different {self.shape} and {__o.shape}")
            return False
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                if round(self.mat[i][j], self.precision) != round(__o.mat[i][j], self.precision):
                    print(f"elements are different self.mat[{i}][{j}]={self.mat[i][j]} and __o.mat[{i}][{j}]={__o.mat[i][j]}")
                    return False
        return True
            

    def __getitem__(self, s):
        if isinstance(s, tuple) and len(s) > len(self.shape):
            raise IndexError("Too many indices.")
        if isinstance(s, int):
            return Matrix([self.mat[s]], f"{self.name}[{s}]", self.precision)
        if isinstance(s, slice):
            start = s.start
            end = s.stop
            step = s.step
            a = self.row(start, end, step)
            a.name = f"{self.name}[{print_slice(s)}]"
            return a
        if isinstance(s, tuple):
            s1 = [s[1].start, s[1].stop, s[1].step] if isinstance(s[1], slice) else [
                s[1], s[1]+1, 1]
            # print(s1)
            a = self[s[0]].col(*s1)
            a.name = self.name + \
                "[" + print_slice(s[0]) + "," + print_slice(s[1]) + "]"
            return a

    def __setitem__(self, key, value):
        # print(key, type(key))
        if not isinstance(value, Matrix):
            if isinstance(value, (list, tuple)):
                value = Matrix(value)
                # print(value.shape)
            elif isinstance(value, (int, float)) and self[key].shape == [1, 1]:
                # print("b is an integer")
                if isinstance(key, int):
                    if self.shape[0] == 1:
                        key = (0, key)
                    elif self.shape[1] == 1:
                        key = (key, 0)
                    else:
                        raise IndexError("Too less indices.")
                self.mat[key[0]][key[1]] = value
                return
            else:
                raise TypeError(
                    "Value must be a matrix or a list/tuple of lists/tuples.")
        if self.__getitem__(key).shape != value.shape:
            raise ValueError("Cannot set matrix to matrix of different shape.")
        if isinstance(key, int):
            key = (slice(key, key+1, 1), slice(0, self.shape[1], 1))
        if isinstance(key, slice):
            key = (key, slice(0, self.shape[1], 1))
        if isinstance(key, tuple):
            key = [key[i] if isinstance(key[i], slice) else slice(
                key[i], key[i]+1, 1) for i in range(len(key))]
        # print(key, type(key))
        start_x = key[0].start if key[0].start != None else 0
        end_x = key[0].stop if key[0].stop != None else self.shape[0]
        step_x = key[0].step if key[0].step != None else 1
        start_y = key[1].start if key[1].start != None else 0
        end_y = key[1].stop if key[1].stop != None else self.shape[1]
        step_y = key[1].step if key[1].step != None else 1
        # print(start_x, end_x, step_x, start_y, end_y, step_y)

        i = 0
        for i_ in range(start_x, end_x, step_x):
            j = 0
            for j_ in range(start_y, end_y, step_y):
                # print(i_, j_, i, j)
                self.mat[i_][j_] = value.mat[i][j]
                j += 1
            i += 1


def print_slice(s):
    if isinstance(s, int):
        return str(s)
    else:
        _start_ = s.start if s.start != None else ""
        _end_ = s.stop if s.stop != None else ""
        _step_ = s.step if s.step != None else ""
        return f"{_start_}:{_end_}:{_step_}"


def zeros(shape, name=None, precision=3):
    return Matrix([[0 for j in range(shape[1])] for i in range(shape[0])], name, precision)


def ones(shape, name=None, precision=3):
    return Matrix([[1 for j in range(shape[1])] for i in range(shape[0])], name, precision)


def identity(n, name=None, precision=3):
    return Matrix([[1 if i == j else 0 for j in range(n)] for i in range(n)], name, precision)


def randmat(shape, seed=0.15, name=None, precision=3):
    r = Random(seed)
    return Matrix([[r.LCG() for j in range(shape[1])] for i in range(shape[0])], name, precision)


if __name__ == "__main__":
    A = Matrix(
        [[1, 2, 3],
         [4, 5, 6],
         [7, 8, 9]], "A", 3
    )
    print(A)
    A[1, 1:] = [[66, 67]]
    print(A)
