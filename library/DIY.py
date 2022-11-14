import matplotlib.pyplot as plt
import numpy as np
from matplotlib.image import imread

# Defining function to show image
def show_image(array, title = None, dpi = 100, save = False, show = True):
    plt.figure(dpi = dpi, figsize = (0.3*array.shape[1]/dpi, 0.3*array.shape[0]/dpi))
    plt.axis(False)
    if title: plt.title(title)
    if show: plt.imshow(array, cmap='gray')
    if save!=False: plt.imsave(save, array, cmap='gray')

def comp_percent(original_shape, terms):
    original_size = original_shape[0]*original_shape[1]
    compressed_size = (sum(original_shape)+1)*terms
    return compressed_size/1024, (compressed_size)/original_size*100

def get_rms_error(A, A_):
    error = ((A-A_)**2).sum()
    return ((error/(A.shape[0]*A.shape[1]))**0.5)/2.55

def get_image(terms, U, S, V, A):
    A_ = U[:, :terms]@S[:terms, :terms]@V[:terms, :]
    A_ = A_.astype("int")

    size, comp = comp_percent(A.shape, terms)
    rms_error = get_rms_error(A, A_)

    text = f"{terms = }, size = {int(size)} KB ({round(comp, 1 if comp>1 else 2)}%)"
    return A_, comp, rms_error, text

def reduce(terms, A):
    U, S, V = np.linalg.svd(A)
    S = np.diag(S)
    return U[:, :terms]@S[:terms, :terms]@V[:terms, :]

def scale(A):
    u = A.max()
    l = A.min()
    if u == l or (u==255 and l==0): return A
    A = ((A-l)/(u-l)*255).astype(np.uint8)

class GrayscaleImageSVD:
    def __init__(self, location):
        self.location = location
        self.A = imread(location)
        if len(self.A.shape)!=2: raise ValueError("Image is not grayscale")
        self.U, S, self.V = np.linalg.svd(self.A)
        self.S = np.diag(S)
        
    def reduce(self, terms, type = np.uint8):
        return scale(self.U[:, :terms]@self.S[:terms, :terms]@self.V[:terms, :]).astype(type)
        

class ColourImageSVD:
    def __init__(self, location):
        self.location = location
        self.A = imread(location)
        if len(self.A.shape)!=3: raise ValueError("Image is not colour")
        self.R = self.A[:, :, 0]
        self.G = self.A[:, :, 1]
        self.B = self.A[:, :, 2]
        
        print("SVD for Red Channel")        
        U, S, V = np.linalg.svd(self.R)
        S = np.diag(S)
        self.RU, self.RS, self.RV = U, S, V
        
        print("SVD for Green Channel")
        U, S, V = np.linalg.svd(self.G)
        S = np.diag(S)
        self.GU, self.GS, self.GV = U, S, V
        
        print("SVD for Blue Channel")
        U, S, V = np.linalg.svd(self.B)
        S = np.diag(S)
        self.BU, self.BS, self.BV = U, S, V

    def reduce(self, terms, type = np.uint8):
        R = scale(self.RU[:, :terms]@self.RS[:terms, :terms]@self.RV[:terms, :])
        G = scale(self.GU[:, :terms]@self.GS[:terms, :terms]@self.GV[:terms, :])
        B = scale(self.BU[:, :terms]@self.BS[:terms, :terms]@self.BV[:terms, :])
        return np.dstack((R, G, B)).astype(type)
    