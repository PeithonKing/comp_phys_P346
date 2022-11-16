import matplotlib.pyplot as plt
from library.DIY import *

pic = "crab"

print("creating image object")
img = ColourImageSVD(f"static/{pic}.jpg")

print("reducing image")
A = []
comp = []
errors = []
terms = [1, 5, 10, 20, 50, 100, 200, 500]
places = 2
for n_terms in terms:
    A2, ratio, error = img.reduce(terms = n_terms)
    comp.append(round(ratio, places))
    errors.append(round(error, places))
    A.append(A2)

fig = plt.figure(dpi = 500)
gs = fig.add_gridspec(2, 4, wspace=0.05, hspace=-0.55)
axs = gs.subplots(sharex=True, sharey=True)

for i in range(8):
    axs.flat[i].imshow(A[i], cmap='gray')
    axs.flat[i].axis("off")
    axs.flat[i].set_title(f"{terms[i]} terms, size = {comp[i]}%\nof original image", fontsize = 5)

plt.savefig(f"static/{pic}_evolution.png", bbox_inches='tight')
# plt.show()
