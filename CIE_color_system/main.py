import colour as c
import numpy as np
from numpy.linalg import norm
from CIE_color_fast import *

d = np.random.random([3, 100])
d = d / norm(d, axis=0)
y1 = XYZ_to_sRGB(d)
test = True
for i in range(100):
    # y1 = my_XYZ_to_sRGB(d[i])
    y1i = y1[:, i]
    y2 = c.XYZ_to_sRGB(d[:, i])
    ratio = y2/y1i
    # print(y1i, ratio)
    if any(np.abs(ratio - 1) > 0.01):
        test = False
print("Tests passed?", test)
