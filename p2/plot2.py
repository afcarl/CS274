import numpy as np
import matplotlib.pyplot as plt
from pylab import *
# Generate data...
x=[-5, -3, -4, -6, -5, -3, 3, 4, 5]
y=[-3, -4, -5, 3,  5,  6,  7, 5, 5]

color = [1,1,1,2,2,2,3,3,3]
# Plot...
plt.scatter(x, y, c=color, s=500)
savefig('question3.jpg')
plt.show()