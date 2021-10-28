import matplotlib.pyplot as plt
import numpy as np

xy = [list(range(10)), np.sin(range(10))]
#
# plt.plot(*xy)
# plt.show()
def dayin(x,y):
    print(x)
    print(y)

print(*xy)
dayin(*xy)
print(abs(12-123))