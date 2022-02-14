#!/bin/env python3


from create_random_disturb import gen_random_emat

import numpy as np

cell0 = np.array([[10, 0, 0], [0, 10, 0], [0, 0, 10]])

print(np.dot(cell0, gen_random_emat(1)))
