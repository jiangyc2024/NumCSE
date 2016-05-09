import numpy as np
# array creation
A = np.array([[1, 2], [3, 4]])  # default (row major) storage
B = np.array([[1, 2], [3, 4]], order='F')  # column major storage

# show internal storage
np.ravel(A, 'K')  # array elements as stored in memory: [1, 2, 3, 4]
np.ravel(B, 'K')  # array elements as stored in memory: [1, 3, 2, 4]

# nothing happens to the data on transpose, just the storage order changes
np.ravel(A.T, 'K')  # array elements as stored in memory: [1, 2, 3, 4]
np.ravel(B.T, 'K')  # array elements as stored in memory: [1, 3, 2, 4]

# storage order can be accessed by checking the array's flags
A.flags['C_CONTIGUOUS']  # True
B.flags['F_CONTIGUOUS']  # True
A.T.flags['F_CONTIGUOUS']  # True
B.T.flags['C_CONTIGUOUS']  # True
