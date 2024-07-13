import numpy as np

# Problem 1
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])
print("Problem 1")
print(a + b)
print(a - b)


# Problem 2
print("Problem 2")
c = np.array([[1, 2], [3, 4]])
d = np.array([[5, 6], [7, 8]])
print(c + d)
print(c - d)

# Problem 3
print("Problem 3")
print(np.dot(a, b))

# Problem 4
print("Problem 4")
e = np.array([[1, 2, 3], [4, 5, 6]])
f = np.array([[7, 8, 9, 10], [11, 12, 13, 14], [15, 16, 17, 18]])
print(np.matmul(e, f))

# Problem 5
print("Problem 5")
g = np.array([1, 1, 2])
print(np.linalg.norm(g))

# Problem 6
print("Problem 6")
print(np.transpose(c))
