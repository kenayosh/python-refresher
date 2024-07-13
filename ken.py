import numpy as np

A = [1,2,3,4]
B = [5,6,7,8]
C = np.zeros((4,2))
print(C)

for i in range (4):
    for j in range (2):
        C[i,j] = A[i] + B[i]
        C[i,j+1] = 6
        
        
print(C)
    
    
