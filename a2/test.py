import numpy as np

def calculate_s_matrix():
    M = np.array([[0.16, 0, 0], [0, 0.16, 0], [0, 0, 0.16]])

    G = np.array([[-1, 1, 0], [-1, 0, 1], [0, 0, 0]])
    
    # Compute G^T * M
    intermediate_result = np.dot(G.T, M)
    
    # Compute S = (G^T * M) * G
    S = np.dot(intermediate_result, G)
    
    # The result matrix S now stores the cotangent matrix S
    return S

print(calculate_s_matrix())