from .common import *

########################################################################################################
## Solver Methods

def Central_Difference_Step(P, dPdt, previous_step_info, i, Delta_t):
    P_before_update = P
    if i < 2:
        P = P + Delta_t * dPdt(P)
    else:
        P_2_steps_prior = previous_step_info["P_2_steps_prior"]
        P = P_2_steps_prior + 2 * Delta_t * dPdt(P)
    
    previous_step_info["P_2_steps_prior"] = P_before_update  # Corrected variable name
    
    return P, previous_step_info


# def Backwards_Difference_Step(P, dPdt, previous_step_info, i, Delta_t):
#     """
#     Implements the implicit Euler (backward difference) method:
    
#     (I - delta_t * dPdt) P_new = P_old
    
#     We solve for P_new by inverting (I - delta_t * dPdt).
#     """
#     I = np.eye(P.shape[0])  # Identity matrix
#     A = I - Delta_t * dPdt(P)  # Matrix for implicit method
    
#     # Solve the linear system A P_new = P
#     P_new = np.linalg.solve(A, P)
    
#     return P_new, previous_step_info

# def Backwards_Difference_Step(P, dPdt, previous_step_info, i, Delta_t):
#     I = np.eye(P.shape[0])
    
#     # Function to compute Ax for solver
#     def matvec(X):
#         return (I - Delta_t * dPdt(P)) @ X  # Use matrix-vector multiplication
    
#     # Use GMRES solver to solve (I - Δt dPdt) P_new = P
#     P_new, _ = gmres(matvec, P.flatten())  
#     P_new = P_new.reshape(P.shape)  # Reshape back to matrix form
    
#     return P_new, previous_step_info


# def Crank_Nicolson_Step(P, dPdt, previous_step_info, i, Delta_t):
#     I = np.eye(P.shape[0])
    
#     # Compute dPdt(P) only once to save computation
#     dPdt_P = dPdt(P)
    
#     # Matrix for implicit step
#     A = I - (Delta_t / 2) * dPdt_P
#     B = I + (Delta_t / 2) * dPdt_P
    
#     # Solve A P_new = B P
#     P_new = np.linalg.solve(A, B @ P)
    
#     return P_new, previous_step_info



# def Backwards_Difference_Step(P, dPdt, previous_step_info, i, Delta_t, tol=1e-6, max_iters_sub=100):
#     """
#     Implements the implicit (backward) Euler method using fixed-point iteration.
    
#     Solves: P_new = P_old + Delta_t * dPdt(P_new)
    
#     Parameters:
#     - P: Current P matrix
#     - dPdt: Function computing dP/dt given P
#     - previous_step_info: Dictionary storing iteration history
#     - i: Current time step index
#     - Delta_t: Time step size
#     - tol: Convergence tolerance
#     - max_iters: Maximum iterations for fixed-point method
    
#     Returns:
#     - P_new: Updated P matrix
#     - previous_step_info: Updated iteration history
#     """
#     P_new = P.copy()  # Initial guess (fixed-point iteration starts at P_old)
    
#     for _ in range(max_iters_sub):
#         P_next = P + Delta_t * dPdt(P_new)  # Fixed-point iteration update
        
#         # Check for convergence
#         if np.linalg.norm(P_next - P_new) < tol:
#             break
            
#         P_new = P_next  # Update for next iteration

#     return P_new, previous_step_info


# def Crank_Nicholson_Step(P, dPdt, previous_step_info, i, Delta_t, tol=1e-6, max_iters_sub=100):
#     """
#     Implements the Crank-Nicholson (trapezoidal) method using fixed-point iteration.
    
#     Solves: P_new = P_old + (Δt/2) * (dPdt(P_old) + dPdt(P_new))
    
#     Parameters:
#     - P: Current P matrix
#     - dPdt: Function computing dP/dt given P
#     - previous_step_info: Dictionary storing iteration history
#     - i: Current time step index
#     - Delta_t: Time step size
#     - tol: Convergence tolerance
#     - max_iters: Maximum iterations for fixed-point method
    
#     Returns:
#     - P_new: Updated P matrix
#     - previous_step_info: Updated iteration history
#     """
#     P_new = P.copy()  # Initial guess (fixed-point iteration starts at P_old)
#     dPdt_P = dPdt(P)  # Compute once for efficiency
    
#     for _ in range(max_iters_sub):
#         P_next = P + (Delta_t / 2) * (dPdt_P + dPdt(P_new))  # Fixed-point update
        
#         # Check for convergence
#         if np.linalg.norm(P_next - P_new) < tol:
#             break
            
#         P_new = P_next  # Update for next iteration

#     return P_new, previous_step_info

# def Backwards_Difference_Step(P, dPdt, previous_step_info, i, Delta_t):
#     """
#     Implements the stable backward Euler method via matrix inversion.

#     Solves: (I - Δt * D) P_new = P
    
#     Parameters:
#     - P: Current P matrix
#     - dPdt: Function computing dP/dt given P
#     - previous_step_info: Dictionary storing iteration history
#     - i: Current time step index
#     - Delta_t: Time step size
    
#     Returns:
#     - P_new: Updated P matrix
#     - previous_step_info: Updated iteration history
#     """
#     size = P.shape[0] * P.shape[1]  # Flattened size
#     I = np.eye(size)  # Identity matrix

#     # Compute Jacobian D numerically if needed
#     def D_operator(P_vec):
#         P_mat = P_vec.reshape(P.shape)  # Convert back to matrix
#         return dPdt(P_mat).flatten()  # Apply dPdt and flatten result

#     D = np.array([D_operator(np.eye(size)[j].reshape(P.shape)) for j in range(size)]).T

#     # Solve the linear system (I - Δt * D) P_new = P
#     P_new_vec = np.linalg.solve(I - Delta_t * D, P.flatten())
#     P_new = P_new_vec.reshape(P.shape)

#     return P_new, previous_step_info


# def Crank_Nicholson_Step(P, dPdt, previous_step_info, i, Delta_t):
#     """
#     Implements the stable Crank-Nicholson method via matrix inversion.

#     Solves: (I - Δt/2 * D) P_new = (I + Δt/2 * D) P
    
#     Parameters:
#     - P: Current P matrix
#     - dPdt: Function computing dP/dt given P
#     - previous_step_info: Dictionary storing iteration history
#     - i: Current time step index
#     - Delta_t: Time step size
    
#     Returns:
#     - P_new: Updated P matrix
#     - previous_step_info: Updated iteration history
#     """
#     size = P.shape[0] * P.shape[1]  # Flattened size
#     I = np.eye(size)  # Identity matrix

#     # Compute Jacobian D numerically
#     def D_operator(P_vec):
#         P_mat = P_vec.reshape(P.shape)
#         return dPdt(P_mat).flatten()

#     D = np.array([D_operator(np.eye(size)[j].reshape(P.shape)) for j in range(size)]).T

#     # Solve the linear system (I - Δt/2 * D) P_new = (I + Δt/2 * D) P
#     A = I - (Delta_t / 2) * D
#     B = I + (Delta_t / 2) * D
#     P_new_vec = np.linalg.solve(A, B @ P.flatten())
#     P_new = P_new_vec.reshape(P.shape)

#     return P_new, previous_step_info




# import numpy as np
# from scipy.sparse.linalg import gmres, LinearOperator

# def Crank_Nicholson_Step(P, dPdt, previous_step_info, i, Delta_t):
#     """
#     Implements the Crank-Nicholson method:
#     (I - Δt/2 * D) P_new = (I + Δt/2 * D) P_old
#     """
#     P_shape = P.shape
#     P_size = P.size  # Total number of elements in P

#     # Convert P to 1D array (flatten to 1D) if it is a np.matrix
#     if isinstance(P, np.matrix):
#         P_vec = P.A1  # .A1 gives the flattened 1D array from np.matrix
#     else:
#         P_vec = P.flatten()  # Use .flatten() for np.ndarray

#     # Compute RHS: (I + Δt/2 * D) P_old
#     RHS = P + (Delta_t / 2) * dPdt(P)
#     if isinstance(RHS, np.matrix):
#         RHS = RHS.A1  # Flatten if it's a np.matrix
#     else:
#         RHS = RHS.flatten()  # Flatten if it's a np.ndarray

#     # Define the matrix-vector product for A_linop
#     def matvec(P_vec_flat):
#         if P_vec_flat.size != P_size:
#             raise ValueError(f"Expected size {P_size}, got {P_vec_flat.size}")
#         P_mat = P_vec_flat.reshape(P_shape)
#         return (P_mat - (Delta_t / 2) * dPdt(P_mat)).flatten()  # Ensure output is 1D array

#     A_linop = LinearOperator((P_size, P_size), matvec=matvec)

#     # Solve for P_new using GMRES with adjusted tolerances and maxiter
#     P_new_vec, exit_code = gmres(A_linop, RHS, atol=1e-5, rtol=1e-5, maxiter=1000)

#     if exit_code != 0:
#         print(f"Warning: GMRES did not converge at step {i}, exit code {exit_code}")

#     return P_new_vec.reshape(P_shape), None


# def Backwards_Difference_Step(P, dPdt, previous_step_info, i, Delta_t):
#     """
#     Implements the Backward Difference method:
#     P_new = P_old + Δt * D P_new
#     Rearranged as: (I - Δt * D) P_new = P_old
#     """
#     P_shape = P.shape
#     P_size = P.size  # Total number of elements in P

#     # Convert P to 1D array (flatten to 1D) if it is a np.matrix
#     if isinstance(P, np.matrix):
#         P_vec = P.A1  # .A1 gives the flattened 1D array from np.matrix
#     else:
#         P_vec = P.flatten()  # Use .flatten() for np.ndarray

#     # Compute RHS: P_old
#     RHS = P_vec

#     # Define the matrix-vector product for A_linop (I - Δt * D)
#     def matvec(P_vec_flat):
#         if P_vec_flat.size != P_size:
#             raise ValueError(f"Expected size {P_size}, got {P_vec_flat.size}")
#         P_mat = P_vec_flat.reshape(P_shape)
#         return (P_mat - Delta_t * dPdt(P_mat)).flatten()  # Ensure output is 1D array

#     A_linop = LinearOperator((P_size, P_size), matvec=matvec)

#     # Solve for P_new using GMRES with adjusted tolerances and maxiter
#     P_new_vec, exit_code = gmres(A_linop, RHS, atol=1e-5, rtol=1e-5, maxiter=1000)

#     if exit_code != 0:
#         print(f"Warning: GMRES did not converge at step {i}, exit code {exit_code}")

#     return P_new_vec.reshape(P_shape), None


# import numpy as np
# from scipy.sparse.linalg import gmres, LinearOperator

# def Crank_Nicholson_Step(P, dPdt, previous_step_info, i, Delta_t):
#     """
#     Implements the Crank-Nicholson method:
#     (I - Δt/2 * D) P_new = (I + Δt/2 * D) P_old
#     """
#     P_shape = P.shape
#     P_size = P.size  # Total number of elements in P

#     # Convert P to 1D array (flatten to 1D) if it is a np.matrix
#     if isinstance(P, np.matrix):
#         P_vec = P.A1  # .A1 gives the flattened 1D array from np.matrix
#     else:
#         P_vec = P.flatten()  # Use .flatten() for np.ndarray

#     # Compute RHS: (I + Δt/2 * D) P_old
#     RHS = P + (Delta_t / 2) * dPdt(P)
    
#     if isinstance(RHS, np.matrix):
#         RHS = RHS.A1  # Flatten if it's a np.matrix
#     else:
#         RHS = RHS.flatten()  # Flatten if it's a np.ndarray

#     # Define the matrix-vector product for A_linop
#     def matvec(P_vec_flat):
#         if P_vec_flat.size != P_size:
#             raise ValueError(f"Expected size {P_size}, got {P_vec_flat.size}")
#         P_mat = P_vec_flat.reshape(P_shape)
#         return (P_mat - (Delta_t / 2) * dPdt(P_mat)).flatten()  # Ensure output is 1D array

#     A_linop = LinearOperator((P_size, P_size), matvec=matvec)

#     # Solve for P_new using GMRES with adjusted tolerances and maxiter
#     P_new_vec, exit_code = gmres(A_linop, RHS, atol=1e-5, rtol=1e-5, maxiter=1000)

#     if exit_code != 0:
#         print(f"Warning: GMRES did not converge at step {i}, exit code {exit_code}")

#     # Return P_new as np.matrix to maintain consistency with the input type
#     return np.matrix(P_new_vec.reshape(P_shape)), None


# def Backward_Difference_Step(P, dPdt, previous_step_info, i, Delta_t):
#     """
#     Implements the Backward Difference method:
#     P_new = P_old + Δt * D P_new
#     Rearranged as: (I - Δt * D) P_new = P_old
#     """
#     P_shape = P.shape
#     P_size = P.size  # Total number of elements in P

#     # Convert P to 1D array (flatten to 1D) if it is a np.matrix
#     if isinstance(P, np.matrix):
#         P_vec = P.A1  # .A1 gives the flattened 1D array from np.matrix
#     else:
#         P_vec = P.flatten()  # Use .flatten() for np.ndarray

#     # Compute RHS: P_old
#     RHS = P_vec

#     # Define the matrix-vector product for A_linop (I - Δt * D)
#     def matvec(P_vec_flat):
#         if P_vec_flat.size != P_size:
#             raise ValueError(f"Expected size {P_size}, got {P_vec_flat.size}")
#         P_mat = P_vec_flat.reshape(P_shape)
#         return (P_mat - Delta_t * dPdt(P_mat)).flatten()  # Ensure output is 1D array

#     A_linop = LinearOperator((P_size, P_size), matvec=matvec)

#     # Solve for P_new using GMRES with adjusted tolerances and maxiter
#     P_new_vec, exit_code = gmres(A_linop, RHS, atol


import numpy as np
from scipy.sparse.linalg import gmres, LinearOperator

def Crank_Nicholson_Step(P, dPdt, previous_step_info, i, Delta_t):
    """
    Implements the Crank-Nicholson method:
    (I - Δt/2 * D) P_new = (I + Δt/2 * D) P_old
    """
    P_shape = P.shape
    P_size = P.size  # Total number of elements in P

    # Convert P to 1D array (flatten to 1D) if it is a np.matrix
    if isinstance(P, np.matrix):
        P_vec = P.A1  # .A1 gives the flattened 1D array from np.matrix
    else:
        P_vec = P.flatten()  # Use .flatten() for np.ndarray

    # Compute RHS: (I + Δt/2 * D) P_old
    RHS = P + (Delta_t / 2) * dPdt(P)
    
    if isinstance(RHS, np.matrix):
        RHS = RHS.A1  # Flatten if it's a np.matrix
    else:
        RHS = RHS.flatten()  # Flatten if it's a np.ndarray

    # Define the matrix-vector product for A_linop
    def matvec(P_vec_flat):
        if P_vec_flat.size != P_size:
            raise ValueError(f"Expected size {P_size}, got {P_vec_flat.size}")
        P_mat = P_vec_flat.reshape(P_shape)
        return (P_mat - (Delta_t / 2) * dPdt(P_mat)).flatten()  # Ensure output is 1D array

    A_linop = LinearOperator((P_size, P_size), matvec=matvec)

    # Solve for P_new using GMRES with adjusted tolerances and maxiter
    P_new_vec, exit_code = gmres(A_linop, RHS, atol=1e-5, rtol=1e-5, maxiter=1000)

    if exit_code != 0:
        print(f"Warning: GMRES did not converge at step {i}, exit code {exit_code}")

    # Return P_new as np.matrix to maintain consistency with the input type
    return np.matrix(P_new_vec.reshape(P_shape)), None


def Backwards_Difference_Step(P, dPdt, previous_step_info, i, Delta_t):
    """
    Implements the Backward Difference method:
    P_new = P_old + Δt * D P_new
    Rearranged as: (I - Δt * D) P_new = P_old
    """
    P_shape = P.shape
    P_size = P.size  # Total number of elements in P

    # Convert P to 1D array (flatten to 1D) if it is a np.matrix
    if isinstance(P, np.matrix):
        P_vec = P.A1  # .A1 gives the flattened 1D array from np.matrix
    else:
        P_vec = P.flatten()  # Use .flatten() for np.ndarray

    # Compute RHS: P_old
    RHS = P_vec

    # Define the matrix-vector product for A_linop (I - Δt * D)
    def matvec(P_vec_flat):
        if P_vec_flat.size != P_size:
            raise ValueError(f"Expected size {P_size}, got {P_vec_flat.size}")
        P_mat = P_vec_flat.reshape(P_shape)
        return (P_mat - Delta_t * dPdt(P_mat)).flatten()  # Ensure output is 1D array

    A_linop = LinearOperator((P_size, P_size), matvec=matvec)

    # Solve for P_new using GMRES with adjusted tolerances and maxiter
    P_new_vec, exit_code = gmres(A_linop, RHS, atol=1e-5, rtol=1e-5, maxiter=1000)

    if exit_code != 0:
        print(f"Warning: GMRES did not converge at step {i}, exit code {exit_code}")

    # Return P_new as np.matrix to maintain consistency with the input type
    return np.matrix(P_new_vec.reshape(P_shape)), None








Solver_Func_Dict = {
    "Forwards Difference":   lambda P, dPdt, previous_step_info, i, Delta_t: (P + Delta_t * dPdt(P), previous_step_info),
    "Central Difference":    Central_Difference_Step,
    "Backwards Difference":  Backwards_Difference_Step,
    "Crank-Nicholson":       Crank_Nicholson_Step,
}