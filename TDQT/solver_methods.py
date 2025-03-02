from .common import *

########################################################################################################
## Solver Methods

def Central_Difference_Step(P, dPdt, previous_step_info, i, Delta_t,Nbas):
    P_before_update = P
    if i < 2:
        P = P + Delta_t * dPdt(P)
    else:
        P_2_steps_prior = previous_step_info["P_2_steps_prior"]
        P = P_2_steps_prior + 2 * Delta_t * dPdt(P)
    
    previous_step_info["P_2_steps_prior"] = P_before_update  # Corrected variable name
    
    return P, previous_step_info,Delta_t


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


# import numpy as np
# from scipy.sparse.linalg import gmres, LinearOperator

def Crank_Nicholson_Step(P, dPdt, previous_step_info, i, Delta_t,Nbas):
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
    return np.matrix(P_new_vec.reshape(P_shape)), None,Delta_t


def Backwards_Difference_Step(P, dPdt, previous_step_info, i, Delta_t,N_bas):
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
    return np.matrix(P_new_vec.reshape(P_shape)), None,Delta_t




# def Implicit_Euler_Variable_Timesteps(P, dPdt, previous_step_info, i, Delta_t):
#     """
#     Implements the Implicit Euler method:
#     P_new = P_old + Δt * dPdt(P_new)
#     Rearranged as: P_new - Δt * dPdt(P_new) = P_old
#     Uses a fixed-point iteration (Picard iteration) to solve the implicit equation.
    
#     Adaptive time-stepping is based on the convergence of the iteration.
#     """
#     P_shape = P.shape
#     P_size = P.size  # Total number of elements in P

#     # Convert P to 1D array (flatten to 1D) if it is a np.matrix
#     if isinstance(P, np.matrix):
#         P_vec = P.A1  # .A1 gives the flattened 1D array from np.matrix
#     else:
#         P_vec = P.flatten()  # Use .flatten() for np.ndarray

#     # Initial guess for P_new
#     P_new_vec = P_vec.copy()
    
#     # Fixed-point iteration (Picard iteration) for solving the implicit equation
#     max_iter = 100
#     tol = 1e-8
#     iteration = 0
#     converged = False

#     while iteration < max_iter and not converged:
#         # Compute dP/dt at the current guess
#         dPdt_new = dPdt(P_new_vec.reshape(P_shape))
        
#         # Update P_new_vec using the implicit Euler equation
#         P_new_vec_next = P_vec + Delta_t * dPdt_new.flatten()
        
#         # Check convergence (we stop if the change is smaller than the tolerance)
#         if np.linalg.norm(P_new_vec_next - P_new_vec) < tol:
#             converged = True
#         else:
#             P_new_vec = P_new_vec_next
        
#         iteration += 1

#     if not converged:
#         print(f"Warning: Implicit Euler did not converge in {max_iter} iterations at step {i}")
    
#     # Adaptive time-stepping: If convergence was slow, decrease Delta_t
#     if iteration > 10:  # If many iterations were needed, reduce time step
#         Delta_t *= 0.5
#     else:  # If convergence was quick, slightly increase time step
#         Delta_t *= 1.1
    
#     # Return the new P as a np.matrix to maintain consistency with input type
#     return np.matrix(P_new_vec_next.reshape(P_shape)), Delta_t

# def Implicit_Euler_Variable_Timesteps(P, dPdt, previous_step_info, i, Delta_t, NBas):
#     """
#     Implements the Implicit Euler method with the fixed-point iteration and adaptive time steps.
    
#     The method follows the equation:
#         P(t + Δt) = P(t) + Δt * f(t + Δt, P(t + Δt))
    
#     Convergence is determined by:
#         ||P_{i+1}(t+Δt) - P_i(t+Δt)|| / NBas < 10^-4
    
#     Time step is adjusted according to the convergence behavior:
#         - If converges in the first iteration, double the time step.
#         - If fails to converge after 5 iterations, halve the time step.
#     """
#     P_shape = P.shape
#     P_size = P.size  # Total number of elements in P

#     # Convert P to 1D array (flatten to 1D) if it is a np.matrix
#     if isinstance(P, np.matrix):
#         P_vec = P.A1  # .A1 gives the flattened 1D array from np.matrix
#     else:
#         P_vec = P.flatten()  # Use .flatten() for np.ndarray

#     # Initial guess for P_new
#     P_new_vec = P_vec.copy()
    
#     # Fixed-point iteration (Picard iteration) for solving the implicit equation
#     max_iter = 5  # Max 5 iterations for convergence
#     tol = 1e-4  # Convergence criterion
#     iteration = 0
#     converged = False

#     previous_P_vec = P_vec.copy()  # Store the initial guess for convergence check

#     while iteration < max_iter and not converged:
#         # Compute f(t+Δt, P(t+Δt)) using the previous guess for P_new
#         dPdt_new = dPdt(P_new_vec.reshape(P_shape))
        
#         # Update P_new_vec using the implicit Euler equation
#         P_new_vec_next = P_vec + Delta_t * dPdt_new.flatten()
        
#         # Check convergence (we stop if the change is smaller than the tolerance)
#         change = np.linalg.norm(P_new_vec_next - P_new_vec) / NBas
#         if change < tol:
#             converged = True
#         else:
#             P_new_vec = P_new_vec_next
        
#         iteration += 1

#     # Handle adaptive time-stepping based on convergence
#     if not converged:
#         print(f"Warning: Implicit Euler did not converge in {max_iter} iterations at step {i}")
#         # Halve time step and roll back the density matrix
#         Delta_t *= 0.5
#         P_new_vec = previous_P_vec  # Roll back to the previous state
#     elif iteration == 1:
#         # If converged in the first iteration, double the time step
#         Delta_t *= 2

#     # Return the new P as a np.matrix to maintain consistency with input type
#     return np.matrix(P_new_vec.reshape(P_shape)), previous_step_info,Delta_t


def Implicit_Euler_Variable_Timesteps(P, dPdt, previous_step_info, step, Delta_t, NBas):
    """
    Perform one step of the Implicit Euler method with a convergence check, keeping P(t) fixed.
    """
    P_fixed = P.copy()  # This is the fixed P(t) from the previous time step
    P_guess = P.copy()  # Start with an initial guess for P(t + Delta_t)
    converged = False
    max_iterations = 5  # Max iterations for convergence check (as per the provided method)
    tolerance = 1e-4  # Convergence tolerance (as per the provided method)

    # Start iterative process for implicit update
    for i in range(max_iterations):
        # Update P_guess using implicit Euler: P_guess = P_fixed + Delta_t * dPdt(P_guess)
        P_new = P_fixed + Delta_t * dPdt(P_guess)

        # Check convergence using Euclidean norm
        if np.linalg.norm(P_new - P_guess) / NBas < tolerance:
            converged = True
            # If convergence is achieved in the first iteration, double Delta_t
            if i == 0:
                Delta_t *= 2
            break
        else:
            P_guess = P_new  # Update P_guess for the next iteration

    # If converged, update P, otherwise restore P to the original input
    if converged:
        previous_step_info["converged"] = True
        return P_new, previous_step_info, Delta_t
    else:
        # Restore the original P (no update to time step or t_vec_scaled)
        print(f"Convergence failed at step {step}, restoring P.")
        # If convergence fails after max_iterations, halve Delta_t and roll back to the previous time step
        Delta_t /= 2
        previous_step_info["converged"] = False
        return P_fixed, previous_step_info, Delta_t  # Return the original P to keep time unchanged







#########################################################################
Solver_Func_Dict = {
    "Forwards Difference":   lambda P, dPdt, previous_step_info, i, Delta_t, NBas: (P + Delta_t * dPdt(P), previous_step_info,Delta_t),
    "Central Difference":    Central_Difference_Step,
    "Backwards Difference":  Backwards_Difference_Step,
    "Crank-Nicholson":       Crank_Nicholson_Step,
    "Implicit Euler":        Implicit_Euler_Variable_Timesteps,
}