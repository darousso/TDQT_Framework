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
                print(f"Convergence achieved in first iteration at step {step}, doubling Delta_t.")
            break
        else:
            P_guess = P_new  # Update P_guess for the next iteration

    # If converged, update P, otherwise restore P to the original input
    if converged:
        previous_step_info["converged"] = True
        return P_new, previous_step_info, Delta_t
    else:
        # Restore the original P (no update to time step or t_vec_scaled)
        print(f"Convergence failed at step {step}, restoring P, halving Delta_t.")
        # If convergence fails after max_iterations, halve Delta_t and roll back to the previous time step
        Delta_t /= 2
        previous_step_info["converged"] = False
        return P_fixed, previous_step_info, Delta_t  # Return the original P to keep time unchanged


def RK4_Step(P, dPdt, previous_step_info, i, Delta_t,N_bas):
    """
    Implements the 4th-order Runge-Kutta method for evolving P.

    Parameters:
    - P: Current P matrix
    - dPdt: Function computing dP/dt given P
    - previous_step_info: Dictionary storing iteration history
    - i: Current time step index
    - Delta_t: Time step size
    
    Returns:
    - P_new: Updated P matrix
    - previous_step_info: Updated iteration history
    """
    k1 = dPdt(P)
    k2 = dPdt(P + (Delta_t / 2) * k1)
    k3 = dPdt(P + (Delta_t / 2) * k2)
    k4 = dPdt(P + Delta_t * k3)

    P_new = P + (Delta_t / 6) * (k1 + 2*k2 + 2*k3 + k4)

    # Store intermediate steps if needed (e.g., for adaptive stepping later)
    previous_step_info[i] = {
        "k1": k1, "k2": k2, "k3": k3, "k4": k4, "P_new": P_new
    }

    return P_new, previous_step_info,Delta_t



from scipy.integrate import solve_ivp

def SciPy_Adaptive_Step(P, dPdt, previous_step_info, i, Delta_t,N_bas):
    """
    Uses SciPy's adaptive ODE solver (solve_ivp) with stiff-compatible solvers.
    
    Parameters:
    - P: Current matrix
    - dPdt: Function computing dP/dt
    - previous_step_info: Dictionary storing iteration history
    - i: Current time step index
    - Delta_t: Time step size (used as initial step)
    
    Returns:
    - P_new: Updated matrix
    - previous_step_info: Updated history
    """
    
    # Flatten P into a vector because solve_ivp expects 1D input
    P_shape = P.shape
    P_vec = P.flatten()

    # Define a wrapper for dPdt that operates on a vector
    def dPdt_vec(t, P_flat):
        P_mat = P_flat.reshape(P_shape)  # Reshape back into matrix form
        dPdt_mat = dPdt(P_mat)  # Compute dP/dt in matrix form
        return dPdt_mat.flatten()  # Return as a flat vector

    # Solve using an implicit method (Radau is stable for stiff problems)
    sol = solve_ivp(dPdt_vec, (0, Delta_t), P_vec, method="Radau", atol=1e-8, rtol=1e-6)

    # Extract the final solution and reshape it back into a matrix
    P_new = sol.y[:, -1].reshape(P_shape)

    # Store solution history if needed
    previous_step_info[i] = {"P_new": P_new, "solver_status": sol.message}

    return P_new, previous_step_info, Delta_t



def RK4_Adaptive_Step(P, dPdt, previous_step_info, i, Delta_t, N_bas, tolerance=1e-4, max_factor=2, min_factor=0.5):
    """
    Implements the 4th-order Runge-Kutta method for evolving P with adaptive time-stepping.

    Parameters:
    - P: Current P matrix
    - dPdt: Function computing dP/dt given P
    - previous_step_info: Dictionary storing iteration history
    - i: Current time step index
    - Delta_t: Current time step size
    - tolerance: Desired tolerance for the error estimate (default is 1e-6)
    - max_factor: Factor by which the step size can be increased (default is 2)
    - min_factor: Factor by which the step size can be decreased (default is 0.5)
    
    Returns:
    - P_new: Updated P matrix
    - previous_step_info: Updated iteration history
    - Delta_t: New time step size (adapted)
    """
    
    P_orig = P.copy()
    
    # First Runge-Kutta step (with current Delta_t)
    k1 = dPdt(P)
    k2 = dPdt(P + (Delta_t / 2) * k1)
    k3 = dPdt(P + (Delta_t / 2) * k2)
    k4 = dPdt(P + Delta_t * k3)

    P_new = P + (Delta_t / 6) * (k1 + 2*k2 + 2*k3 + k4)

    # Second Runge-Kutta step (with half of Delta_t, for error estimation)
    Delta_t_half = Delta_t / 2
    k1_half = dPdt(P)
    k2_half = dPdt(P + (Delta_t_half / 2) * k1_half)
    k3_half = dPdt(P + (Delta_t_half / 2) * k2_half)
    k4_half = dPdt(P + Delta_t_half * k3_half)

    P_half_step = P + (Delta_t_half / 6) * (k1_half + 2*k2_half + 2*k3_half + k4_half)

    # Estimate the error between the full step and the half-step result
    error_estimate = np.linalg.norm(P_new - P_half_step)/N_bas

    # Check if the error estimate is within tolerance
    if error_estimate < tolerance:
        # Accept the full step
        previous_step_info = {
#             "P_new": P_new,
#             "error_estimate": error_estimate,
            "converged":True,
        }

        # Increase the step size if error is sufficiently small
        Delta_t_new = Delta_t * max_factor
    else:
        # Reject the step and try again with a smaller time step
        Delta_t_new = Delta_t * min_factor
        P_new = P_orig  # Restore previous step if rejected

        previous_step_info = {
#             "P_new": P_new,
#             "error_estimate": error_estimate,
            "converged": False
        }

    return P_new, previous_step_info, Delta_t_new


# def RK45_Step(P, dPdt, previous_step_info, i, Delta_t, N_bas,tol=1e-4):
#     """
#     Implements the 4(5)th-order Runge-Kutta-Fehlberg method (RK45) for evolving P.

#     Parameters:
#     - P: Current P matrix
#     - dPdt: Function computing dP/dt given P
#     - previous_step_info: Dictionary storing iteration history
#     - i: Current time step index
#     - Delta_t: Time step size
#     - tol: Desired accuracy (tolerance)
    
#     Returns:
#     - P_new: Updated P matrix
#     - previous_step_info: Updated iteration history
#     - Delta_t_new: Adjusted time step size
#     """
#     # Perform RK45 method using two orders (4th and 5th)
#     # Compute k1, k2, k3, k4, k5, k6 using the RK45 coefficients

#     k1 = dPdt(P)
#     k2 = dPdt(P + (Delta_t / 4) * k1)
#     k3 = dPdt(P + (Delta_t / 8) * k1 + (Delta_t / 8) * k2)
#     k4 = dPdt(P + (Delta_t / 2) * k2)
#     k5 = dPdt(P - (Delta_t / 2) * k2 + Delta_t * k3)
#     k6 = dPdt(P + Delta_t * k3)

#     # Compute 4th and 5th order estimates for P
#     P4 = P + (Delta_t / 6) * (k1 + 2*k2 + 2*k3 + k4)
#     P5 = P + (Delta_t / 6) * (k1 + 4*k2 + 2*k3 + k4)

#     # Estimate the error between the 4th and 5th order solutions
#     error = np.linalg.norm(P5 - P4)  / N_bas

#     # Adjust time step size based on error estimate
#     if error < tol:
#         Delta_t_new = Delta_t * min(5, (tol / error) ** 0.2)
#         previous_step_info["converged"] = True
#         P_new = P5  # Accept the 5th order solution
#     else:
#         Delta_t_new = Delta_t * max(0.2, (tol / error) ** 0.25)
#         previous_step_info["converged"] = False
#         P_new = P  # Reject this step, P remains the same

#     # Store the immediate previous step information in the dictionary
#     previous_step_info["P_prev"] = P  # Store the previous step

#     return P_new, previous_step_info, Delta_t_new

def RK45_Step(P, dPdt, previous_step_info, i, Delta_t, N_bas, tol=1e-4):
    """
    Implements the 4(5)th-order Runge-Kutta-Fehlberg method (RK45) for evolving P.

    Parameters:
    - P: Current P matrix
    - dPdt: Function computing dP/dt given P
    - previous_step_info: Dictionary storing iteration history
    - i: Current time step index
    - Delta_t: Time step size
    - tol: Desired accuracy (tolerance)
    
    Returns:
    - P_new: Updated P matrix
    - previous_step_info: Updated iteration history
    - Delta_t_new: Adjusted time step size
    """
    # Standard RK45 coefficients
    a2, a3, a4, a5, a6 = 1/4, 3/8, 12/13, 1, 1/2
    b21 = 1/4
    b31, b32 = 3/32, 9/32
    b41, b42, b43 = 1932/2197, -7200/2197, 7296/2197
    b51, b52, b53, b54 = 439/216, -8, 3680/513, -845/4104
    b61, b62, b63, b64, b65 = -8/27, 2, -3544/2565, 1859/4104, -11/40

    c1, c3, c4, c5, c6 = 16/135, 6656/12825, 28561/56430, -9/50, 2/55  # 5th-order
    d1, d3, d4, d5, d6 = 25/216, 1408/2565, 2197/4104, -1/5, 0  # 4th-order

    # Compute Runge-Kutta stages
    k1 = dPdt(P)
    k2 = dPdt(P + Delta_t * (b21 * k1))
    k3 = dPdt(P + Delta_t * (b31 * k1 + b32 * k2))
    k4 = dPdt(P + Delta_t * (b41 * k1 + b42 * k2 + b43 * k3))
    k5 = dPdt(P + Delta_t * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4))
    k6 = dPdt(P + Delta_t * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5))

    # Compute 5th and 4th order estimates
    P5 = P + Delta_t * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 + c6 * k6)
    P4 = P + Delta_t * (d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6)

    # Compute relative error
    norm_P = max(np.linalg.norm(P), np.linalg.norm(P5), 1e-10)  # Avoid division by zero
    error = np.linalg.norm(P5 - P4) / norm_P

    # Adjust time step size based on error estimate
    safety_factor = 0.9  # Standard safety factor
    if error < tol:
        Delta_t_new = Delta_t * min(5, safety_factor * (tol / error) ** 0.2)
        P_new = P5  # Accept step
        previous_step_info["converged"] = True
    else:
        Delta_t_new = Delta_t * max(0.2, safety_factor * (tol / error) ** 0.25)
        P_new = P  # Reject step

    # Store the immediate previous step information in the dictionary
    previous_step_info["P_prev"] = P
    previous_step_info["error"] = error
    previous_step_info["Delta_t_new"] = Delta_t_new

    return P_new, previous_step_info, Delta_t_new


#########################################################################
Solver_Func_Dict = {
    "Forwards Difference":   lambda P, dPdt, previous_step_info, i, Delta_t, NBas: (P + Delta_t * dPdt(P), previous_step_info,Delta_t),
    "Central Difference":    Central_Difference_Step,
    "Backwards Difference":  Backwards_Difference_Step,
    "Crank-Nicholson":       Crank_Nicholson_Step,
    "Implicit Euler":        Implicit_Euler_Variable_Timesteps,
    "RK4":                   RK4_Step,
    "SciPy_Adaptive":        SciPy_Adaptive_Step,
    "RK4_Adaptive_Step":     RK4_Adaptive_Step,
    "RK45":                  RK45_Step,
}