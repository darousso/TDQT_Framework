from .common import *
from .utils import *
from .solver_methods import *

class Solver:
    
    def __init__(self, device, Gamma, Delta_t, max_iters, list_of_outputs,max_t=None):
        """Initializes functions needed for solver"""
        
        ####################################################
        ## importing stuff from device
        
        self.device = device  
        self.Gamma  = Gamma
        self.Delta_t = Delta_t
        self.max_iters = max_iters
        self.list_of_outputs= list_of_outputs
        self.max_t=max_t
        
        H=device.H
        S=device.S
        
        Ub=device.Ub
        Hs=device.Hs
        Ss=device.Ss
        
        U=device.U
        Hss=device.Hss
        
        Target_P0ss=device.Target_P0ss
        Target_P0s=device.Target_P0s
        Target_P0=device.Target_P0
        self.Initial_P0ss=device.Initial_P0ss
        self.Initial_P0s=device.Initial_P0s
        self.Initial_P0=device.Initial_P0
        
        Hs_REM=device.Hs_REM
        Hs_LEM=device.Hs_LEM
        
        size_full=device.size_full
        size_EM=device.size_EM
        size_lead_L=device.size_lead_L
        size_lead_R=device.size_lead_R
        
        self.size_EM=device.size_EM
        self.size_lead_L=device.size_lead_L
        self.size_lead_R=device.size_lead_R
        
        clean_matrices=device.clean_matrices
        
        ####################################################
        ## Other Relevant Matrices

        UbU=Ub@U
        UdUbd=U.H@Ub.H
        UbdUd=Ub.H@U.H
        Ud=U.H
        Ubd=Ub.H
            
        ####################################################
        ## Hadamard Masks for Density Matrix Target
        
        CONTACTS=0*Target_P0
        CONTACTS[:size_lead_L,:size_lead_L]=1
        CONTACTS[-size_lead_R:,-size_lead_R:]=1

        CORNERS=0*Target_P0
        CORNERS[:size_lead_L,:size_lead_L]=1
        CORNERS[-size_lead_R:,-size_lead_R:]=1
        CORNERS[:size_lead_L,-size_lead_R:]=1
        CORNERS[-size_lead_R:,:size_lead_L]=1

        CROSS=0*Target_P0
        CROSS[:size_lead_L,:]=1
        CROSS[:,:size_lead_L]=1
        CROSS[-size_lead_R:,:]=1
        CROSS[:,-size_lead_R:]=1
        CROSS=CROSS-CORNERS
        
        ####################################################
        ## Deal with scaling
        
        t_scale = hbar / abs(H).max()
        H_scaled = H / abs(H).max()
        Hss_scaled = Hss / abs(H).max()
        Gamma_scaled = Gamma * t_scale
        Delta_t_scaled = Delta_t / t_scale
        
        self.Gamma_scaled=Gamma_scaled
        self.Delta_t_scaled=Delta_t_scaled 
        self.t_scale=t_scale
        
        
        
        

        ####################################################
        ## Define derivative site
        
        # Precompute S^{-1} H_scaled and S^{-1} H_scaled.T using solve
        invS_H_scaled = np.linalg.solve(S, H_scaled)  # Equivalent to self.invS @ H_scaled
        H_scaled_invS = np.linalg.solve(S.T, H_scaled.T).T  # Equivalent to H_scaled @ self.invS

        def dPdt_Site_Scaled(P):
            # Compute commutator using precomputed S^{-1} * H
            COMMUTE = -1j * (invS_H_scaled @ P - P @ H_scaled_invS)

            if Gamma == 0:
                DAMP = 0
            else:
                # Solve linear system instead of computing inverses explicitly
                Pss = np.linalg.solve(UbU, P @ np.linalg.solve(UbdUd, np.eye(U.shape[0])))
                DAMP = -Gamma_scaled * UbU @ (Pss * (CORNERS + CROSS / 2) - Target_P0ss * CONTACTS) @ UdUbd

            return COMMUTE + DAMP
        
        self.dPdt_Site_Scaled=dPdt_Site_Scaled
        
        ####################################################
        ## Define outputs site
        
        def Get_J_All_Elements_Site(P):
            return np.array([np.imag(P[n, n+1]) * 2 * e / hbar * (H[n, n+1]) for n in range(size_full - 1)])
        
        def Get_J_Central_Site(P):
            return np.array([np.imag(P[int(size_full/2)-1, int(size_full/2)]) * 2 * e / hbar * (H[int(size_full/2)-1, int(size_full/2)]) for n in range(size_full - 1)])
        
        def Get_J_Leads_Only_Site(P):
            Ps = np.linalg.solve(Ub, P @ np.linalg.solve(Ubd, np.eye(U.shape[0])))
            Ps_EML = Ps[size_lead_L:size_lead_L + size_EM, :size_lead_L]
            Ps_EMR = Ps[size_lead_L:size_lead_L + size_EM, -size_lead_R:]
            return np.imag((Ps_EML @ Hs_LEM - Ps_EMR @ Hs_REM).trace())[0,0] / hbar * e
        
        output_funcs_site_dict = {
            "Full P Matrix":  lambda P: P,
            "P@S diagonal":   lambda P: (P @ S).diagonal(),
            "J All Elements": Get_J_All_Elements_Site,
            "J Leads Only":   Get_J_Leads_Only_Site,
            "J Central":      Get_J_Central_Site,
        }    
        self.output_funcs_site_dict=output_funcs_site_dict
        
        
        
        ####################################################
        ## Define derivative state
        
        def dPdt_State_Scaled(Pss):
            COMMUTE = -1j * (Hss_scaled @ Pss - Pss @ Hss_scaled)

            if Gamma == 0:
                DAMP = 0
            else:
                DAMP = -Gamma_scaled * (Pss * (CORNERS + CROSS / 2) - Target_P0ss * CONTACTS) 

            return COMMUTE + DAMP

        self.dPdt_State_Scaled=dPdt_State_Scaled
        
        ####################################################
        ## Define outputs state
        
        def Get_J_Leads_Only_State(Pss):
            Ps=U@Pss@Ud
            Ps_EML = Ps[self.size_lead_L:self.size_lead_L + self.size_EM, :self.size_lead_L]
            Ps_EMR = Ps[self.size_lead_L:self.size_lead_L + self.size_EM, -self.size_lead_R:]
            return np.imag((Ps_EML @ Hs_LEM - Ps_EMR @ Hs_REM).trace())[0,0] / hbar * e
        
        UdUbdS=UdUbd@ S
        
        def Get_J_Central_State(Pss):
            P=(UbU@Pss@UdUbd)
            return np.array([np.imag(P[int(size_full/2)-1, int(size_full/2)]) * 2 * e / hbar * (H[int(size_full/2)-1, int(size_full/2)]) for n in range(size_full - 1)])

        output_funcs_state_dict = {
            "Full Pss Matrix"         : lambda Pss: Pss,
            "Pss diagonal"            : lambda Pss: Pss.diagonal(),
            "P@S diagonal From State" : lambda Pss: (UbU@Pss@UdUbdS).diagonal(),
            "J Leads Only From State" : Get_J_Leads_Only_State,
            "J Central From State"    : Get_J_Central_State,
        }    
        self.output_funcs_state_dict=output_funcs_state_dict
        
        
        
        
        
        
    ########################################################################################################

#     def Propagate_In_Site(self, solver_type="Forwards Difference"):
        
# #         self.Initial_P0
# #         self.output_funcs_site_dict
# #         self.dPdt_Site_Scaled
# #         self.Delta_t_scaled
# #         self.max_iters
# #         self.list_of_outputs
        
#         ####################################################
#         ## Actually loop
        
#         start = time.time()
        
#         P = self.Initial_P0
#         previous_step_info = {}
#         outputs_dict = {output_type: [self.output_funcs_site_dict[output_type](P)] for output_type in self.list_of_outputs}
        
#         for i in range(1, self.max_iters):
#             P, previous_step_info = Solver_Func_Dict[solver_type](P, self.dPdt_Site_Scaled, previous_step_info, i, self.Delta_t_scaled)
            
#             for output_type in self.list_of_outputs:
#                 outputs_dict[output_type].append(self.output_funcs_site_dict[output_type](P))
                
#             if i % 10 == 0:
#                 print(f"Step {i} out of {self.max_iters} - Time: {time.time() - start:.2f}s")
                
#         self.outputs_dict=outputs_dict
#         self.t_vec_in_fs=np.array(range(0,self.max_iters))*self.Delta_t/fs
#         return outputs_dict
    
    
    
#     ########################################################################################################

#     def Propagate_In_State(self, solver_type="Forwards Difference"):
        
# #         self.Initial_P0ss
# #         self.output_funcs_state_dict
# #         self.dPdt_State_Scaled
# #         self.Delta_t_scaled
# #         self.max_iters
# #         self.list_of_outputs
        
#         ####################################################
#         ## Actually loop
        
#         start = time.time()
        
#         Pss = self.Initial_P0ss
#         previous_step_info = {}
#         outputs_dict = {output_type: [self.output_funcs_state_dict[output_type](Pss)] for output_type in self.list_of_outputs}
        
#         for i in range(1, self.max_iters):
#             Pss, previous_step_info = Solver_Func_Dict[solver_type](Pss, self.dPdt_State_Scaled, previous_step_info, i, self.Delta_t_scaled)
            
#             for output_type in self.list_of_outputs:
#                 outputs_dict[output_type].append(self.output_funcs_state_dict[output_type](Pss))
                
#             if i % 10 == 0:
#                 print(f"Step {i} out of {self.max_iters} - Time: {time.time() - start:.2f}s")
                
#         self.outputs_dict=outputs_dict
#         self.t_vec_in_fs=np.array(range(0,self.max_iters))*self.Delta_t/fs
#         return outputs_dict
    
    
    
    
    
    ########################################################################################################
    
    def Propagate_In_Site(self, solver_type="RK45"):
        """
        Propagates the system using the selected solver with adaptive time-stepping for the Site loop.
        Updates outputs_dict only if the time step is successful.
        """
        start = time.time()

        P = self.Initial_P0
        previous_step_info = {"converged":True}
        Delta_t_scaled = self.Delta_t_scaled
        NBas = P.size  # Total number of basis functions (dimension of the matrices)
        outputs_dict = {output_type: [self.output_funcs_site_dict[output_type](P)] for output_type in self.list_of_outputs}
        t_vec_scaled = [0]

        for i in range(1, self.max_iters):
            # Call solver with updated previous_step_info which includes convergence status
            P, previous_step_info, Delta_t_scaled = Solver_Func_Dict[solver_type](P, self.dPdt_Site_Scaled, previous_step_info, i, Delta_t_scaled, NBas)

            # Convergence check from previous_step_info
            if previous_step_info["converged"]:
                # Convergence succeeded, update t_vec_scaled and append to outputs_dict
                t_vec_scaled.append(t_vec_scaled[-1] + Delta_t_scaled)
                for output_type in self.list_of_outputs:
                    outputs_dict[output_type].append(self.output_funcs_site_dict[output_type](P))
            else:
                print(f"Step {i}: Convergence failed, restoring P from the previous time step.")

            # Output progress every 10 steps
            if i % 10 == 0:
                print(f"Step {i} out of {self.max_iters} - Time: {time.time() - start:.2f}s")
                print("Delta_t, current_t [fs]=",Delta_t_scaled* self.t_scale / fs,t_vec_scaled[-1]* self.t_scale / fs)
                
            if not self.max_t is None:
                if t_vec_scaled[-1]>(self.max_t/self.t_scale):
                    break

        self.outputs_dict = outputs_dict
        t_vec_scaled = np.array(t_vec_scaled)
        self.t_vec_scaled = t_vec_scaled
        self.t_vec_in_fs = t_vec_scaled * self.t_scale / fs
        return outputs_dict



    def Propagate_In_State(self, solver_type="RK45"):
        """
        Propagates the system using the selected solver with adaptive time-stepping for the State loop.
        Updates outputs_dict only if the time step is successful.
        """
        start = time.time()

        Pss = self.Initial_P0ss
        previous_step_info = {"converged":True}
        Delta_t_scaled = self.Delta_t_scaled
        NBas = Pss.size  # Total number of basis functions (dimension of the matrices)
        outputs_dict = {output_type: [self.output_funcs_state_dict[output_type](Pss)] for output_type in self.list_of_outputs}
        t_vec_scaled = [0]

        for i in range(1, self.max_iters):
            # Call solver with updated previous_step_info which includes convergence status
            Pss, previous_step_info, Delta_t_scaled = Solver_Func_Dict[solver_type](Pss, self.dPdt_State_Scaled, previous_step_info, i, Delta_t_scaled, NBas)

            # Convergence check from previous_step_info
            if previous_step_info["converged"]:
                # Convergence succeeded, update t_vec_scaled and append to outputs_dict
                t_vec_scaled.append(t_vec_scaled[-1] + Delta_t_scaled)
                for output_type in self.list_of_outputs:
                    outputs_dict[output_type].append(self.output_funcs_state_dict[output_type](Pss))
            else:
                print(f"Step {i}: Convergence failed, restoring P from the previous time step.")

            # Output progress every 10 steps
            if i % 10 == 0:
                print(f"Step {i} out of {self.max_iters} - Time: {time.time() - start:.2f}s")
                print("Delta_t, current_t [fs]=",Delta_t_scaled* self.t_scale / fs,t_vec_scaled[-1]* self.t_scale / fs)
                
            if not self.max_t is None:
                if t_vec_scaled[-1]>(self.max_t/self.t_scale):
                    break

        self.outputs_dict = outputs_dict
        t_vec_scaled = np.array(t_vec_scaled)
        self.t_vec_scaled = t_vec_scaled
        self.t_vec_in_fs = t_vec_scaled * self.t_scale / fs
        return outputs_dict



    
    
    
    
#     from scipy.integrate import solve_ivp
#     from scipy.integrate import solve_ivp
#     import numpy as np

#     def Propagate_In_Site_SciPy(self):
#         """
#         Propagates the system using SciPy's adaptive ODE solver (solve_ivp) for the Site loop.
#         Handles complex matrices by solving real and imaginary parts separately.
#         Stores only extracted outputs to save memory.
#         """
#         start = time.time()

#         P = np.array(self.Initial_P0)  # Convert np.matrix to np.array
#         Delta_t_scaled = self.Delta_t_scaled
#         P_shape = P.shape
#         P_real = np.array(P.real).flatten()
#         P_imag = np.array(P.imag).flatten()

#         T_final = self.max_t / self.t_scale if self.max_t is not None else self.max_iters * Delta_t_scaled

#         def dPdt_real_imag(t, P_flat_real_imag):
#             """ Computes dP/dt split into real and imaginary parts. """
#             P_real, P_imag = np.split(P_flat_real_imag, 2)
#             P_mat = (P_real + 1j * P_imag).reshape(P_shape)
#             dPdt = np.array(self.dPdt_Site_Scaled(P_mat))  # Convert dPdt to np.array
#             return np.concatenate([dPdt.real.flatten(), dPdt.imag.flatten()])

#         # Output storage
#         outputs_dict = {output_type: [] for output_type in self.list_of_outputs}
#         t_vec_scaled = []

#         def store_outputs(t, P_flat_real_imag):
#             """ Callback to store extracted outputs. """
#             P_real, P_imag = np.split(P_flat_real_imag, 2)
#             P_mat = (P_real + 1j * P_imag).reshape(P_shape)
#             t_vec_scaled.append(t)
#             for output_type in self.list_of_outputs:
#                 outputs_dict[output_type].append(self.output_funcs_site_dict[output_type](P_mat))

#         # Solve real-valued ODE
#         P0_real_imag = np.concatenate([P_real, P_imag])  # Ensure 1D
#         sol = solve_ivp(dPdt_real_imag, (0, T_final), P0_real_imag, method="LSODA",  # Change method to LSODA for stability
#                         atol=1e-8, rtol=1e-6, vectorized=False, dense_output=False)

#         for i in range(len(sol.t)):  
#             store_outputs(sol.t[i], sol.y[:, i])
#             if i % 10 == 0:
#                 print(f"Step {i} out of {self.max_iters} - Time: {time.time() - start:.2f}s")
#                 print("Delta_t, current_t [fs]=",Delta_t_scaled* self.t_scale / fs,t_vec_scaled[-1]* self.t_scale / fs)

#         self.outputs_dict = outputs_dict
#         self.t_vec_scaled = np.array(t_vec_scaled)
#         self.t_vec_in_fs = self.t_vec_scaled * self.t_scale / fs

#         return outputs_dict


#     def Propagate_In_State_SciPy(self):
#         """
#         Propagates the system using SciPy's adaptive ODE solver (solve_ivp) for the State loop.
#         Handles complex matrices by solving real and imaginary parts separately.
#         Stores only extracted outputs to save memory.
#         """
#         start = time.time()

#         Pss = np.array(self.Initial_P0ss)  # Convert np.matrix to np.array
#         Delta_t_scaled = self.Delta_t_scaled
#         Pss_shape = Pss.shape
#         Pss_real = np.array(Pss.real).flatten()
#         Pss_imag = np.array(Pss.imag).flatten()

#         T_final = self.max_t / self.t_scale if self.max_t is not None else self.max_iters * Delta_t_scaled

#         def dPdt_real_imag(t, Pss_flat_real_imag):
#             """ Computes dP/dt split into real and imaginary parts. """
#             Pss_real, Pss_imag = np.split(Pss_flat_real_imag, 2)
#             Pss_mat = (Pss_real + 1j * Pss_imag).reshape(Pss_shape)
#             dPdt = np.array(self.dPdt_State_Scaled(Pss_mat))  # Convert dPdt to np.array
#             return np.concatenate([dPdt.real.flatten(), dPdt.imag.flatten()])

#         # Output storage
#         outputs_dict = {output_type: [] for output_type in self.list_of_outputs}
#         t_vec_scaled = []

#         def store_outputs(t, Pss_flat_real_imag):
#             """ Callback to store extracted outputs. """
#             Pss_real, Pss_imag = np.split(Pss_flat_real_imag, 2)
#             Pss_mat = (Pss_real + 1j * Pss_imag).reshape(Pss_shape)
#             t_vec_scaled.append(t)
#             for output_type in self.list_of_outputs:
#                 outputs_dict[output_type].append(self.output_funcs_state_dict[output_type](Pss_mat))

#         # Solve real-valued ODE
#         Pss0_real_imag = np.concatenate([Pss_real, Pss_imag])  # Ensure 1D
#         sol = solve_ivp(dPdt_real_imag, (0, T_final), Pss0_real_imag, method="LSODA",  # Change method to LSODA for stability
#                         atol=1e-8, rtol=1e-6, vectorized=False, dense_output=False)

#         for i in range(len(sol.t)):  
#             store_outputs(sol.t[i], sol.y[:, i])
#             if i % 10 == 0:
#                 print(f"Step {i} out of {self.max_iters} - Time: {time.time() - start:.2f}s")
#                 print("Delta_t, current_t [fs]=",Delta_t_scaled* self.t_scale / fs,t_vec_scaled[-1]* self.t_scale / fs)

#         self.outputs_dict = outputs_dict
#         self.t_vec_scaled = np.array(t_vec_scaled)
#         self.t_vec_in_fs = self.t_vec_scaled * self.t_scale / fs

#         return outputs_dict


    
    
    ########################################################################################################
    
    def Visualize_Results(self,figsize=(5,5),dpi=200):
        plots_dict={}

        outputs_dict=self.outputs_dict
        t_vec_in_fs=self.t_vec_in_fs
        size_lead_L=self.size_lead_L
        size_lead_R=self.size_lead_R
        size_EM=self.size_EM

        for output_type in outputs_dict.keys():

            if output_type=="J All Elements":
                plots_dict[output_type]=  Plot_J_All_Elements(  outputs_dict[output_type],t_vec_in_fs,size_lead_L,size_lead_R,size_EM,figsize,dpi)

            elif output_type=="J Leads Only" or output_type=="J Leads Only From State":
                plots_dict[output_type]=  Plot_J_Leads_Only(    outputs_dict[output_type],t_vec_in_fs,size_lead_L,size_lead_R,size_EM,figsize,dpi)

            elif output_type=="P@S diagonal" or output_type=="P@S diagonal From State":
                plots_dict[output_type]=  Plot_PS_Diag(         outputs_dict[output_type],t_vec_in_fs,size_lead_L,size_lead_R,size_EM,figsize,dpi)

            elif output_type=="Pss diagonal":
                plots_dict[output_type]=  Plot_Pss_Diag(        outputs_dict[output_type],t_vec_in_fs,size_lead_L,size_lead_R,size_EM,figsize,dpi)
                
            else:
                plots_dict[output_type]=()

        self.plots_dict=plots_dict
        return plots_dict




########################################################################################################
## Results Plotting

def Plot_J_All_Elements(J_vec,t_vec,size_lead_L,size_lead_R,size_EM,figsize,dpi):

    elements_to_plot_dict={}

    n=0
    elements_to_plot_dict["L Lead End"+" ({0}-{1})".format(n,n+1)]=n

    n=int((size_lead_L-1)/2)
    elements_to_plot_dict["L Lead Mid"+" ({0}-{1})".format(n,n+1)]=n

    n=size_lead_L-1
    elements_to_plot_dict["L Lead - EM"+" ({0}-{1})".format(n,n+1)]=n

    n=size_lead_L+int((size_EM-1)/2)
    elements_to_plot_dict["EM Mid"+" ({0}-{1})".format(n,n+1)]=n

    n=size_lead_L+size_EM-1
    elements_to_plot_dict["R Lead - EM"+" ({0}-{1})".format(n,n+1)]=n

    n=size_lead_L+size_EM+int((size_lead_R-1)/2)
    elements_to_plot_dict["R Lead Mid"+" ({0}-{1})".format(n,n+1)]=n

    n=size_lead_L+size_lead_R+size_EM-2
    elements_to_plot_dict["R Lead End"+" ({0}-{1})".format(n,n+1)]=n


    fig1=plt.figure(figsize=figsize,dpi=dpi,facecolor="white")
    colors=plt.cm.twilight(np.linspace(0,1,len(J_vec[0])))
    for n in range(len(J_vec[0])):
        plt.plot(t_vec,[J_vec[i][n]*1e3 for i in range(len(J_vec))],label="_",color=colors[n])
    for label,n in elements_to_plot_dict.items():
        plt.plot([],[],label=label,color=colors[n])
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title("Device Current (All Nodes)")
    plt.xlabel("Time [fs]")
    plt.ylabel("Current [mA]")
    # plt.ylim(-0.03,0.03)
    plt.show()
    plt.close()


    fig2 = plt.figure(figsize=figsize,dpi=dpi,facecolor="white")
    colors = plt.cm.plasma(np.linspace(0, 1, len(elements_to_plot_dict.keys())))
    base_linewidth = 2  # Starting linewidth
    linewidth_step = 0.2  # Decrease each line slightly
    for color_ind, (label, n) in enumerate(elements_to_plot_dict.items()):
        plt.plot(
            t_vec, 
            [J_vec[i][n] * 1e3 for i in range(len(J_vec))], 
            label=label, 
            color=colors[color_ind], 
            linewidth=max(0.5, base_linewidth - color_ind * linewidth_step)  # Ensure linewidth doesn't go too low
        )
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title("Device Current (All Nodes)")
    plt.xlabel("Time [fs]")
    plt.ylabel("Current [mA]")
    # plt.ylim(-0.03,0.03)
    plt.show()
    plt.close()

    fig3=plt.figure(figsize=figsize,dpi=dpi,facecolor="white")
    plt.axvline(size_lead_L-1,color="k")
    plt.axvline(size_lead_L-1+size_EM,color="k")
    plt.plot(np.arange(len(J_vec[0]))+0.5,np.array([max([J_vec[i][n] for i in range(len(J_vec))]) for n in range(len(J_vec[0]))])*1e3)
    plt.plot(np.arange(len(J_vec[0]))+0.5,np.array([min([J_vec[i][n] for i in range(len(J_vec))]) for n in range(len(J_vec[0]))])*1e3)
    plt.title("Device Current (All Nodes)")
    plt.xlabel("Element")
    plt.ylabel("Max/Min Current [mA]")
    # plt.ylim(-0.03,0.03)
    plt.show()
    plt.close()
    
    return (fig1,fig2,fig3)


def Plot_J_Leads_Only(J_vec,t_vec,size_lead_L,size_lead_R,size_EM,figsize,dpi):

    fig=plt.figure(figsize=figsize,dpi=dpi,facecolor="white")
    plt.plot(t_vec,np.array(J_vec)*1e3)
    plt.title("Device Current (LEM-REM)")  
    plt.xlabel("Time [fs]")
    plt.ylabel("Current [mA]")
    plt.show()
    plt.close()

    return (fig,)

def Plot_PS_Diag(PS_vec,t_vec,size_lead_L,size_lead_R,size_EM,figsize,dpi):

    # Plotting absolute diagonals of P
    range_thing=range(size_lead_L,size_lead_L+size_EM)
    colors=plt.cm.plasma(np.linspace(0,1,len(range_thing)))
    fig1=plt.figure(figsize=figsize,dpi=dpi,facecolor="white")
    for color_ind,n in enumerate(range_thing):
        plt.plot(t_vec,[abs(PS[0,n]) for PS in PS_vec], label="_",color=colors[color_ind])
    for color_ind,n in enumerate(range_thing):
        if n in [int(n) for n in np.linspace(range_thing[0],range_thing[-1],10)]:
            plt.plot([],[],label="State "+str(int(n)),color=colors[int(color_ind)])
    plt.plot(t_vec,[abs(sum([PS[0,n] for n in range_thing])) for PS in PS_vec], label='sum',color="g")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title("Abs Diagonals of P@S EM")
    plt.axhline([1],color="k")
    plt.xlabel("Time [fs]")
    plt.ylim(0,2)
    plt.show()
    plt.close()
    
    range_thing=range(size_lead_L)
    colors=plt.cm.plasma(np.linspace(0,1,len(range_thing)))
    fig2=plt.figure(figsize=figsize,dpi=dpi,facecolor="white")
    for color_ind,n in enumerate(range_thing):
        plt.plot(t_vec,[abs(PS[0,n]) for PS in PS_vec], label="_",color=colors[color_ind])
    for color_ind,n in enumerate(range_thing):
        if n in [int(n) for n in np.linspace(range_thing[0],range_thing[-1],10)]:
            plt.plot([],[],label="State "+str(int(n)),color=colors[int(color_ind)])
    plt.plot(t_vec,[abs(sum([PS[0,n] for n in range_thing])) for PS in PS_vec], label='sum',color="g")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title("Abs Diagonals of P@S L Lead")
    plt.axhline([1],color="k")
    plt.xlabel("Time [fs]")
    plt.ylim(0,2)
    plt.show()
    plt.close()

    range_thing=range(size_lead_L+size_EM,size_lead_L+size_EM+size_lead_R)
    colors=plt.cm.plasma(np.linspace(0,1,len(range_thing)))
    fig3=plt.figure(figsize=figsize,dpi=dpi,facecolor="white")
    for color_ind,n in enumerate(range_thing):
        plt.plot(t_vec,[abs(PS[0,n]) for PS in PS_vec], label="_",color=colors[color_ind])
    for color_ind,n in enumerate(range_thing):
        if n in [int(n) for n in np.linspace(range_thing[0],range_thing[-1],10)]:
            plt.plot([],[],label="State "+str(int(n)),color=colors[int(color_ind)])
    plt.plot(t_vec,[abs(sum([PS[0,n] for n in range_thing])) for PS in PS_vec], label='sum',color="g")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title("Abs Diagonals of P@S R Lead")
    plt.axhline([1],color="k")
    plt.xlabel("Time [fs]")
    plt.ylim(0,2)
    plt.show()
    plt.close()

    return (fig1,fig2,fig3)

def Plot_Pss_Diag(Pss_vec,t_vec,size_lead_L,size_lead_R,size_EM,figsize,dpi):

    fig=plt.figure(figsize=figsize,dpi=dpi,facecolor="white")
    colors=plt.cm.RdBu(np.linspace(0,1,len(t_vec)))
    for t_ind in range(len(t_vec)):
        plt.plot(abs(np.array(Pss_vec[t_ind])[0].flatten()), label="_",color=colors[t_ind])
    for t_ind in [int(t) for t in np.linspace(0,len(t_vec)-1,10)]:
        plt.plot([],[], label="time_ind="+str(t_ind),color=colors[t_ind])
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title("State Occupancy")  
    plt.xlabel("Element")
    plt.ylabel("Occupancy")
    plt.ylim(-0.1,1.1)
    plt.show()
    plt.close()

    return (fig,)