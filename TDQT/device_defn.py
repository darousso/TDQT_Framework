from .common import *
from .utils import *

class Device:
    def __init__(self,
                H, S, 
                size_lead_L,size_lead_R,
                V_bias,
                T_L=0, T_R=0, T_EM=0,
                clean_matrices=False,
            ):
        
        """Initializes the device and the relevant matrices"""
        
        ####################################################
        
        size_full=H.shape[0]
        size_EM=size_full-size_lead_L-size_lead_R
        
        if clean_matrices:
            H=Clean_Matrix(H)
            S=Clean_Matrix(S)
        
        ####################################################
        
        S_L=S[:size_lead_L,:size_lead_L]
        S_EM=S[size_lead_L:size_lead_L+size_EM,size_lead_L:size_lead_L+size_EM]
        S_R=S[size_lead_L+size_EM:,size_lead_L+size_EM:]
        S_LEM=S[:size_lead_L,size_lead_L:size_lead_L+size_EM]
        S_EML=S[size_lead_L:size_lead_L+size_EM,:size_lead_L]
        S_REM=S[size_lead_L+size_EM:,size_lead_L:size_lead_L+size_EM]
        S_EMR=S[size_lead_L:size_lead_L+size_EM,size_lead_L+size_EM:]
        S_LR=S[:size_lead_L,size_lead_L+size_EM:]
        S_RL=S[size_lead_L+size_EM:,:size_lead_L]
        
        ####################################################
        ## Ub
        
        invS_L=np.linalg.inv(S_L)
        invS_R=np.linalg.inv(S_R)

        Ub=np.matrix(np.eye(size_full))+0j
        Ub[:size_lead_L,size_lead_L:size_lead_L+size_EM]=-invS_L@S_LEM
        Ub[-size_lead_R:,size_lead_L:size_lead_L+size_EM]=-invS_R@S_REM
        
        if clean_matrices:
            Ub=Clean_Matrix(Ub)
        
        ####################################################
        ## Ss, Hs

        Ss=Ub.H@S@Ub
        Hs=Ub.H@H@Ub
        
        if clean_matrices:
            Ss=Clean_Matrix(Ss)
            Hs=Clean_Matrix(Hs)

        Ss_L=Ss[:size_lead_L,:size_lead_L]
        Ss_EM=Ss[size_lead_L:size_lead_L+size_EM,size_lead_L:size_lead_L+size_EM]
        Ss_R=Ss[size_lead_L+size_EM:,size_lead_L+size_EM:]

        Hs_L=Hs[:size_lead_L,:size_lead_L]
        Hs_EM=Hs[size_lead_L:size_lead_L+size_EM,size_lead_L:size_lead_L+size_EM]
        Hs_R=Hs[size_lead_L+size_EM:,size_lead_L+size_EM:]
        
        Hs_LEM=Hs[:size_lead_L,size_lead_L:size_lead_L+size_EM]
        Hs_REM=Hs[-size_lead_R:,size_lead_L:size_lead_L+size_EM]

        ####################################################
        ## U

        U=np.matrix(np.zeros((size_full,size_full)))+0j
        
        E_EM,U_EM=scipy.linalg.eigh(Hs_EM,b=Ss_EM)
        U[size_lead_L:size_lead_L+size_EM,size_lead_L:size_lead_L+size_EM]=U_EM
        U_EM=np.matrix(U_EM)

        E_L,U_L=scipy.linalg.eigh(Hs_L,b=Ss_L)
        U[:size_lead_L,:size_lead_L]=U_L
        U_L=np.matrix(U_L)

        E_R,U_R=scipy.linalg.eigh(Hs_R,b=Ss_R)
        U[-size_lead_R:,-size_lead_R:]=U_R
        U_R=np.matrix(U_R)

        U=np.matrix(U)
        
        if clean_matrices:
            U=Clean_Matrix(U)

        ####################################################
        ## Hss

        Hss=U.H@Hs@U

        if clean_matrices:
            Hss=Clean_Matrix(Hss)
            
        ####################################################  
        ## initial and target density matrix

        E_F_L=np.mean([sorted(E_L)[int(len(E_L)/2)-1],sorted(E_L)[int(len(E_L)/2)]])
        E_F_R=np.mean([sorted(E_R)[int(len(E_R)/2)-1],sorted(E_R)[int(len(E_R)/2)]])
        E_F_EM=np.mean([sorted(E_EM)[int(len(E_EM)/2)-1],sorted(E_EM)[int(len(E_EM)/2)]])
        E_F=np.real(np.mean([sorted(np.diag(Hss))[int(len(np.diag(Hss))/2)-1],sorted(np.diag(Hss))[int(len(np.diag(Hss))/2)]]))

        mu_L=E_F_R+V_bias/2 #eV
        mu_R=E_F_L-V_bias/2 #eV
        mu_EM=E_F_EM
        Target_P0ss=np.zeros((size_full,size_full))*0j
        Target_P0ss[:size_lead_L,:size_lead_L]=np.diag(expit(- (E_L - mu_L) / (kB * T_L)))
        Target_P0ss[-size_lead_R:,-size_lead_R:]=np.diag(expit(- (E_R - mu_R) / (kB * T_R)))
        Target_P0ss[size_lead_L:size_lead_L+size_EM,size_lead_L:size_lead_L+size_EM]=np.diag(expit(- (E_EM - mu_EM) / (kB * T_EM)))
        Target_P0s=U@Target_P0ss@U.H
        Target_P0=Ub@Target_P0s@Ub.H

#         mu_L=E_F
#         mu_R=E_F
#         mu_EM=E_F
#         Initial_P0ss=np.diag(expit(- (np.real(np.diag(Gss)) - E_F) / (kB * T_EM)))
#         Initial_P0s=U@Initial_P0ss@Ud
#         Initial_P0=Ub@Initial_P0s@Ubd

        Initial_P0ss=Target_P0ss
        Initial_P0s=Target_P0s
        Initial_P0=Target_P0
        
        if clean_matrices:
            Target_P0ss  = Clean_Matrix(Target_P0ss )
            Target_P0s   = Clean_Matrix(Target_P0s  )
            Target_P0    = Clean_Matrix(Target_P0   )
            Initial_P0ss = Clean_Matrix(Initial_P0ss)
            Initial_P0s  = Clean_Matrix(Initial_P0s )
            Initial_P0   = Clean_Matrix(Initial_P0  )
        
        ####################################################
        ## Assigning to self
        
        self.H=H
        self.S=S
        
        self.Ub=Ub
        self.Hs=Hs
        self.Ss=Ss
        
        self.U=U
        self.Hss=Hss
        
        self.Target_P0ss=Target_P0ss
        self.Target_P0s=Target_P0s
        self.Target_P0=Target_P0
        self.Initial_P0ss=Initial_P0ss
        self.Initial_P0s=Initial_P0s
        self.Initial_P0=Initial_P0
        
        self.E_F_L=E_F_L
        self.E_F_EM=E_F_EM
        self.E_F_R=E_F_R
        self.E_F=E_F
        
        self.Hs_REM=Hs_REM
        self.Hs_LEM=Hs_LEM
        
        self.size_full=size_full
        self.size_EM=size_EM
        self.size_lead_L=size_lead_L
        self.size_lead_R=size_lead_R
        
        self.clean_matrices=clean_matrices
        
        ####################################################
        
    def Transform_State_To_Site_Density(self, Pss):
        """Transforms a given density matrix from the state to the site representation."""
        return self.Ub@self.U@Pss@self.U.H@self.Ub.H

    def Transform_Site_To_Space_Density(self, P):
        """Transforms a given density matrix from the site to the state representation."""
        return np.linalg.solve(self.Ub@self.U, P @ np.linalg.solve(self.Ub.H@self.U.H, np.eye(self.U.shape[0])))

    def Visualize(self):
        """Visualizes the state of the device."""
        
        Plot_Matrix(self.H            ,"H"            ,size_lead_L=self.size_lead_L,size_lead_R=self.size_lead_R)
        Plot_Matrix(self.S            ,"S"            ,size_lead_L=self.size_lead_L,size_lead_R=self.size_lead_R)
        
        Plot_Matrix(self.Ub           ,"Ub"           ,size_lead_L=self.size_lead_L,size_lead_R=self.size_lead_R)
        Plot_Matrix(self.Hs           ,"Hs"           ,size_lead_L=self.size_lead_L,size_lead_R=self.size_lead_R)
        Plot_Matrix(self.Ss           ,"Ss"           ,size_lead_L=self.size_lead_L,size_lead_R=self.size_lead_R)
        
        Plot_Matrix(self.U            ,"U"            ,size_lead_L=self.size_lead_L,size_lead_R=self.size_lead_R)
        Plot_Matrix(self.Hss          ,"Hss"          ,size_lead_L=self.size_lead_L,size_lead_R=self.size_lead_R)

        Plot_Matrix(self.Target_P0ss  ,"Target_P0ss"  ,size_lead_L=self.size_lead_L,size_lead_R=self.size_lead_R)
        Plot_Matrix(self.Target_P0s   ,"Target_P0s"   ,size_lead_L=self.size_lead_L,size_lead_R=self.size_lead_R)
        Plot_Matrix(self.Target_P0    ,"Target_P0"    ,size_lead_L=self.size_lead_L,size_lead_R=self.size_lead_R)
        
        Plot_Matrix(self.Initial_P0ss ,"Initial_P0ss" ,size_lead_L=self.size_lead_L,size_lead_R=self.size_lead_R)
        Plot_Matrix(self.Initial_P0s  ,"Initial_P0s"  ,size_lead_L=self.size_lead_L,size_lead_R=self.size_lead_R)
        Plot_Matrix(self.Initial_P0   ,"Initial_P0"   ,size_lead_L=self.size_lead_L,size_lead_R=self.size_lead_R)
        
        plt.figure(figsize=(5,5),dpi=200)
        plt.plot(np.diag(self.Hss),".")
        plt.axvline(self.size_lead_L-0.5)
        plt.axvline(self.size_lead_L+self.size_EM-0.5)
        plt.axhline(self.E_F)
        plt.xlabel("Element")
        plt.ylabel("Energy (diag(Hss))")
        plt.show()
        plt.close()
        
        print("E_F_L=",self.E_F_L)
        print("E_F_R=",self.E_F_R)
        print("E_F_EM=",self.E_F_EM)
        print("E_F=",self.E_F)