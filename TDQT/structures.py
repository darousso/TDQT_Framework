from .common import *

def Get_Structure(structure="Oded",args_to_pass={}):
    
    ######################################################################################################
    if structure=="Oded":
    
        SizeLeadL=300
        SizeLeadR=300

        SizeML=50
        SizeMR=50
        SizeM=6

        N_s=SizeML+SizeM+SizeMR

        alpha_L=0
        alpha_M=0
        alpha_R=0

        beta_L=-0.2
        beta_M=-0.2
        beta_R=-0.2
        beta_LM=-0.2
        beta_RM=-0.2
        
        
        ###################################################
        
        G_EM_LR=np.zeros((SizeML,SizeMR))
        G_EM_L=matrix_block(alpha_L,beta_L,SizeML)
        G_EM_R=matrix_block(alpha_R,beta_R,SizeMR)
        G_EM_LEM=np.zeros((SizeML,SizeM))
        G_EM_LEM[-1,0]=beta_LM
        G_EM_REM=np.zeros((SizeML,SizeM))
        G_EM_REM[0,-1]=beta_RM


        G_EM_RL=G_EM_LR.T
        G_EM_EML=G_EM_LEM.T
        G_EM_EMR=G_EM_REM.T

        G_EM_M=matrix_block(alpha_M,beta_M,SizeM)


        G_EM=np.block([[G_EM_L,G_EM_LEM,G_EM_LR],[G_EM_EML,G_EM_M,G_EM_EMR],[G_EM_RL,G_EM_REM,G_EM_R]])


        ########################################

        G_LR=np.zeros((SizeLeadL,SizeLeadR))
        G_L=alpha_L*np.diag(np.ones(SizeLeadL),k=0)+beta_L*np.diag(np.ones(SizeLeadL-1),k=1)+beta_L*np.diag(np.ones(SizeLeadL-1),k=-1)
        G_R=alpha_R*np.diag(np.ones(SizeLeadR),k=0)+beta_R*np.diag(np.ones(SizeLeadR-1),k=1)+beta_R*np.diag(np.ones(SizeLeadR-1),k=-1)
        G_LEM=np.zeros((SizeLeadL,N_s))
        G_LEM[-1,0]=beta_L
        G_REM=np.zeros((SizeLeadR,N_s))
        G_REM[0,-1]=beta_R


        G_RL=G_LR.T
        G_EML=G_LEM.T
        G_EMR=G_REM.T

        ########################################

        G=np.block([[G_L,G_LEM,G_LR],[G_EML,G_EM,G_EMR],[G_RL,G_REM,G_R]])

        ########################################
        
        eigenvalues, eigenvectors = scipy.linalg.eigh(G)  # Sorted eigenvalues

        S = eigenvectors.T @ eigenvectors  # Should be identity if eigenvectors are orthonormal
        
        
    ######################################################################################################
    elif structure=="TinyH":
        
        with open(os.path.join(os.path.dirname(__file__),"..","StructuresLibrary","TinyH","H.txt"),"r") as f:
            HINFORAW=f.read()
        HINFO=[[np.complex128(i) for i in thing.split(' ')] for thing in HINFORAW.split('\n') if thing!=""]
        SizeM=int(max(max([i[0] for i in HINFO]),max([i[1] for i in HINFO])))+1
        G_EM_M=np.zeros((SizeM,SizeM))
        for line in HINFO:
            G_EM_M[int(line[0]),int(line[1])]=line[2]
            G_EM_M[int(line[1]),int(line[0])]=line[2]
        G_EM_M=G_EM_M*hartree_to_eV +0j
    
        ########################################
        
        SizeLeadL=300
        SizeLeadR=300

        SizeML=50
        SizeMR=50
        
        N_s=SizeML+SizeM+SizeMR

        alpha_L=0
        alpha_M=0
        alpha_R=0

        beta_L=-0.2
        beta_M=-0.2
        beta_R=-0.2
        beta_LM=-0.2
        beta_RM=-0.2


        ###################################################

        G_EM_LR=np.zeros((SizeML,SizeMR))
        G_EM_L=matrix_block(alpha_L,beta_L,SizeML)
        G_EM_R=matrix_block(alpha_R,beta_R,SizeMR)
        G_EM_LEM=np.zeros((SizeML,SizeM))
        G_EM_LEM[-1,0]=beta_LM
        G_EM_REM=np.zeros((SizeML,SizeM))
        G_EM_REM[0,-1]=beta_RM


        G_EM_RL=G_EM_LR.T
        G_EM_EML=G_EM_LEM.T
        G_EM_EMR=G_EM_REM.T

        
        G_EM=np.block([[G_EM_L,G_EM_LEM,G_EM_LR],[G_EM_EML,G_EM_M,G_EM_EMR],[G_EM_RL,G_EM_REM,G_EM_R]])


        ########################################

        G_LR=np.zeros((SizeLeadL,SizeLeadR))
        G_L=alpha_L*np.diag(np.ones(SizeLeadL),k=0)+beta_L*np.diag(np.ones(SizeLeadL-1),k=1)+beta_L*np.diag(np.ones(SizeLeadL-1),k=-1)
        G_R=alpha_R*np.diag(np.ones(SizeLeadR),k=0)+beta_R*np.diag(np.ones(SizeLeadR-1),k=1)+beta_R*np.diag(np.ones(SizeLeadR-1),k=-1)
        G_LEM=np.zeros((SizeLeadL,N_s))
        G_LEM[-1,0]=beta_L
        G_REM=np.zeros((SizeLeadR,N_s))
        G_REM[0,-1]=beta_R


        G_RL=G_LR.T
        G_EML=G_LEM.T
        G_EMR=G_REM.T

        ########################################

        G=np.block([[G_L,G_LEM,G_LR],[G_EML,G_EM,G_EMR],[G_RL,G_REM,G_R]])

        ########################################

        eigenvalues, eigenvectors = scipy.linalg.eigh(G)  # Sorted eigenvalues

        S = eigenvectors.T @ eigenvectors  # Should be identity if eigenvectors are orthonormal

        
    ######################################################################################################
    else:
        print("structure=",structure," not recognized as an option")
        return None
    
    return np.matrix(G),np.matrix(S),SizeLeadL,SizeLeadR


def matrix_block(alpha,beta,Size):
    return alpha*np.diag(np.ones(Size),k=0)+beta*np.diag(np.ones(Size-1),k=1)+beta*np.diag(np.ones(Size-1),k=-1)  