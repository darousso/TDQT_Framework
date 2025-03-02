from .common import *

def Plot_Matrix(G,label="",norm=True,figsize=(10,10),dpi=200,size_lead_L=None,size_lead_R=None):
    if type(G)!=np.ndarray and type(G)!=np.matrix:
        G=G.todense()
        
    if norm:     
        Splot=G.copy()
        Splot=np.abs(Splot)
        Plot_Real_Matrix(Splot,label=label+" Norm",figsize=figsize,dpi=dpi,size_lead_L=size_lead_L,size_lead_R=size_lead_R)
    else:        
        Splot=G.copy()
        Splot=np.real(Splot)
        Plot_Real_Matrix(Splot,label=label+" Real",figsize=figsize,dpi=dpi,size_lead_L=size_lead_L,size_lead_R=size_lead_R)
        
        Splot=G.copy()
        Splot=np.imag(Splot)
        Plot_Real_Matrix(Splot,label=label+" Imag",figsize=figsize,dpi=dpi,size_lead_L=size_lead_L,size_lead_R=size_lead_R)
        
def Plot_Real_Matrix(G,label="",figsize=(10,10),dpi=200,size_lead_L=None,size_lead_R=None):
    plt.figure(figsize=figsize,dpi=dpi)
    plt.title(label)
    Splot=G.copy()
    Splot[Splot==0]=np.nan
    plt.imshow(Splot,interpolation="nearest")
    if not size_lead_L is None and not size_lead_R is None:
        size_EM=G.shape[0]-size_lead_L-size_lead_R
        plt.axhline(size_lead_L-0.5,color="k")
        plt.axhline(size_lead_L+size_EM-0.5,color="k")
        plt.axvline(size_lead_L-0.5,color="k")
        plt.axvline(size_lead_L+size_EM-0.5,color="k")
    plt.colorbar()
    plt.show()
    plt.close()
    

def Clean_Matrix(G, eps=1e-12):
    """
    Cleans the matrix G by setting small real parts to zero if they are below eps
    and small imaginary parts to zero if they are below eps.
    
    Ensures that the matrix remains complex when modified.
    """
    G = G.astype(np.complex128)  # Ensure complex type
    newG = G.copy()

    mask_real = np.abs(np.real(newG)) < eps
    mask_imag = np.abs(np.imag(newG)) < eps

    # Use np.where to avoid shape mismatch
    newG = np.where(mask_real, 1j * np.imag(newG), newG)  # Preserve only the imaginary part
    newG = np.where(mask_imag, np.real(newG), newG)        # Preserve only the real part

    return np.matrix(newG)
