#Gaussian Kernel
import numpy as np
from scipy.special import erf


def gausskernel(x,y,n,sigma):
    """Creates a gauss kernel map of a grid of n rows and n columns
        given a standard deviation sigma. Returns its matrix.
        Computational cost: O(n*n*N)
    """
    N=len(x)
    b=sigma*sigma
    mat=np.zeros((n,n),dtype=float)
    for i in range(n):
        xo=(2*i+1)/(2*n)
        for j in range(n):
            k=0
            yo=(2*j+1)/(2*n)
            while k<N:
                #this r is actually squared
                r=((x[k]-xo)**2+(y[k]-yo)**2)
                mat[i,j]=mat[i,j]+np.exp(-r/b)
                k=k+1
            mat[i,j]=mat[i,j]-(N/4)*np.pi*sigma*sigma*(erf((1-xo)/sigma)+erf(xo/sigma))*(erf((1-yo)/sigma)+erf(yo/sigma))
    return mat


def weightedker(matrix,n,sigma):
    """Creates a gauss kernel map of a grid of n rows and n columns
      given a standard deviation sigma. Returns its matrix.
      Computational cost: O(n*n*N)
    """
    x=len(matrix[:,1])
    y=len(matrix[1,:])
    N=x*y
    b=sigma*sigma
    mat=np.zeros((n,n),dtype=float)
    for i in range(n):
        print i
        xo=(2*i+1)/(2*n)
        for j in range(n):
            k=0
            yo=(2*j+1)/(2*n)
            for k in range(x):
              xi=(2*k+1)/(2*x)
              for l in range(y):
                yi=(2*l+1)/(2*y)
                #this r is actually squared
                r=((xi-xo)**2+(yi-yo)**2)
                mat[i,j]=mat[i,j]+matrix[k,l]*np.exp(-r/b)
            mat[i,j]=mat[i,j]-(N/4)*np.pi*sigma*sigma*(erf((1-xo)/sigma)+erf(xo/sigma))*(erf((1-yo)/sigma)+erf(yo/sigma))
    return mat
    
    
def vecgaussker(x,y,v,theta,n,sigma):
    """Creates a vectorial gauss kernel map of a grid of n rows and n columns
        given a standard deviation sigma, given the vector v and its angle theta.
        Returns a 3-dimensional matrix whose third component is the two components of the vector.
        Computational cost: O(2n*n*N)
    """
    vx=v*np.cos(theta)
    vy=v*np.sin(theta)
    vector=np.matrix([0.2,0.3])
    
    N=len(x)
    mat=np.zeros((n,n,2),dtype=float)
    for i in range(n):
        for j in range(n):
            k=0
            while k<N:
                #this r is actually squared
                r=((x[k]-i/(n-1))**2+(y[k]-j/(n-1))**2)
                if r<0.1:
                    w=np.exp(-r/(2*sigma*sigma))
                    mat[i,j,0]=mat[i,j,0]+w*vx[k]
                    mat[i,j,1]=mat[i,j,1]+w*vy[k]
                k=k+1
                    
            #mat[i,j]=(np.sqrt(mat[i,j][0]**2+mat[i,j][1]**2),arctan2(mat[i,j][1]/mat[i,j][0]))
                
    return mat


def vecgaussker_cartesian(x,y,vx,vy,n,sigma):
    """Creates a vectorial gauss kernel map of a grid of n rows and n columns
        given a standard deviation sigma, given the components of the vector.
        Returns a 3-dimensional matrix whose third component is the two components of the vector.
        Computational cost: O(2n*n*N)
    """
    
    N=len(x)
    mat=np.zeros((n,n,2),dtype=float)
    mat2=np.zeros((n,n),dtype=float)
    #another choice: b=2*sigma*sigma
    b=sigma*sigma
    for i in range(n):
        for j in range(n):
            k=0
            while k<N:
                #this r is actually squared
                r=((x[k]-i/(n-1))**2+(y[k]-j/(n-1))**2)
                w=np.exp(-r/b)
                mat[i,j,0]=mat[i,j,0]+w*vx[k]
                mat[i,j,1]=mat[i,j,1]+w*vy[k]
                #the squared root term can be replaced by 1 if the input is normalized
                mat2[i,j]=mat2[i,j]+w*np.sqrt(vx[k]*vx[k]+vy[k]*vy[k])
                k=k+1
                    
            #mat[i,j]=(np.sqrt(mat[i,j][0]**2+mat[i,j][1]**2),arctan2(mat[i,j][1]/mat[i,j][0]))
    mat[:,:,0]=mat[:,:,0]/mat2
    mat[:,:,1]=mat[:,:,1]/mat2
    return mat


def gaussker_gaia(ra,dec,ravel,decvel,vtang,n,sigma):
    """Creates a vectorial gauss kernel map of a grid of n rows and n columns
        given a standard deviation sigma, given the components of the vector.
        The velocities are treated as unitary vectors.
        Returns a 3-dimensional matrix whose third component is the two components of the vector.
        Computational cost: O(2n*n*N)
    """
    maxra=np.amax(ra)
    minra=np.amin(ra)
    maxdec=np.amax(dec)
    mindec=np.amin(dec)

    sizera=maxra-minra
    sizedec=maxdec-mindec

    x=(ra-minra)/sizera
    y=(dec-mindec)/sizedec

    vx=np.asarray(ravel)/vtang
    vy=np.asarray(decvel)/vtang
    
    N=len(x)
    mat=np.zeros((n,n,2),dtype=float)
    mat2=np.zeros((n,n),dtype=float)
    #another choice: b=2*sigma*sigma
    b=2*sigma*sigma
    for i in range(n):
        for j in range(n):
            k=0
            while k<N:
                #this r is actually squared
                r=((x[k]-i/(n-1))**2+(y[k]-j/(n-1))**2)
                w=np.exp(-r/b)
                mat[i,j,0]=mat[i,j,0]+w*vx[k]
                mat[i,j,1]=mat[i,j,1]+w*vy[k]
                mat2[i,j]=mat2[i,j]+w
                k=k+1
                    
            #mat[i,j]=(np.sqrt(mat[i,j][0]**2+mat[i,j][1]**2),arctan2(mat[i,j][1]/mat[i,j][0]))
    mat[:,:,0]=mat[:,:,0]/mat2
    mat[:,:,1]=mat[:,:,1]/mat2
    return mat

def tri_gaussker_gaia(ra,dec,ravel,decvel,radvel,vtot,n,sigma):
    """Creates a vectorial gauss kernel map of a grid of n rows and n columns
        given a standard deviation sigma, given the components of the vector.
        The velocities are treated as unitary vectors.
        Returns a 3-dimensional matrix whose third component is the three components of the vector.
        Computational cost: O(2n*n*N)
    """
    maxra=np.amax(ra)
    minra=np.amin(ra)
    maxdec=np.amax(dec)
    mindec=np.amin(dec)

    sizera=maxra-minra
    sizedec=maxdec-mindec

    x=(ra-minra)/sizera
    y=(dec-mindec)/sizedec

    vx=np.asarray(ravel)/vtot
    vy=np.asarray(decvel)/vtot
    vz=np.asarray(radvel)/vtot
    
    N=len(x)
    mat=np.zeros((n,n,3),dtype=float)
    mat2=np.zeros((n,n),dtype=float)
    #another choice: b=sigma*sigma
    b=2*sigma*sigma
    for i in range(n):
        for j in range(n):
            k=0
            while k<N:
                #this r is actually squared
                r=((x[k]-i/(n-1))**2+(y[k]-j/(n-1))**2)
                w=np.exp(-r/b)
                mat[i,j,0]=mat[i,j,0]+w*vx[k]
                mat[i,j,1]=mat[i,j,1]+w*vy[k]
                mat[i,j,2]=mat[i,j,2]+w*vz[k]
                #normalization factor equals w if unit velocities are used
                mat2[i,j]=mat2[i,j]+w
                k=k+1
                    
            #mat[i,j]=(np.sqrt(mat[i,j][0]**2+mat[i,j][1]**2),arctan2(mat[i,j][1]/mat[i,j][0]))
    mat[:,:,0]=mat[:,:,0]/mat2
    mat[:,:,1]=mat[:,:,1]/mat2
    mat[:,:,2]=mat[:,:,2]/mat2
    return mat

def totaltri_gaussker_gaia(ra,dec,dist,ravel,decvel,radvel,vtot,n,sigma):
    """Creates a vectorial gauss kernel map of a grid of n rows and n columns
        given a standard deviation sigma, given the components of the vector.
        The velocities are treated as unitary vectors.
        Returns a 3-dimensional matrix whose third component is the three components of the vector.
        Computational cost: O(2n*n*N)
    """






