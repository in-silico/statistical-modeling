import numpy as np
import pylab as pb
pb.ion()

#########################
## Question 1

# the input points X and Y are always arrays with d columns

def mattern32(X,Y,sigma2=1.,theta=.2):
    d = np.sqrt(np.sum((X[:,None,:]-Y[None,:,:])**2/theta**2,2))
    k = sigma2*(1+np.sqrt(3)*d)*np.exp(-np.sqrt(3)*d)
    return(k)

def mattern52(X,Y,sigma2=1.,theta=.2):
    d = np.sqrt(np.sum((X[:,None,:]-Y[None,:,:])**2/theta**2,2))
    k = sigma2*(1 + np.sqrt(5)*d + (5.0/3.0)*d**2)*np.exp(-np.sqrt(5)*d)
    return(k)

def constant(X,Y,sigma2=1.):
    k = sigma2*np.ones((X.shape[0],Y.shape[0]))
    return(k)

def brownian(X,Y,sigma2=1.):
    k = sigma2*np.fmin(X,Y.T)
    return(k)

def dirac_delta(X,Y,sigma2=1.):
    k = sigma2*np.all(X[:,None,:]==Y[None,:,:],axis=2)
    return(k)

def gaussian(X,Y,sigma2=1.,theta=2.3):
    d2 = np.sum((X[:,None,:]-Y[None,:,:])**2/theta**2,2)
    k = sigma2*np.exp(-d2/2.)
    return(k + np.eye(X.shape[0])*1e-6)

## plot kernel
#x = np.linspace(-1,1,200)[:,None]
#y = kern(x,np.zeros((1,1)),1,.2)
#pb.plot(x,y,linewidth=2)

#########################
## Question 2
def sampleGP(x,mu,kern,n,**kwargs):
    # return n sample paths from a centred GP N(mu,kern) evaluated at x
    L = np.linalg.cholesky(kern(x,x))
    X = np.random.normal(0,1,(n,mu.shape[0]))

    return np.dot(X,L.T) + mu


#########################
## Question 3

#########################
## Question 4
def predGPR(x,X,F,kern,**kwargs):
    #function returning predicted mean and variance
    #input:	  x, np.array with d columns representing m prediction points
    #		  X, np.array with d columns representing the DoE
    #		  F, np.array with 1 column representing the observations
    # 		  kern, a kernel function
    #         **kwargs, arguments that can be passed to the kernel function
    #output:  m, predicted mean at x, np.array of shape (m,1)
    #		  v, predicted variance matrix, np.array of shape (m,m)

    # ... to be completed ...

    return(m,v)

#########################
## Question 5
def ftest(x):
    return(np.sin(3*np.pi*x)+2*x)

def plotModel(x,m,v,**kwargs):
    x = x.flatten()
    upper=m+2*np.sqrt(v)
    lower=m-2*np.sqrt(v)
    pb.plot(x,m,color="#204a87",linewidth=2,**kwargs)
    pb.fill(np.hstack((x,x[::-1])),np.hstack((upper,lower[::-1])),color="#729fcf",alpha=0.3)
    pb.plot(x,upper,color="#204a87",linewidth=0.2)
    pb.plot(x,lower,color="#204a87",linewidth=0.2)

########################
## Our code

x = np.linspace(-5,5,100)[:,None]
mu = np.zeros(100)

y = sampleGP(x,mu,gaussian,1)
pb.plot(x[:,0],y[0,:])
