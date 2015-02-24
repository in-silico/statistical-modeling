import numpy as np
import pylab as pb
import scipy.stats as stats
import time
import matplotlib.pyplot as plt

#pb.ion()

#########################
## Question 1
# first 4 columns are input: Wing-length, Wing-width, Tail-length, Arm-length
# last 2 columns are outputs in seconds (flight time for two trials)
data = np.genfromtxt('lab1_data.csv',delimiter=',')

X = data[:,0:4]
F = np.mean(data[:,4:6],axis=1)[:,None]

#########################
## Question 2

def B(x):
    # function returning the matrix of basis functions evaluated at x
    #input:	  x, np.array with d columns
    #output:  a matrix of geberal term B_{i,j} = b_j(x_i)
    b0 = np.ones((x.shape[0],1))
    b1 = (x[:,3])[:,None]
    b2 = b1*b1
    B = np.hstack((b0,b1,b2))
    return(B)

def LR(X,F,B):
    #input:	  X, np.array with d columns representing the DoE
    #		  F, np.array with 1 column representing the observations
    #		  B, a function returning the (p) basis functions evaluated at x
    # 		  tau2, noise variance
    #output:  beta, estimate of coefficients np.array of shape (p,1)
    #		  covBeta, cov matrix of beta, np.array of shape (p,p)

    phi = B(X)
    tmp = np.linalg.inv(np.dot(phi.T,phi))
    beta = np.dot(tmp,np.dot(phi.T,F))
    N,D = X.shape
    diff = F - np.dot(phi,beta)
    sigma = (1 / (N - D)) * np.dot(diff,diff)
    covBeta = sigma**2 * tmp

    return(beta,covBeta)

#########################
## Question 3

def predLR(x,B,beta,covBeta):
    #function returning predicted mean and variance
    #input:	  x, np.array with d columns representing m prediction points
    #		  B, a function returning the (p) basis functions evaluated at x
    #		  beta, estimate of the regression coefficients
    # 		  covBeta, covariance matrix of beta
    #output:  m, predicted mean at x, np.array of shape (m,1)
    #		  v, predicted variance matrix, np.array of shape (m,m)

    phi = B(x)
    m = np.dot(phi,beta)
    v = np.dot(phi,np.dot(covBeta,phi.T))

    return(m,v)

def plotModel(x,m,v):
    #input:	  x, np.array with d columns representing m prediction points
    #		  m, predicted mean at x, np.array of shape (m,1)
    #		  v, predicted variance matrix, np.array of shape (m,m)
    x = x.flatten()
    m = m.flatten()
    v = np.diag(v)
    upper=m+2*np.sqrt(v)
    lower=m-2*np.sqrt(v)
    plt.plot(x,m,color="#204a87",linewidth=2)
    plt.fill(np.hstack((x,x[::-1])),np.hstack((upper,lower[::-1])),color="#729fcf",alpha=0.3)
    plt.plot(x,upper,color="#204a87",linewidth=0.2)
    plt.plot(x,lower,color="#204a87",linewidth=0.2)

def R2(X,F,B,beta):
    return(1-sum((F-np.dot(B(X),beta))**2)/sum((F-np.mean(F))**2))

#########################
## Question 4

def pvalue(beta,covBeta,X):
    df = X.shape[0] - len(beta)
    cdf = stats.t.cdf(np.abs(beta)/np.sqrt(np.diag(covBeta)),df)
    return(2*(1 - cdf))

#########################
## Question 4

#########################
## Question 6



#########################
## Our code

(beta, covBeta) = LR(X,F[:,0],B)
phi = B(X)
phi_x = phi[:,1]
x = (np.linspace(phi_x.min() - 1,phi_x.max() + 1,200))[:,None]
xmat = np.hstack((x,x,x,x))
m,v = predLR(xmat,B,beta,covBeta)

#plt.plot(xmat[:,0],m)
#plt.plot(phi_x,F[:,0],'rx')
#plt.show()

plotModel(phi_x,m,v)
