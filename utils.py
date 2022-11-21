
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm 


#---------Define the function to compute the matrix A for a given space step h and N the number of points----------------

def get_A(sigma,h,r,x):
    #----Inputs: sigma is the volatility, h is the space step, x is the space grid,r is the interest rate
    #----Output: A is the matrix A = A1+A2+A3
    N = len(x)
    A3 = r *np.eye(N)
    
    #A2 = np.diag(x)
    #A2 +=np.diag(-x[1:],k=+1)
    #A2[N-1,N-1] -= (x[N-1]*(1+h)) 
    #A2 *= r/h
    
   
    A2 = np.zeros((N,N))
    A2 += np.diag(x[:N-1],k=+1)
    A2 -= np.diag(x[1:],k=-1)
    
    x_square = np.power(x,2)
    A1 = np.diag(2*x_square)
    A1 += np.diag(-x_square[1:],k=-1)
    A1 += np.diag(-x_square[:N-1],k=1)
    A1[N-1,N-1] -= x_square[N-1]*(1+h)
    A1 *= (sigma**2)/(2*h**2)
    
    return A1+A2+A3
    

def phi_tilde(x):
    return 1-x


def d_m_plus(sigma,r,x,T,t):
    d = np.log(x)+(r+(sigma**2)/2)*(T-t)
    d /= sigma*np.sqrt(T-t)
    return d

def d_m_minus(sigma,r,x,T,t):
    d = np.log(x)+(r-(sigma**2)/2)*(T-t)
    d /= sigma*np.sqrt(T-t)
    return d
       
def solu_exact_2(sigma,r,t,T,x):
    if (x==0):
        return 1 #If 1>2r/sigma**2
    
    else:
        if (t==T):
            a = ((1/x)**(2*r/(sigma**2)))+np.exp(-r*T)
            a *= (sigma**2/(2*r))
            a *= np.exp(-r*T)
            a *= x
            return a 
        
        else : 
            
<<<<<<< HEAD
            d_m_plus = (np.log(x) + (r + (sigma**2/2))*(T-t))/(sigma*np.sqrt(T-t))
            d_m_moin = (np.log(x) + (r - (sigma**2/2))*(T-t))/(sigma*np.sqrt(T-t))
            d_m_r = d_m_plus - 2*(r/sigma)*(T-t)
            first_term = np.exp(-r*(T-t))*norm.cdf(-d_m_moin)
            second_term = -x*norm.cdf(-d_m_plus)
            third_term_1 = np.exp(r*(T-t))*norm.cdf(d_m_plus)
=======
            d_m_plus = (np.log(x) + (r + (sigma**2/2))*t)/(sigma*np.sqrt(t))
            d_m_moin = (np.log(x) + (r - (sigma**2/2))*t)/(sigma*np.sqrt(t))
            d_m_r = d_m_plus - 2*(r/sigma)*t
            first_term = np.exp(-r*t)*norm.cdf(-d_m_moin)
            second_term = -x*norm.cdf(-d_m_plus)
            third_term_1 = np.exp(r*t)*norm.cdf(d_m_plus)
>>>>>>> 72714aae1c427f664f6fe1e02495f05741755b16
            third_term_2 = -((1/x)**(2*r/(sigma**2)))*norm.cdf(d_m_r)
            third_term = third_term_1 + third_term_2
            third_term *= x*np.exp(-r*(T-t))*((sigma**2)/(2*r))  
            
            return first_term + second_term + third_term
        
        
def solu_exact(sigma,r,t,T,x):
    
    d1 = (np.log(.0000001 + x) + (r + .5*sigma*sigma)*T)/(sigma*np.sqrt(T))

    a1 = norm.cdf(-d1+sigma*np.sqrt(T)) #*np.exp(-r*T) enlevé
    a2 = - x*norm.cdf(-1.*d1)*np.exp(r*T)
    a3 = x*sigma*sigma/(2*r) #*np.exp(-r*T) enlevé
    f1 = -np.power(.0000001 + x,-2*r/(sigma*sigma))*norm.cdf(d1-(2*r/sigma)*np.sqrt(T))
    f2 = np.exp(r*T) *norm.cdf(d1)

    sol_ex = a1+a2+a3*(f1+f2)
    return sol_ex

solu_exact_vect = np.vectorize(solu_exact)
