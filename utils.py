import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm 


#---------Define the function to compute the matrix A for a given space step h and N the number of points----------------

def get_A(sigma,h,r,x):
    #----Inputs: sigma is the volatility, h is the space step, x is the space grid,r is the interest rate
    #----Output: A is the matrix A = A1+A2+A3
    N = len(x)
    A3 = r *np.eye(N)
    
    A2 = np.diag(x)
    A2 +=np.diag(-x[1:],k=-1)
    A2 *= -r/h
    
    x_square = np.power(x,2)
    A1 = np.diag(-2*x_square)
    A1 += np.diag(x_square[1:],k=-1)
    A1 += np.diag(x_square[:N-1],k=1)
    A1[N-1,N-1] += x_square[N-1]/(1-h)
    A1 *= -(sigma**2)/(2*h**2)
    
    return A1+A2+A3
    



def phi_tilde(x):
    return 1-x


def solu_exact(sigma,r,t,T,x):
    if (x==0):
        return np.exp(-r*(T)) #If 1>2r/sigma**2
    
    else:
        if (t==0):
            a = ((1/x)**(2*r/(sigma**2)))+np.exp(-r*T)
            a *= (sigma**2/(2*r))
            a *= np.exp(-r*T)
            a *= x
            return a 
        
        else : 
            d_1_tilde = (np.log(x) + (r + (sigma**2/2)*t))/(sigma*np.sqrt(t))
            first_term = np.exp(r*t)*norm.cdf(-d_1_tilde+sigma*np.sqrt(t))
            second_term = x*norm.cdf(-d_1_tilde)
            second_term = -second_term
            third_term_1 = np.exp(r*t)*norm.cdf(d_1_tilde)
            third_term_2 = ((1/x)**(2*r/(sigma**2)))*norm.cdf(d_1_tilde-(2*r/sigma)*np.sqrt(t))
            third_term = third_term_1 + third_term_2
            third_term *= x*np.exp(r*t)*((sigma**2)/(2*r))  
            
            return first_term + second_term + third_term  
    
    
solu_exact_vect = np.vectorize(solu_exact)