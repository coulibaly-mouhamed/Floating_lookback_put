from price_simu import *


if __name__ == '__main__':
    
    # Parametres principaux
    K = 1.0        # strike
    T = 1.0        # echeance
    r = 1.0        # taux de l'actif sans risque
    sigma = 2.0    # volatilite du sous-jacent
    M = 2    #Stike
    
    N = 50
    J = 500 
    
    method = 'EI'
    compute_solution(sigma,r,N,J,M,T,method)
  
