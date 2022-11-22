from price_simu import *
import scipy.stats as st


if __name__ == '__main__':
    
    # Parametres principaux
    K = 1.0        # strike
    T = 1.0        # echeance
    r = 1.0        # taux de l'actif sans risque
    sigma = 2.0    # volatilite du sous-jacent
    M = 2    #Stike
    
    N = 10
    J = 2000
    
    method = 'EE'
    #compute_solution(sigma,r,N,J,M,T,method)
    study_error = True
    if (study_error):
        
        #Compute error for a fixed time step J=200
        J_cte = 2000
        #Choose N such that k/h**2 <= 0.2 with h = 1/N and k = T/J
        k = T/J_cte
        N_space_st = np.arange(10,20,1) #20
        error_list = []
        for N_ in N_space_st:
            error = compute_solution(sigma,r,N_,J_cte,M,T,method)
            error_list.append(error)
        error_list = np.array(error_list)
        #Plot the error in log scale
        plt.figure()
        plt.plot(np.log(1/N_space_st),np.log(error_list))
        plt.xlabel("log(h)")
        plt.ylabel("log(error)")
        plt.title("Error in function of h for J = 2000 for "+method)
        plt.savefig("error_N_J2000"+method+".png")
        #linear regression between log(error) and log(h)
        lr_h = st.linregress(np.log(1/N_space_st),np.log(error_list))
        print("Linear regression for h:\n ",lr_h)
        #print("Polynomial regressionfor h:\n",np.polyfit(np.log(1/N_space_st),np.log(error_list),2,full=True))
        
        
        
        #Compute error for a fixed space step N=100
        N_cte = 12
        error_list_2 = []
        J_time_st = np.arange(1000,2000,5)
        for J_ in J_time_st:
            error = compute_solution(sigma,r,N_cte,J_,M,T,method)
            error_list_2.append(error)
        error_list_2 = np.array(error_list_2)
        #Plot the error in log scale
        plt.figure()
        plt.plot(np.log(T/J_time_st),np.log(error_list_2))
        plt.xlabel("log(k)")
        plt.ylabel("log(error)")
        plt.title("Error in function of k for N=10 for "+method)
        plt.savefig("error_J_N10"+method+".png")
        
        #linear regression between log(error) and log(k)
        lr_k = st.linregress(np.log(T/J_time_st),np.log(error_list_2))
        print("Linear regression for k: \n",lr_k)
        #print("Polynomial regression for k:\n",np.polyfit(np.log(T/J_time_st),np.log(error_list_2),2,full=True))
  
