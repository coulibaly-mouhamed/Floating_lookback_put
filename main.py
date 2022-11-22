from turtle import color
from price_simu import *
import scipy.stats as st


if __name__ == '__main__':
    
    # Parametres principaux
    K = 1.0        # strike
    T = 2.0        # echeance
    r = 1.0        # taux de l'actif sans risque
    sigma = 2.0    # volatilite du sous-jacent
    M = 2    #Stike
    
    N = 20
    J = 2000
    
    method = 'EI'
    #compute_solution(sigma,r,N,J,M,T,method)
    study_error = False
    if (study_error):
        
        #Compute error for a fixed time step J=200
        J_cte = 100
        print("Erreur sur h.....\n")
        #Choose N such that k/h**2 <= 0.2 with h = 1/N and k = T/J
        k = T/J_cte
        N_space_st = np.arange(100,1000,70) #20
        error_list = []
        for N_ in N_space_st:
            print("N = ",N_)
            error = compute_solution(sigma,r,N_,J_cte,M,T,method)
            error_list.append(error)
        error_list = np.array(error_list)
        #Plot the error in log scale
        plt.figure()
        val_err, = plt.plot(np.log(1/N_space_st),np.log(error_list),color='blue', linestyle='solid', linewidth=2)
        #val_ref, = plt.plot(np.log(1/N_space_st),np.log(1/N_space_st**2),color='red')
        plt.xlabel("log(h)")
        plt.ylabel("log(error)")
        plt.title("Error in function of h for J = 2000 for "+method)
        #plt.legend([val_err,val_ref],["Error","2log(h)"])
        plt.savefig("error_N_J2000"+method+".png")
        #linear regression between log(error) and log(h)
        lr_h = st.linregress(np.log(1/N_space_st),np.log(error_list))
        print("Linear regression for h:\n ",lr_h)
        #print("Polynomial regressionfor h:\n",np.polyfit(np.log(1/N_space_st),np.log(error_list),2,full=True))
        
        
        
        #Compute error for a fixed space step N=12
        N_cte = 12
        print("Erreur sur k.....\n")
        error_list_2 = []
        J_time_st = np.arange(10,200,20)
        for J_ in J_time_st:
            print("J = ",J_)
            error = compute_solution(sigma,r,N_cte,J_,M,T,method)
            error_list_2.append(error)
        error_list_2 = np.array(error_list_2)
        #Plot the error in log scale
        plt.figure()
        val_err, =plt.plot(np.log(T/J_time_st),np.log(error_list_2),color='blue', linestyle='solid', linewidth=2)
        #val_ref,  = plt.plot(np.log(T/J_time_st),np.log(T/J_time_st),color='red')
        
        plt.xlabel("log(k)")
        plt.ylabel("log(error)")
        plt.title("Error in function of k for N=10 for "+method)
        #plt.legend([val_err,val_ref],["Error","log(k)"])
        plt.savefig("error_J_N10"+method+".png")
        
        #linear regression between log(error) and log(k)
        lr_k = st.linregress(np.log(T/J_time_st),np.log(error_list_2))
        print("Linear regression for k: \n",lr_k)
        #print("Polynomial regression for k:\n",np.polyfit(np.log(T/J_time_st),np.log(error_list_2),2,full=True))
  
    test_stability = True
    if (test_stability):
        N_list = np.arange(10,100,10)
        T=2
        J = 100
        norm_inf=[]
        for N_ in N_list:
            h = 1/N_
            x = np.linspace(0,1,N_)
            A = get_A(sigma,h,r,x)
            k= T/J
            if (method == 'EE'):
                B = (np.eye(N_) - k*A)
                #Calcul norme infinie de B
                norm_inf_b = np.linalg.norm(B,np.inf)
                norm_inf.append(norm_inf_b)
            elif (method == 'EI'):
                B = (np.linalg.inv(np.eye(N_) + k*A))
                #Calcul norme infinie de B
                norm_inf_b = np.linalg.norm(B,np.inf)
                norm_inf.append(norm_inf_b)
            
            elif (method == 'CN'):
                B = (np.linalg.inv(np.eye(N_) + k*A/2).dot(np.eye(N_) - k*A/2))
                #Calcul norme infinie de B
                norm_inf_b = np.linalg.norm(B,np.inf)
                norm_inf.append(norm_inf_b)
        plt.figure()
        plt.plot(1/N_list,norm_inf)
        plt.xlabel("h")
        plt.ylabel("Norme infinie de B")
        plt.title("Norme infinie de B en fonction de N pour "+method)
        plt.savefig("norm_inf_"+method+".png")