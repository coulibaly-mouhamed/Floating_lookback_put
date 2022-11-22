from utils import *


def compute_solution(sigma,r,N,J,M,T,method='EI'):
    """
    Compute the solution of a Floating lookback Put 
    Inputs: 
        sigma = volatility of the model
        N : the number of points for space discretization 
        J: Number of points for time discretization
        M: The stike the option 
        T: The maturity of the option
        method : method of time discretization used (EE,EI or CN)
    Outputs : 
        sol: The solution at time 0
        x: The space grid
        t: The time grid
    """
    plot_option = False
    h = 1/(N)
    k = T/(J)
    x = np.linspace(0,1,N)
    t = np.linspace(0,T,J)
    A = get_A(sigma,h,r,x)
    
    #----Compute the solution at time 0
    sol = solu_exact_vect(sigma,r,0,T,x)
    
    
    
    if (method == 'EE'):
        B = (np.eye(N) - k*A)
    
    elif (method == 'EI'):
        B = (np.linalg.inv(np.eye(N) + k*A))
    
    elif (method == 'CN'):
        B = (np.linalg.inv(np.eye(N) + k*A/2).dot(np.eye(N) - k*A/2))
        
    #----Initial condition
    W_0 = phi_tilde(x)
    
    #----Loop over time
    Uj = W_0
    err = np.linalg.norm(Uj-solu_exact_vect(sigma,r,T,T,x))
    max_err = err
    for j in range (1,J):
        Uj = B.dot(Uj)
        err = np.linalg.norm(Uj-solu_exact_vect(sigma,r,T-t[j],T,x))
        if (err>max_err):
            max_err = err
        if (plot_option):
            if (j + 1) % 200 == 0:
                plt.plot(x, Uj, color='blue', linestyle='dashed', linewidth=1)
    
    if (plot_option):      
        valin,  = plt.plot(x, W_0,  color='orange', linestyle='solid', linewidth=2)
        valfin, = plt.plot(x, Uj,  color='red', linestyle='solid', linewidth=2)
        valex, = plt.plot(x, sol, color='green', linestyle='solid', linewidth=2)
        
        plt.xlabel("x")
        plt.ylabel("Prix du contrat")

        plt.xlim((0.,1))
        
        plt.legend([valin, valfin,valex], ["Valeur à l'écheance", "Valeur (approx)","Valeur exacte"])
        #plt.legend([valin, valfin], ["Valeur à l'écheance", "Valeur (approx)"])
        #plt.legend([ valfin,valex], [ "Valeur (approx)","Valeur exacte"])

        plt.title(method)
        nom_fichier = "opt_floating_" + method + ".png"
        plt.savefig(nom_fichier)
    
    
    
    return max_err


    