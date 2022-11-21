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
    plot_option = True
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
    for j in range (1,J):
        Uj = B.dot(Uj)
        if (plot_option):
            if (j + 1) % 50 == 0:
                plt.plot(x, np.exp(r*T)*Uj, color='blue', linestyle='dashed', linewidth=1)
    
    if (plot_option):      
        valin,  = plt.plot(x, W_0,  color='orange', linestyle='solid', linewidth=2)
        valfin, = plt.plot(x, np.exp(r*T)*Uj,  color='red', linestyle='solid', linewidth=2)
        valex, = plt.plot(x, sol, color='green', linestyle='solid', linewidth=2)
        
        plt.xlabel("Actif sous-jacent")
        plt.ylabel("Valeur de l'option")

        plt.xlim((0.,1))
        
        plt.legend([valin, valfin,valex], ["Valeur à l'écheance", "Valeur (approx)","Valeur exacte"])
        #plt.legend([valin, valfin], ["Valeur à l'écheance", "Valeur (approx)"])
        #plt.legend([ valfin,valex], [ "Valeur (approx)","Valeur exacte"])

        plt.title(method)
        nom_fichier = "opt_floating_" + method + ".png"
        plt.savefig(nom_fichier)
    
    erreur = sol-np.exp(r*T)*Uj
    erreur = np.linalg.norm(erreur)
    
    return erreur


    