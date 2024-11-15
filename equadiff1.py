from Imports import *  # Importe tous les modules nécessaires (assumé que 'Imports' contient les modules nécessaires)

'''
On prend ici le modèle d'une sphère de rayon r soumise uniquement au poids et à des frottemnts linéaires à la vitesse pour commencer.
'''


# Fonction décrivant les équations du mouvement balistique
def balistic_equ(t, val, r, mu, g, M):
    x, y, z , x_p, y_p, z_p = val
    # Équations du mouvement pour la position et la vitesse
    dydt = [x_p, y_p, z_p, 6 * np.pi * mu * r * x_p / M, 6 * np.pi * mu * r * y_p / M, 6 * np.pi * mu * r * z_p / M - g]
    return np.array(dydt)

# Fonction pour résoudre les équations du mouvement
def find_solution(conditions_init, constantes, t_span):
    r = constantes['r']
    mu = constantes['mu']
    g = constantes['g']
    M = constantes['M']
    # Résoudre les équations différentielles
    solution = solve_ivp(
        fun=lambda t, val: balistic_equ(t, val, r, mu, g, M),
        t_span=t_span,
        y0=conditions_init,
        method='RK45',
        dense_output=True)
    return solution

# Définition des constantes
constantes = {'r': 0.1,
              'mu': 18.5e-6,
              'g': 9.81,
              'M': 50}

# Conditions initiales et paramètres d'angle et de vitesse
def condition_init(v,alpha,beta=0,origin=(0,0,0)):
    conditions_init = [origin[0],origin[1] ,origin[2], v * np.cos(beta) * np.cos(alpha), v * np.cos(alpha) * np.sin(beta), v * np.sin(alpha)]
    return conditions_init

conditions_init=condition_init(500,np.radians(45))


def tronquer_liste(liste, seuil):
    nouvelle_liste = []
    for valeur in liste:
        if valeur < seuil:
            break  
        nouvelle_liste.append(valeur)
    return nouvelle_liste


def resultats(conditions_init, constantes, seuil, temps):
    results = find_solution(conditions_init, constantes, (0, temps))
    t = np.linspace(0, temps, 1000)
    # Obtient les valeurs de position et de vitesse à partir des résultats
    sol_t = results.sol(t)
    x, y, z, x_p, y_p, z_p = sol_t[0], sol_t[1], sol_t[2], sol_t[3], sol_t[4], sol_t[5]
    
    # Filtrage des valeurs
    z = tronquer_liste(z, seuil)
    y = y[:len(z)]
    x = x[:len(z)]
    x_p = x_p[:len(z)]
    y_p = y_p[:len(z)]
    z_p = z_p[:len(z)]
    t = t[:len(z)]
    return t, x, y, z, x_p, y_p, z_p

t, x, y, z , x_p, y_p, z_p = resultats(conditions_init, constantes, 0, 200)

