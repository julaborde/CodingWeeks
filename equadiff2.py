from Imports import *  # Importe tous les modules nécessaires (assumé que 'Imports' contient les modules nécessaires)

'''
On garde ici le modèle d'une sphère comme dans le 1er fichier equadiff_basique auquel on rajoute l'impacte de la force de Coriolis
en coordonnées sphériques, on garde pour le moment des frottements prortionnels à la vitesse.
'''

def balistic_equ(t, val, D, g, M, omega):
    r, theta, phy , r_p, theta_p, phy_p = val
    # Équations du mouvement pour la position et la vitesse
    dydt = [r_p,
            phy_p,
            theta_p,
            -g + 2*omega*phy_p*np.sin(theta)**2 - D/M*r_p +r*phy_p**2 - r*(theta_p*np.sin(phy))**2,
            1/r * (2*r*omega*phy_p*np.sin(theta)*np.cos(theta) - D/M*(r*phy_p) + r*theta_p**2*np.sin(phy)*np.cos(phy) - 2*r_p*phy_p),
            1/(r*np.sin(phy)) * (-2*omega*(r_p*np.sin(theta)+np.cos(theta)*r*theta_p) - D/M*(r*theta_p*np.sin(phy)) - 2*r_p*theta_p*np.sin(phy) - 2*r*phy_p*theta_p*np.cos(phy))]
    return dydt

# Fonction pour résoudre les équations du mouvement
def find_solution(conditions_init, constantes, t_span):
    D = constantes['D']         # D = 6*np.pi*r*mu  comme dans le fichier equadiff_basique
    g = constantes['g']
    M = constantes['M']
    omega = constantes['omega']
    # Résoudre les équations différentielles
    solution = solve_ivp(
        fun=lambda t, val: balistic_equ(t, val, D, g, M, omega),
        t_span=t_span,
        y0=conditions_init,
        method='RK45',
        dense_output=True)
    return solution

constantes = {'D' : 3.49e-5,
              'g': 9.81,
              'M': 150,
              'omega':7.292115e-5}

lat = 2.16218           # degrès
long = 48.71041         # degrès
altitude = 61           # m
rayon_terre = 6.371e3   # m

conditions_init = [altitude + rayon_terre, np.radians(lat), np.radians(long), 800, 0, 0]


def tronquer_liste(liste, seuil):
    nouvelle_liste = []
    for valeur in liste:
        if valeur < seuil:
            break  
        nouvelle_liste.append(valeur)
    return nouvelle_liste


def resultats(conditions_init,constantes, seuil, temps):
    results = find_solution(conditions_init, constantes, (0, temps))
    t = np.linspace(0, temps, 1000)
    # Obtient les valeurs de position et de vitesse à partir des résultats
    r = results.sol(t)[0]
    phy = results.sol(t)[1]
    theta = results.sol(t)[2]
    r_p = results.sol(t)[3]
    phy_p = results.sol(t)[4]
    theta_p = results.sol(t)[5]
    v = np.sqrt(r_p**2 + (r*theta_p)**2 + (r*phy_p*np.sin(theta)**2))
    z = r - rayon_terre
    z = tronquer_liste(z,seuil)
    theta = theta[:len(z)]
    phy = phy[:len(z)]
    v = v[:len(z)]
    theta_p = theta_p[:len(z)]
    phy_p = phy_p[:len(z)]
    t = t[:len(z)]
    return t, z, np.degrees(theta), np.degrees(phy), v

t, z, theta, phy, v = resultats(conditions_init, constantes, altitude, 170)
