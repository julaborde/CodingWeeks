from Imports import *  # Importe tous les modules nécessaires (assumé que 'Imports' contient les modules nécessaires)


def balistic_equ(t, val, D, g, M, omega):
    r, theta, phy , r_p, theta_p, phy_p = val
    # Équations du mouvement pour la position et la vitesse
    dydt = [r_p,
            phy_p,
            theta_p,
            -g + 2*omega*phy_p*np.sin(theta)**2 - D/M*r_p**2 +r*phy_p**2 - r*(theta_p*np.sin(phy))**2,
            1/r * (2*r*omega*phy_p*np.sin(theta)*np.cos(theta) - D/M*(r*phy_p)**2 + r*theta_p**2*np.sin(phy)*np.cos(phy) - 2*r_p*phy_p),
            1/(r*np.sin(phy)) * (-2*omega*(r_p*np.sin(theta)+np.cos(theta)*r*theta_p) - D/M*(r*theta_p*np.sin(phy))**2 - 2*r_p*theta_p*np.sin(phy) - 2*r*phy_p*theta_p*np.cos(phy))]
    return dydt

# Fonction pour résoudre les équations du mouvement
def find_solution(conditions_init, constantes, t_span):
    D = constantes['D']         # D = 1/2*rho*S*Cd  avec S la surface projetée et Cd le coeffidcient de trainée (pour un sphère, on prendra Cd = 0.47 et S = np.pi*r**2 (r = 1m))
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

constantes = {'D' : 0.21,
              'g': 9.81,
              'M': 1000,
              'omega':7.292115e-5}

long = 2.16218           # degrès
lat = 48.71041         # degrés
altitude = 61           # m
rayon_terre = 6.371e6   # m


conditions_init = [altitude + rayon_terre, np.radians(long), np.pi/2 - np.radians(lat), 2050, 1e-4 , 1e-4 ]

def tronquer_liste(liste, seuil):
    nouvelle_liste = []
    for valeur in liste:
        if valeur < seuil:
            break  
        nouvelle_liste.append(valeur)
    return nouvelle_liste

def resultats(conditions_init,constantes, seuil, temps):
    results = find_solution(conditions_init, constantes, (0, temps))
    t = np.linspace(0, temps, 10000)
    # Obtient les valeurs de position et de vitesse à partir des résultats
    r = results.sol(t)[0]
    phy = results.sol(t)[1]
    theta = results.sol(t)[2]
    r_p = results.sol(t)[3]
    phy_p = results.sol(t)[4]
    theta_p = results.sol(t)[5]
    v = np.sqrt(r_p**2 + (r*theta_p)**2 + (r*phy_p*np.sin(theta)**2))
    z = r - rayon_terre
    z = tronquer_liste(z, seuil)
    r = r[:len(z)]
    theta = theta[:len(z)]
    phy = phy[:len(z)]
    v = v[:len(z)]
    theta_p = theta_p[:len(z)]
    phy_p = phy_p[:len(z)]
    t = t[:len(z)]
    return t, r, np.degrees(np.pi/2 - theta), np.degrees(phy), v


def distance(lon1, lat1, lon2, lat2):
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])

    # Calcul des différences de latitude et de longitude
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    # Formule de Haversine
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    # Rayon moyen de la Terre en kilomètres (pour une sphère)
    radius_earth = 6371.0

    # Distance géodésique
    distance = radius_earth * c
    return distance


if __name__ == "__main__":
    t, r, theta, phy, v = resultats(conditions_init, constantes, 0,  1000)
    distance = distance(long, lat, phy[-1], theta[-1])
    print("La distance parcourue est de : {} de km".format(distance))
    print(phy[-1])
    print(theta[-1])
    plt.plot(t, r-rayon_terre)
    plt.show()

