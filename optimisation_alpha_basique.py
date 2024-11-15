import equadiff_basique as equa
import numpy as np

# Paramètres pour l'optimisation
v = 500.0
temps = 200.0
seuil = 0
constantes = {'r': 0.1,
              'mu': 18.5e-6,
              'g': 9.81,
              'M': 50}

def optimisation_angle_portée(res, v, constantes, seuil, temps):

    def fonction_objectif(alpha, v, constantes, seuil, temps):
        # Cette fonction renvoie la portée pour la maximisation, alpha est en radians
        conditions_init = equa.condition_init(v, alpha)
        portee = equa.resultats(conditions_init, constantes, seuil, temps)[1][-1]
        # Supposons que la portée soit la distance horizontale maximale atteinte
        return portee

    fonction_obj = np.vectorize(fonction_objectif)

    # Créez une liste d'angles de 0 à pi/2 avec un pas de res
    alpha_liste = np.arange(0, np.radians(90), res)
    alpha_liste = np.append(alpha_liste, np.radians(90))  # Ajoute pi/2 à la liste

    # Calculez les portées correspondantes
    portee_liste = fonction_obj(alpha_liste, v, constantes, seuil, temps )

    # Trouvez l'indice de l'angle qui donne la portée maximale
    i_max = np.argmax(portee_liste)

    # Obtenez l'angle et la portée maximales
    alphamax = alpha_liste[i_max]
    portee_max = portee_liste[i_max]

    return np.degrees(alphamax), portee_max

# Exemple d'utilisation
res = np.radians(1)  # Résolution en radians
alphamax, portee_max = optimisation_angle_portée(res, v, constantes, seuil, temps)

print(f"L'angle optimal est {alphamax} degrés avec une portée maximale de {portee_max} mètres.")