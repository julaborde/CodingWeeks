from Imports import *
import equadiff2 as equa 
#from equadiff_basique import resultats
#from equadiff_cor_frottements import resultats
#from equadiff_coriolis import resultats

def create_earth_sphere(image_file):
    # On crée une figure de taille 800 par 800 avec une couleur de background noire et une couleur de figure blanche
    fig = mlab.figure(size=(800, 800), bgcolor=(0, 0, 0), fgcolor=(1, 1, 1))

    # On charge la texture
    img = tvtk.JPEGReader()
    img.file_name = image_file
    texture = tvtk.Texture(input_connection=img.output_port, interpolate=1)

    # # On crée une sphère source avec comme rayon le rayon de la Terre et comme résolution angulaire Nrad
    r = 6371e3
    Nrad = 1800
    sphere = tvtk.TexturedSphereSource(radius=r, theta_resolution=Nrad, phi_resolution=Nrad)

    # On assemble la pipeline et on assigne la texture
    sphere_mapper = tvtk.PolyDataMapper(input_connection=sphere.output_port)
    sphere_actor = tvtk.Actor(mapper=sphere_mapper, texture=texture)
    fig.scene.add_actor(sphere_actor)


    #On souhaite visualiser les latitudes et longitudes donc on définit une grille d'angles 
    phi, theta = np.mgrid[-np.pi/2:np.pi/2:180j, 0.0:2.0 * np.pi:360j]

    # on fait des boucles for pour tracer les latittudes et longitudes selon un intervalle de degrés.
    for lat in range(-90, 91, 30):
        lat_rad = np.radians(lat)
        x_line = r * np.cos(lat_rad) * np.cos(theta)
        y_line = r * np.cos(lat_rad) * np.sin(theta)
        z_line = r * np.sin(lat_rad) * np.ones_like(theta)

        #Ici, le terminal renvoyait une erreur de shape sur z, on a donc eu recours à .flatten().
        mlab.plot3d(x_line.flatten(), y_line.flatten(), z_line.flatten(), color=(1, 1, 1), tube_radius=5000)
        #print("Shapes - x_line:", x_line.shape, "y_line:", y_line.shape, "z_line:", z_line.shape)
        #mlab.plot3d(x_line, y_line, z_line, color=(1, 1, 1), tube_radius=10)

    for lon in range(0, 361, 30):
        lon_rad = np.radians(lon)
        x_line = r * np.cos(phi) * np.cos(lon_rad)
        y_line = r * np.cos(phi) * np.sin(lon_rad)
        z_line = r * np.sin(phi) * np.ones_like(lon_rad)

        mlab.plot3d(x_line.flatten(), y_line.flatten(), z_line.flatten(), color=(1, 1, 1), tube_radius=5000)
        #print("Shapes - x_line:", x_line.shape, "y_line:", y_line.shape, "z_line:", z_line.shape)
        #mlab.plot3d(x_line, y_line, z_line, color=(1, 1, 1), tube_radius=10)

    # On réajuste l'angle de vue initial
    mlab.view(azimuth=0, elevation=90, distance=50000)

    return fig

def gravitational_acceleration(position):
    # Fonction pour calculer l'accélération gravitationnelle en fonction de la position
    r = np.linalg.norm(position)
    g0 = 9.81  # Accélération gravitationnelle à la surface de la Terre en m/s^2
    g = g0 * (6371e3 / r)**2  # Variation de l'accélération gravitationnelle avec l'altitude
    return -g * position / r

def ballistic_trajectory(conditions_init, temps_total, intervalle_temps=1):
   
    r_terre,long_init, lat_init, vitesse_init, angle_inclinaison, angle_azimuth = conditions_init
    omega = 7.292115e-5
    # Conversion des angles en radians
    angle_inclinaison_rad = np.radians(angle_inclinaison)
    angle_azimuth_rad = np.radians(angle_azimuth)

    # Conversion des coordonnées initiales en cartésiennes
    x_init = r_terre * np.cos(np.radians(lat_init)) * np.cos(np.radians(long_init) + np.pi)
    y_init = r_terre * np.cos(np.radians(lat_init)) * np.sin(np.radians(long_init)+ np.pi)
    z_init = r_terre * np.sin(np.radians(lat_init))

    # Calcul des composantes de la vitesse initiale
    vitesse_init_x = vitesse_init * np.cos(angle_inclinaison_rad) * np.cos(angle_azimuth_rad)
    vitesse_init_y = vitesse_init * np.cos(angle_inclinaison_rad) * np.sin(angle_azimuth_rad)
    vitesse_init_z = vitesse_init * np.sin(angle_inclinaison_rad)

    # Initialisation des listes pour stocker les positions au fil du temps
    positions_x = [x_init]
    positions_y = [y_init]
    positions_z = [z_init]

    # Intégration numérique pour calculer la trajectoire
    temps = 0
    while temps <= temps_total:
        temps += intervalle_temps

        position = np.array([positions_x[-1], positions_y[-1], positions_z[-1]])
        acceleration = gravitational_acceleration(position)

        # Mise à jour des positions en utilisant les équations de la dynamique
        positions_x.append(positions_x[-1] + vitesse_init_x * intervalle_temps)
        positions_y.append(positions_y[-1] + vitesse_init_y * intervalle_temps)
        positions_z.append(positions_z[-1] + vitesse_init_z * intervalle_temps)

        # Mise à jour de la vitesse en fonction de l'accélération gravitationnelle
        vitesse_init_x += acceleration[0] * intervalle_temps
        vitesse_init_y += acceleration[1] * intervalle_temps
        vitesse_init_z += acceleration[2] * intervalle_temps

        # Mise à jour de la position en fonction de la rotation de la Terre
        delta_longitude = omega * intervalle_temps
        positions_x[-1], positions_y[-1] = rotate_coordinates(positions_x[-1], positions_y[-1], delta_longitude)

        # Vérification de la collision avec la Terre
        distance_to_earth_center = np.sqrt(positions_x[-1]**2 + positions_y[-1]**2 + positions_z[-1]**2)
        if distance_to_earth_center < r_terre:
            break

    return positions_x, positions_y, positions_z

def rotate_coordinates(x, y, delta_longitude):
     # Fonction pour effectuer une rotation des coordonnées autour de l'axe z
    x_new = x * np.cos(delta_longitude) - y * np.sin(delta_longitude)
    y_new = x * np.sin(delta_longitude) + y * np.cos(delta_longitude)
    return x_new, y_new

def haversine_distance(lat1, lon1, lat2, lon2):
    # Convertir les latitudes et longitudes de degrés à radians
    lat1_rad, lon1_rad = np.radians(lat1), np.radians(lon1)
    lat2_rad, lon2_rad = np.radians(lat2), np.radians(lon2)

    # Calcul des différences entre les latitudes et longitudes
    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad

    # Formule de la distance haversine
    a = np.sin(dlat / 2)**2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2)**2
    c = 2 * np.arcsin(np.sqrt(a))

    # Rayon de la Terre en mètres (approximatif)
    radius_earth = 6371000

    # Calcul de la distance
    distance = radius_earth * c

    return distance

constantes = {'D': 8.121e-3, 'g': 9.81, 'M': 150, 'omega': 7.292115e-5}

def parametre_tir(constantes,temps):
    coordonnees_villes = {
    'Paris': {'latitude': 48.8566, 'longitude': 2.3522, 'altitude': 35},
    'Marseille': {'latitude': 43.2965, 'longitude': 5.3698, 'altitude': 15},
    'Lyon': {'latitude': 45.75, 'longitude': 4.85, 'altitude': 173},
    'Toulouse': {'latitude': 43.6047, 'longitude': 1.4442, 'altitude': 150},
    'Nice': {'latitude': 43.7102, 'longitude': 7.2620, 'altitude': 10},
    'Nantes': {'latitude': 47.2184, 'longitude': -1.5536, 'altitude': 27},
    'Strasbourg': {'latitude': 48.5734, 'longitude': 7.7521, 'altitude': 143},
    'Montpellier': {'latitude': 43.6110, 'longitude': 3.8767, 'altitude': 27},
    'Bordeaux': {'latitude': 44.8374, 'longitude': -0.5792, 'altitude': 9},
    'Lille': {'latitude': 50.6292, 'longitude': 3.0573, 'altitude': 20},
    'Rennes': {'latitude': 48.1173, 'longitude': -1.6778, 'altitude': 30},
    'Reims': {'latitude': 49.2583, 'longitude': 4.0317, 'altitude': 80},
    'Le Havre': {'latitude': 49.4938, 'longitude': 0.1077, 'altitude': 10},
    'Clermont-Ferrand': {'latitude': 45.7772, 'longitude': 3.0822, 'altitude': 330},
    'Aix-en-Provence': {'latitude': 43.5297, 'longitude': 5.4474, 'altitude': 210},
    'Grenoble': {'latitude': 45.1885, 'longitude': 5.7245, 'altitude': 214},
    'Angers': {'latitude': 47.4784, 'longitude': -0.5632, 'altitude': 20},
    'Dijon': {'latitude': 47.3220, 'longitude': 5.0415, 'altitude': 245},
    'Nîmes': {'latitude': 43.8374, 'longitude': 4.3601, 'altitude': 58},
    'Saint-Denis': {'latitude': 48.9362, 'longitude': 2.3574, 'altitude': 33}
}
    rep=str(input("Villes dans le top20 français?(mettre True or False) : "))
    if rep=='True':
        VilleA=input("Ville de tir ?")
        while VilleA not in coordonnees_villes :
            print ("Mettre une majuscule au début")
            VilleA=input("Ville de tir ?")
        VilleB=input("Ville cible ?")
        while VilleB not in coordonnees_villes :
            print ("Mettre une majuscule au début")
            VilleB=input("Ville cible ?")
        lat_A=coordonnees_villes[VilleA]['latitude']
        long_A=coordonnees_villes[VilleA]['longitude']
        altitude_A=coordonnees_villes[VilleA]['altitude']
        lat_B=coordonnees_villes[VilleB]['latitude']
        long_B=coordonnees_villes[VilleB]['longitude']
        altitude_B=coordonnees_villes[VilleB]['altitude']
   
    else:
        # Demander à l'utilisateur les coordonnées de A
        latitude_A = float(input("Entrez la latitude de A : "))
        longitude_A = float(input("Entrez la longitude de A : "))
        altitude_A = float(input("Entrez l'altitude de A : "))

        # Demander à l'utilisateur les coordonnées de B
        latitude_B = float(input("Entrez la latitude de B : "))
        longitude_B = float(input("Entrez la longitude de B : "))
        altitude_B = float(input("Entrez l'altitude de B : "))

    #Définition des constantes
    rayon_terre = 6.371e6

    # Point d'arrivée A
    point_A = [rayon_terre+altitude_A, long_A, lat_A]

    # Point d'arrivée B
    point_B = [rayon_terre+altitude_B, long_B, lat_B]
    ###precision=float(input("Precision (en km) ?"))
    ###precision_m=precision*1000
    vectvitesseA0=np.array([800,0,0])

    # Fonction d'objectif à minimiser
    def objective_function(vectvitesseA, constantes, point_A, point_B, temps):
        long_A=point_A[1]
        lat_A=point_A[2]
        pA=np.array([point_A[0], np.pi/2 - np.radians(lat_A), np.radians(long_A)])
        conditions_init_A=np.concatenate([pA,vectvitesseA])
        # Calcul des résultats avec les conditions initiales fournies pour le point A
        long_arriv= equa.resultats(conditions_init_A, constantes,altitude_B,temps)[2][-1]
        lat_arriv=90-equa.resultats(conditions_init_A, constantes,altitude_B,temps)[3][-1]
        distance=haversine_distance(lat_B,long_B,lat_arriv,long_arriv)
        return distance
    
    vectvitesseA0=np.array([800,0,0])
    # Utiliser la fonction minimize pour trouver les conditions initiales optimales pour A
    objective_partial = partial(objective_function, constantes=constantes, point_A=point_A, point_B=point_B, temps=temps)
    result = minimize(objective_partial, vectvitesseA0, method='Nelder-Mead')
    

    return (point_A,point_B,result.x,result.fun)

print(parametre_tir(constantes,200))
    