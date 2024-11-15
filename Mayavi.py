"""L'objectif de cette fonctionnalité a été de découvrir la bibliothèque mayavi d'affichage scientifique, les sphères texturées de tvtk
et de visualiser par l'animation la trajectoire d'un missile suivant la courbure de la Terre. 

WARNING : il est nécessaire d'installer les bibliothèques mayavi et tvtk à l'aide de la commande 
pip install [package].
pip install vtk==9.2
pip install PyQt5
pip install mayavi
"""

from Imports import *
import equadiff4 as equa2
#from equadiff_basique import resultats
#from equadiff_cor_frottements import resultats
#from equadiff_coriolis import resultats
choice=input("Modele simplifié (equ1) ou modele avancé (autre) : ")
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

if __name__ == "__main__":

    # On génère la planète Terre sous la forme d'une sphère avec un map d'une image sphériquement paramétrée (une projection géographique Plate Carrée, 
    # qui est basée sur une grille de latitude et de longitude égales et non sur une projection de surface égale) de la Terre.
    image_file = "blue_marble_spherical.jpg"
    figure = create_earth_sphere(image_file)

    # Conditions initiales pour le tir balistique
    altitude = 0
    r = 6371e3
    conditions_init = [r + altitude, 2.33333, 48.866669, 7000, 47, 297]  # [rayon, longitude, latitude, vitesse initiale, angle_inclinaison, angle_azimuth]
    
    # Paramètres pour le calcul de la trajectoire
    temps_total = 5000  # en secondes
    intervalle_temps = 1  # en secondes

    # Calcul de la trajectoire balistique
    if choice=="equ1":
        trajectoire_x, trajectoire_y, trajectoire_z = ballistic_trajectory(conditions_init, temps_total, intervalle_temps)
    else:
        trajectoire_x, trajectoire_y, trajectoire_z = equa2.resultats()
    # Récupération des coordonnées du point d'impact
    lat_fin = np.degrees(np.arcsin(trajectoire_z[-1]/r))
    long_fin =  - np.degrees(2 * np.arctan(trajectoire_y[-1] / (np.sqrt(trajectoire_x[-1]**2 + trajectoire_y[-1]**2) - trajectoire_x[-1])) -np.pi) - 180
    print (lat_fin, long_fin)
    
    # Affichage du point de lancement, de la trajectoire et du point d'impact
    lat_init = np.radians(conditions_init[2])
    long_init = np.radians(conditions_init[1]) + np.pi
    x_init = conditions_init[0] * np.cos(lat_init) * np.cos(long_init)
    y_init = conditions_init[0] * np.cos(lat_init) * np.sin(long_init)
    z_init = conditions_init[0] * np.sin(lat_init)

    # Positionnons Pyongyang sur la carte qui était le point visé avant déviation liée à la rotation de la Terre
    lat_piong = np.radians(39.019444)
    long_piong = np.radians(125.738052)+ np.pi
    x_piong = conditions_init[0] * np.cos(lat_piong) * np.cos(long_piong)
    y_piong = conditions_init[0] * np.cos(lat_piong) * np.sin(long_piong)
    z_piong = conditions_init[0] * np.sin(lat_piong)
        

    mlab.points3d(x_init, y_init, z_init, color=(0, 0, 1), scale_factor=50000)
    mlab.plot3d(trajectoire_x, trajectoire_y, trajectoire_z, color=(1, 0, 0), tube_radius=10000)
    mlab.points3d(trajectoire_x[-1], trajectoire_y[-1], trajectoire_z[-1], color=(0, 0, 1), scale_factor=50000)
    mlab.points3d(x_piong, y_piong, z_piong, color=(0, 0, 1), scale_factor=50000)

    

    # Animation
    @mlab.animate(delay=10) # décorateur animation mlab
    
    # Fonction qui réalise l'animation en calculant d'abord une trajectoire et en récupérant les points à chaque instants
    def animation():

        #Récupération des positions successives d'une trajectoire
        positions_x, positions_y, positions_z = ballistic_trajectory(conditions_init, temps_total, intervalle_temps)

        x = np.array(positions_x)
        y = np.array(positions_y)
        z = np.array(positions_z)

        # Initialisation de l'affichage sous la forme d'un plot "vide"
        points = mlab.points3d(0, 0, 0, color=(0, 1, 0), scale_factor=50000)

        # Mise à jour de la trajectoire
        for i in range(len(positions_x)):
            points.mlab_source.reset(x=x[:i], y=y[:i], z=z[:i])
            yield

    # Run the animation
    if choice =="equ1":
        animation()
    mlab.show()