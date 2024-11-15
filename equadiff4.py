import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


                # Force : le Poids
M0=10000 # kg : la masse initial du missile
k=0.002 #kg/s
a= 40
b= 0.029 #s-1
g0= 9.81 # m.s-2 : accéleration gravitationnelle standard


def Mbrulee(t):
    return k*a*np.exp(b*t)*t

def Masse(t):
    return (M0-Mbrulee(t))

def Masse_p(t):
    return (-k*a*np.exp(b*t)*(1+b*t))

def g(r):
    return g0*(6371e3/r)**2

def Poids(r,t):
    return np.array([-Masse(t)*g(r),
                     0,
                     0])



             # rho variable : modèle U.S. Standard atmosphère

T0= 288.15 # K : la température au niveau de la mer
Rho0=1.225 # kg.m-3 : densité de l'air au niveau de la mer
L=0.0065 #  K.m : le gradient de température (lapse rate)
M=0.0289644 # kg.mol-1 : la masse molaire de l'air
P0=101325 # Pa : la pression au niveau de la mer
R=287.05 # J.(kg.K)-1
Rt=6371000 # m : le rayon de la terre

#def Rho(r):
    #return Rho0*(1-((r-Rt)*L)/T0)**((g0*M)/(R*L))



                    # Hypothèses supplémentaires à prendre dans un second temps
# Température variable T(r)=t_0-(3.5/2)*np.log(r-Rt)
# Gradient de température variable L(r)=3.5/r




                # Modélisation d'un vent
# Vent constant
def vent_constant(C1,C2,C3):
    return np.array([C1,
                     C2,
                     C3])

# Vent tourbillonant
def vent_tourb(r,Ct):
    return np.array([0,
                     0,
                     Ct*r])

# Vent polaire
def vent_pol(theta,Cp):
    return np.array([0,
                     Cp*theta,
                     0])

# Vent composite
def vent_comp(r,theta,C1,C2,C3,Cr,Cp,Ct):
    return np.array([C1+Cr*r,
                     C2+Cp*theta,
                     C3+Ct*r])




                # Force de traînée et force de portance

C_d=0.12 # les valeurs admisses sont entre 0.05 à 0.2
C_z=0.1 # les valeurs admisses sont entre 0 et 0.2
S= np.pi

# la norme du vecteur vitesse V
def Norme_V(r,r_p,theta,theta_p,phi,phi_p):
    return np.sqrt((r_p)**2
                   +(r*theta_p)**2
                   +(r*np.sin(theta)*phi_p)**2)

# Force de Trainée sans vent
def Trainee(r,r_p,theta,theta_p,phi,phi_p):
    D=-0.5*C_d*S*Rho0*Norme_V(r,r_p,theta,theta_p,phi,phi_p)
    return np.array([D*r_p,
                     D*r*theta_p,
                     D*r*np.sin(theta)*phi_p])

# force de portance sans vent
def Portance(r,r_p,theta,theta_p,phi,phi_p):
    Z=0.5*C_z*S*Rho0*(Norme_V(r,r_p,theta,theta_p,phi,phi_p))**2
    return np.array([Z,
                     0,
                     0]) # on prend vect(n) colinéaire a er

# aux vues de comment se comporte le système de resolution, il est préférable de définir en amont C1,C2,C3,Cr,Cp,Ct
# à moins que l'on arrive à les faire saisir en ligne de commande

# norme vecteur vitesse V moins le vecteur vent perturbant Vp

C1=0
C2=0
C3=0
Cr=0
Cp=0
Ct=0

def Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p):
    return np.sqrt((r_p-(C1+Cr*r))**2
                   +(r*theta_p-(C2+Cp*theta))**2
                   +(r*np.sin(theta)*phi_p-(C3+Ct*r))**2)



# Force de trainée F_d avec vent
def Trainee_Vp(r,r_p,theta,theta_p,phi,phi_p):
    D=-0.5*C_d*S*Rho0*Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p)
    return np.array([r_p-(C1+Cr*r),
                    r*theta_p-(C2+Cp*theta),
                    r*np.sin(theta)*phi_p-(C3+Ct*r)])

# Force de portance avec un vent
def Portance_Vp(r,r_p,theta,theta_p,phi,phi_p):
    Z=0.5*C_z*S*Rho0*Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p)
    return np.array([1-(C1+Cr*r),
                    -(C2+Cp*theta),
                    -(C3+Ct*r)])




                    # Force de pousée

P0=100000 # Pa
ve=2010 # m.s-1  #  mp.sqrt((2*1,2*8,314*3000)/((1,2-1)*0.07412) (Mprog=74,12 # g.mol-1 ; gamma=1.2 ; Tcomb=3000 # K )
pe=30000000 # Pa  
Ae=2*np.pi*1 # m : la surface d'éjection des gaz (on prend un missile de diamètre de 1 mètre)

# variation de la pression selon l'altitude
#def Pression(r):
   # return P0*(1-(L*(r-Rt)/T0))**((g(r)*M)/(R*L))

def Pousee(r,r_p,theta,theta_p,phi,phi_p,t):
    k=(Mbrulee(t)*ve+(pe-P0)*Ae)
    return np.array([k*np.cos(theta),
                     k*np.sin(theta)*np.cos(phi),
                     k*np.sin(theta)*np.sin(phi)])




                    # Force de Coriolis

w=7.292*10**(-5)# rad/s : rotation de la terre

def wproj(theta):
    return np.array([w*np.cos(theta),
                     (-1)*w*np.sin(theta),
                     0])

def Coriolis(r,r_p,theta,theta_p,phi,phi_p,t):
    wr,wtheta,wphi=wproj(theta)
    return np.array([(-2*Masse(t))*(wtheta*r*np.sin(theta)*phi_p-wphi*r*theta_p),
                     (-2*Masse(t))*(wphi*r_p-wr*r*np.sin(theta)*phi_p),
                     (-2*Masse(t))*(wr*r*theta_p-wtheta*r_p)])


        #### toutes les forces

def systeme1(t, u):
    r, r_p, theta, theta_p, phi, phi_p = u
   
    eq1 = r_p
    eq2 = -(Masse_p(t)/Masse(t))*r_p + (r*(theta_p)**2 + r*np.sin(theta)*(phi_p)**2) - g(r) + (1/Masse(t))*(-0.5*C_d*S*Rho0*Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p)*(r_p-(C1+Cr*r)) + 0.5*C_z*S*Rho0*Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p)*(1-(C1+Cr*r)) + (Mbrulee(t)*ve + (pe-P0)*Ae)*np.cos(theta) - 2*Masse(t)*((-1)*w*np.sin(theta)*r*np.sin(theta)*phi_p - 0*r*theta_p))
    eq3 = theta_p
    eq4 = -(Masse_p(t)/Masse(t))*theta_p + (np.sin(theta)*np.cos(theta)*(phi_p)**2 - (1/r)*2*r_p*theta_p) + (1/(Masse(t)*r))*(-0.5*C_d*S*Rho0*Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p)*(r*theta_p-(C2+Cp*theta)) + 0.5*C_z*S*Rho0*Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p)*(-(C2+Cp*theta)) + (Mbrulee(t)*ve + (pe-P0)*Ae)*np.sin(theta)*np.cos(phi) - 2*Masse(t)*(0*r_p - w*np.cos(theta)*r*np.sin(theta)*phi_p))
    eq5 = phi_p
    eq6 = -(Masse_p(t)/Masse(t))*phi_p - (1/r)*(2*r_p*phi_p + 2*theta_p*phi_p*(np.cos(theta)/np.sin(theta))) + (1/(Masse(t)*r*np.sin(theta)))*(-0.5*C_z*S*Rho0*Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p)*(r*np.sin(theta)*phi_p - (C3+Ct*r)) + 0.5*C_z*S*Rho0*Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p)*(-(C3+Ct*r)) + (Mbrulee(t)*ve + (pe-P0)*Ae)*np.sin(theta)*np.sin(phi) - 2*Masse(t)*(w*np.cos(theta)*r*theta_p - (-1)*w*np.sin(theta)*r_p))

    return [eq1, eq2, eq3, eq4, eq5, eq6]
print("Frappe pour une vitesse initiale de 7000m/s, angle d'inclinaison 47 degre et angle d'azimuth 297 (equ1)")
long = float(input(" Longitude : "))
lat = float(input(" Latitude : "))
# Conditions initiales
r0, r_p0, theta0, theta_p0, phi0, phi_p0 = Rt, 2, np.radians(long), 2.0, np.pi/2 - np.radians(lat), 0.0
u0 = [r0, r_p0, theta0, theta_p0, phi0, phi_p0]

# Temps
t0, t_fin = 0, 100  # Vous pouvez ajuster ces valeurs selon vos besoins
t_span = (t0, t_fin)
nb_points = 1000
temps1 = np.linspace(t0, t_fin, nb_points)

# Résoudre le système d'équations différentielles avec solve_ivp
resultats = solve_ivp(systeme1, t_span, u0, method='RK45', dense_output=True)

# Extraire les résultats

r1 = resultats.sol(temps1)[0]
r_p1 = resultats.sol(temps1)[1]
theta1 = resultats.sol(temps1)[2]
theta_p1 = resultats.sol(temps1)[3]
phi1 = resultats.sol(temps1)[4]
phi_p1 = resultats.sol(temps1)[5]



                #### toutes les forces sauf la pousée

def systeme2(t, u):
    r, r_p, theta, theta_p, phi, phi_p = u
   
    eq1 = r_p
    eq2 =  (r*(theta_p)**2 + r*np.sin(theta)*(phi_p)**2) - g(r) + (1/Masse(100))*(-0.5*C_d*S*Rho0*Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p)*(r_p-(C1+Cr*r)) + 0.5*C_z*S*Rho0*Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p)*(1-(C1+Cr*r))  - 2*Masse(100)*((-1)*w*np.sin(theta)*r*np.sin(theta)*phi_p - 0*r*theta_p))
    eq3 = theta_p
    eq4 =  (np.sin(theta)*np.cos(theta)*(phi_p)**2 - (1/r)*2*r_p*theta_p) + (1/(Masse(100)*r))*(-0.5*C_d*S*Rho0*Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p)*(r*theta_p-(C2+Cp*theta)) + 0.5*C_z*S*Rho0*Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p)*(-(C2+Cp*theta))  - 2*Masse(100)*(0*r_p - w*np.cos(theta)*r*np.sin(theta)*phi_p))
    eq5 = phi_p
    eq6 =  - (1/r)*(2*r_p*phi_p + 2*theta_p*phi_p*(np.cos(theta)/np.sin(theta))) + (1/(Masse(100)*r*np.sin(theta)))*(-0.5*C_z*S*Rho0*Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p)*(r*np.sin(theta)*phi_p - (C3+Ct*r)) + 0.5*C_z*S*Rho0*Norme_V_moins_Vp(r,r_p,theta,theta_p,phi,phi_p)*(-(C3+Ct*r))  - 2*Masse(100)*(w*np.cos(theta)*r*theta_p - (-1)*w*np.sin(theta)*r_p))

    return [eq1, eq2, eq3, eq4, eq5, eq6]

# Conditions initiales
r2, r_p2, theta2, theta_p2, phi2, phi_p2 = r1[-1], r_p1[-1], theta1[-1], theta_p1[-1], phi1[-1], phi_p1[-1]
u2 = [r2, r_p2, theta2, theta_p2, phi2, phi_p2]

# Temps
t0, t_fin = 100, 13000  # ajuster ces valeurs selon besoins
t_span = (t0, t_fin)
nb_points = 10000
temps2 = np.linspace(t0, t_fin, nb_points)

# Résoudre le système d'équations différentielles avec solve_ivp
resultats = solve_ivp(systeme2, t_span, u2, method='RK45', dense_output=True)

# Extraire les résultats

r2 = resultats.sol(temps2)[0]
r_p2 = resultats.sol(temps2)[1]
theta2 = resultats.sol(temps2)[2]
theta_p2 = resultats.sol(temps2)[3]
phi2 = resultats.sol(temps2)[4]
phi_p2 = resultats.sol(temps2)[5]



def resultats():
    t=np.concatenate((temps1,temps2))
    r=np.concatenate((r1,r2))
    theta= np.concatenate((theta1,theta2))
    phi=np.concatenate((phi1,phi2))
    t_f = []
    r_f = []
    theta_f = []
    phi_f = []
    for i in range(len(r)):
        if r[i] - Rt > 0:
            t_f.append(t[i])
            r_f.append(r[i])
            theta_f.append(theta[i])
            phi_f.append(phi[i]+np.pi)

    x = r_f * np.sin(theta_f) * np.cos(phi_f)
    y = r_f * np.sin(theta_f) * np.sin(phi_f)
    z = r_f * np.cos(theta_f)
    return x,y,z





'''
theta = np.radians(theta1)
phi = np.radians(phi1)
x1 = r1 * np.sin(theta1) * np.cos(phi1)
y1 = r1 * np.sin(theta1) * np.sin(phi1)
z1 = r1 * np.cos(theta1)

x2 = r2 * np.sin(theta2) * np.cos(phi2)
y2 = r2 * np.sin(theta2) * np.sin(phi2)
z2 = r2 * np.cos(theta2)

seuil = Rt + 50
z_f1 = [i for i in z1 if i >= seuil]
z_f2 = [i for i in z2 if i >= seuil]


x1 = x1[len(z_f1):]
y1 = y1[len(z_f1):]
z1 = z1[len(z_f1):]


x2 = x2[:len(z_f2)]
y2 = y2[:len(z_f2)]
z2 = z2[:len(z_f2)]


print(np.sqrt(    (x2[-1]-x1[0])**2  + (y2[-1]-y1[0])**2 +(z2[-1]-z1[0])**2             )//1000)
print(max(max(z1),max(z2))//1000)

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Scatter plot of 3D coordinates
ax.scatter(x1, y1, z1, c='b', marker='o')

ax.scatter(x2, y2, z2, c='g', marker='o')
plt.show()
'''
'''

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, c='r', marker='o')

# Set labels for each axis
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')


plt.show()
'''