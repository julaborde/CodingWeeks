import equadiff1 as eq
from Imports import *

constantes = {'r' : 0.1,
              'mu' : 18.5e-6,
              'g' : 9.81,
              'M' : 50}

alpha = 45
beta = 45
v = 800
conditions_init = [0, 0, 0, v*np.cos(beta)*np.cos(alpha), v*np.cos(alpha)*np.sin(beta), v*np.sin(alpha)]
seuil = 0
temps=1000

res = eq.resultats(conditions_init, constantes,seuil,temps)

fig, (ax1,ax2,ax3) = plt.subplots(3)
ax1.plot(res[0], res[1], label = 'x en fonction du temps')
ax2.plot(res[0], res[2], label = 'y en fonction du temps')
ax3.plot(res[0], res[3], label = 'z en fonction du temps')

ax1.set_title("Composantes cartésiennes en fonction du temps")
#ax1.label("X")
#ax2.label("Y")
#ax3.label("Z")
ax1.grid()
ax2.grid()
ax3.grid()
plt.show()

plt.figure("Trajectoire en 3D")
axes = plt.axes(projection="3d")
print(axes, type(axes))


axes.plot(res[1], res[2], res[3])


axes.set_xlabel("X")
axes.set_ylabel("Y")
axes.set_zlabel("Z")
axes.set_title("Trajectoire dans l'espace")
# Afficher le graphique en 3D
plt.show()