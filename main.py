import numpy as np
import matplotlib.pyplot as plt
import gaussker

latmax=43.77
latmin=36.00
lonmin=-9.30
lonmax=4.31
n=368
m=368
#lectura de los textos
latitud=[]
longitud=[]
poblacion=[]
with open('lat.txt', 'r') as f:
    data = f.readlines()
for line in data:
    latitud.append(float(line))
with open('lon.txt', 'r') as f:
    data = f.readlines()
for line in data:
    longitud.append(float(line))
with open('pob.txt', 'r') as f:
    data = f.readlines()
for line in data:
    poblacion.append(int(line))

#Creacion de la matriz de poblacion n*m
matrix=np.zeros((n,m))
for mun in range(len(poblacion)):
    j= int((longitud[mun]-lonmin)*n//(lonmax-lonmin))
    i= int((latitud[mun]-latmin)*m//(latmax-latmin))
    matrix[i,j]=matrix[i,j]+poblacion[mun]

for i in range(n):
  for j in range(m):
    matrix[i,j]=matrix[i,j]

sigma=0.01

plt.figure(1)

mesh=plt.pcolormesh(matrix)
plt.colorbar(mesh)
plt.show()
plt.savefig('fig1.png')

mat=gaussker.weightedker(matrix,368,sigma)
plt.figure(2)
mesh=plt.pcolormesh(mat)
plt.colorbar(mesh)
plt.show()
plt.savefig('fig2.png')

with open('spaindensity.txt','wb') as f:
    for line in matrix:
        np.savetxt(f, line, fmt='%.2f')

with open('convoluted.txt','wb') as f:
    for line in mat:
        np.savetxt(f, line, fmt='%.2f')