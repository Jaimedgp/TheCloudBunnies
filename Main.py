import numpy as np
import gaussker
import wget
from scipy.misc import imread, imresize, imsave
import matplotlib.pyplot as plt
import scipy.signal as signal
import gaussker

def arribaEspana(n,m):
    latmax=43.77
    latmin=36.00
    lonmin=-9.30
    lonmax=4.31
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
    with open('spaindensity.txt','wb') as f:
        for line in matrix:
            np.savetxt(f, line, fmt='%.2f')

def leeSpain(path, m):

    Q = [ [0 for i in range(0,m)]  for i in range(0, m)]

    with open(path, 'r') as fr:
        q = fr.readlines()
        for i in range(0,m):
            for j in range(0,m):
                try:
                    Q[i][j] = float(q[j+(m*i)])
                except ValueError:
                    print 'Error'
    return Q

def getMapLightning(m):
    #url='http://www.aemet.es/imagenes_d/eltiempo/observacion/rayos/201811162200_r79g.gif'
    url = 'http://2.bp.blogspot.com/-MViQnuVJCXg/TcZOE0t5yEI/AAAAAAAACYI/igevZeeZwOI/s1600/rayos.gif'
    image = wget.download(url)
    image = imread(image)

    Q = [ [0 for i in range(0,368)]  for i in range(0, 285)]
    for i in range(64, 349):
        Q[i-64] = image[i][112:480]

    maps = imresize(Q, [m, m])

    for i in range(len(maps)):
        for j in range(len(maps[0])):
            a,b,c = maps[i][j]
            maxC = 80

            if a == b and b == c:
                maps[i][j] = 0
            elif a > maxC and b > maxC and c > maxC:
                maps[i][j] = 0
            else:
                maps[i][j] = 1
    return maps[:,:,0]


def Convolve(matrix1, matrix2, m):

    convolve = np.zeros((m, m))

    for i in range(len(matrix1)):
        for j in range(len(matrix1[0])):
            convolve[i][j] = matrix1[i][j]*matrix2[i][j]

    return convolve

if __name__ == '__main__':
    coste=5112
    m = 90 
    Map = getMapLightning(m)

    arribaEspana(m,m)
    poblacion = leeSpain('/home/jaimedgp/TheCloudBunnies/spaindensity.txt', m)

    P = Convolve(Map, poblacion, m)

    RawCost = P*coste
    colorCost = np.arctan(P/100000)
    imsave('color.png', colorCost)
    imsave('cost.png', RawCost)
