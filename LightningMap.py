import wget
import numpy as np
from scipy.misc import imread, imresize
import matplotlib.pyplot as plt
import scipy.signal as signal


def getMapLightning():
    #url='http://www.aemet.es/imagenes_d/eltiempo/observacion/rayos/201811162200_r79g.gif'
    url = 'http://2.bp.blogspot.com/-MViQnuVJCXg/TcZOE0t5yEI/AAAAAAAACYI/igevZeeZwOI/s1600/rayos.gif'
    image = wget.download(url)
    image = imread(image)

    Q = [ [0 for i in range(0,368)]  for i in range(0, 285)]
    for i in range(64, 349):
        Q[i-64] = image[i][112:480]

    maps = imresize(Q, [368, 368])

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


def Convolve(matrix1, matrix2):

    convolve = np.zeros((368, 368))

    for i in range(len(matrix1)):
        for j in range(len(matrix1[0])):
            convolve[i][j] = matrix1[i][j]*matrix2[i][j]

    return convolve

#############################
#####         MAIN      #####
#############################

Map = getMapLightning()
filters = [ [255 for i in range(0,368)] for j in range(0,368)]

P = Convolve(Map, filters)

plt.imshow(P)
plt.show()
