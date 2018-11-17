import wget
from math import sqrt
import numpy as np
from scipy.misc import imread, imresize#, imshow
import matplotlib.pyplot as plt


def getMapLightning():
    url='http://www.aemet.es/imagenes_d/eltiempo/observacion/rayos/201811162200_r79g.gif'
    image = wget.download(url)
    image = imread(image)

    Q = [ [0 for i in range(0,368)]  for i in range(0, 285)]
    for i in range(64, 349):
        Q[i-64] = image[i][112:480]

    maps = Q
    plt.imshow(maps)
    plt.show()

    for i in range(len(maps)):
        for j in range(len(maps[0])):
            a,b,c = maps[i][j]
            maxC = 100

            if a == b and b == c:
                maps[i][j] = 0
            elif a > maxC and b > maxC and c > maxC:
                maps[i][j] = 0
            else:
                maps[i][j] = 1
    return maps
