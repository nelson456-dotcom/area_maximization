#Maximisation of the area
#par l'algorithme PSO

#created by Mohamed Nizar Driouich.

from scipy import *
from math import *
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import sys
import pyclipper
import numpy as np
import random
import functools
import numpy.random as randnp
from PIL import Image, ImageDraw
import copy
from sympy.geometry import Point, Polygon
import csv

fig  = plt.figure()
canv = fig.add_subplot(1, 1, 1)
canv.set_xlim(0, 500)
canv.set_ylim(0, 500)

Nb_Cycles = 700  
Nb_Indiv  = 20
polygon = ((10,10),(10,300),(250,300),(350,130),(200,10)) 
w      = 0.9
ro_max = 1





def getBounds(polygon):
    minMaxSets = [polygon[0][0], polygon[0][0], polygon[1][0], polygon[0][1]]
    for linePoly in polygon:
            if linePoly[0] < minMaxSets[0]:
                minMaxSets[0] = linePoly[0]
            elif linePoly[0] > minMaxSets[1]:
                minMaxSets[1] = linePoly[0]
            elif linePoly[1] < minMaxSets[2]:
                minMaxSets[2] = linePoly[1]
            elif linePoly[1] > minMaxSets[3]:
                minMaxSets[3] = linePoly[1]
    return minMaxSets


def inpolygon(pos, polygon):     
    np = len(polygon)
    inside = False
    for i in range(len(pos)):
        inside = False
        for i1 in range(np): 
            i2 = (i1+1) % np
            if min(polygon[i1][0], polygon[i2][0]) < pos[i][0] < max(polygon[i1][0], polygon[i2][0]):
                if (polygon[i1][1] + (polygon[i2][1]-polygon[i1][1])/(polygon[i2][0]-polygon[i1][0])*(pos[i][0]-polygon[i1][0]) - pos[i][1]) > 0:
                    inside = not inside

        if inside == 0:
            return 1
    return 0


def makingCenterGravCircums(polygon):
    temp = [0 for i in range(len(polygon[0]))]
    for linePoly in polygon:
        for i in range(len(linePoly)):
            temp[i] += linePoly[i]
    return temp

    
def makingCentGrav2(polygon):
    centGrav = [0, 0]
    for line in polygon:
        centGrav[0] += line[0]/len(polygon)
        centGrav[1] += line[1]/len(polygon)
    return centGrav

def rotate(coord, angle, coord1):
    temp = [coord[0]*math.cos(angle) - coord[1]*math.sin(angle), coord[0]*math.sin(angle) + coord[1]*cos(angle) ]
    temp[0] += coord1[0]
    temp[1] += coord1[1]
    return temp

def initOne(polygon):
    minMaxSets = getBounds(polygon)
    coord = []
    flag = 1
    centGrav = makingCentGrav2(polygon)
    while flag == 1:
        pos = []
        angle = random.uniform(0, math.pi)
        pos.append([random.uniform(minMaxSets[0], minMaxSets[1]), random.uniform(minMaxSets[2], minMaxSets[3])])
        pos.append([random.uniform(minMaxSets[0], minMaxSets[1]), random.uniform(minMaxSets[2], minMaxSets[3])])
        pos.append(angle)
        coord =  sol2rect(pos)    
        flag = inpolygon(coord, polygon)
    return pos


def initPop(nb, polygon):
    return [initOne(polygon) for i in range(nb)]
def calcDist(pop1, pop2):
    return math.sqrt(pow((pop1[0] - pop2[0]), 2) + pow((pop1[1] - pop2[1]), 2))
def modifCoord(pos, minMaxSets, flag, polygon):
    temp4Pos =  sol2rect(pos)
    if inpolygon(temp4Pos, polygon) == 1:
        return flag
    return -1

def sol2rect(pos):
    tempPos = []
    tempPos.append(pos[1])
    tempPos.append(rotate([pos[1][0] - pos[0][0], pos[1][1] - pos[0][1]], -pos[2], pos[0]))
    tempPos.append(rotate([pos[1][0] - pos[0][0], pos[1][1] - pos[0][1]], math.pi, pos[0]))
    tempPos.append(rotate([tempPos[1][0] - pos[0][0], tempPos[1 ][1] - pos[0][1]], math.pi, pos[0]))
    return tempPos


def calcArea(pop):
    temp4Pos =  sol2rect(pop)
    return calcDist(temp4Pos[0], temp4Pos[1]) * calcDist(temp4Pos[0], temp4Pos[3])

def updatePosition(pos, v):
    for i in range(len(pos)):
        if i == 0 or i == 1:
            for j in range(len(v[i])):
                pos[i][j] += v[i][j]
        else:
            pos[i] += v[i]
    return pos

def updateVelocity(pos, v, w, ro_max, indBest,  groupBest):
    for i in range(len(v)):
        if i == 0 or i == 1:
            for j in range(len(v[i])):
                v[i][j] = w*v[i][j] + random.uniform(0, ro_max)*(indBest[i][j] - pos[i][j]) + random.uniform(0, ro_max)*(groupBest[i][j] - pos[i][j]) 
        else:
            v[i] = w*v[i] + random.uniform(0, ro_max)*(indBest[i] - pos[i]) + random.uniform(0, ro_max)*(groupBest[i] - pos[i])
            v[i] = abs(v[i]) % (2*math.pi)
            if v[i] > math.pi:
                v[i] -= math.pi
    return v

def makingInitianBest(pop, centerGrav):
    tempIndBest   = [pop[0][0], pop[0][1], pop[0][2]]
    bestArea      = calcArea(tempIndBest)
    for linePop in pop:
        currArea      = calcArea(linePop)
        if bestArea < currArea:
            bestArea    = currArea
            tempIndBest = copy.deepcopy(linePop)
    tempGroupBest = copy.deepcopy(tempIndBest)
    return tempIndBest, tempGroupBest

def muration(pos, minMaxSets, polygon):
    flag = 1
    while  flag == 1:
        randTemp = random.randint (0, 3)
        if randTemp == 0:
            pos[0] = [random.uniform(minMaxSets[0], minMaxSets[1]), random.uniform(minMaxSets[2], minMaxSets[3])]
        elif randTemp == 1:
            pos[1] = [random.uniform(minMaxSets[0], minMaxSets[1]), random.uniform(minMaxSets[2], minMaxSets[3])]
        else:
            pos[2] = random.uniform(0, math.pi)
        flag = inpolygon(sol2rect(pos), polygon)
    return pos

bestResultArea = []
bestResultCmbi = []
for ite in range(30):
    minMaxSets         = getBounds(polygon)
    centerGrav         = makingCenterGravCircums(polygon)
    pop                = initPop(Nb_Indiv, polygon)
    indBest, groupBest = makingInitianBest(pop, centerGrav)
    bestArea           = calcArea(groupBest)
    velocity           = [[0, 0], [0, 0], 0]


    child = []
    

    for i in range(Nb_Cycles):
        flag2 = 0
        for j in range(Nb_Indiv):
            tempPos  = copy.deepcopy(pop[j])
            flag1 = 1
            while flag1 != 0:
                velocity = updateVelocity(tempPos, velocity, w, ro_max, indBest, groupBest)
                tempPos = updatePosition(tempPos, velocity)
                flag1 = modifCoord(tempPos, minMaxSets, flag1, polygon)
                flag1 += 1  
                if flag1 == 10:               
                    tempPos = muration(tempPos, minMaxSets, polygon)
                    flag1 = 1
            tempArea = calcArea(tempPos)
            if calcArea(pop[j]) < tempArea:
                flag2 = 1
                pop[j] = copy.deepcopy(tempPos)
                del tempPos
                tempBestArea = tempArea
        if flag2 == 1: 
            if bestArea < tempBestArea:
                bestArea  = tempBestArea
                groupBest = copy.deepcopy(pop[j])
        tempBestArea = 0
    print(ite)

    bestResultArea.append(calcArea(groupBest))
    bestResultCmbi.append(groupBest)


print("Best Result Area")
for i in range(30):
    print(bestResultArea[i])
print(bestResultCmbi)
for i in range(30):
    print(bestResultCmbi[i])




for ite in range(30):
    temp4Pos =  sol2rect(bestResultCmbi[ite])
    for i in range(len(temp4Pos)):
        temp4Pos[i] = tuple(temp4Pos[i])
    im = Image.new('RGB', (600, 600), "white")
    draw = ImageDraw.Draw(im)
    draw.polygon((polygon), fill=200, outline=(255, 255, 255))
    draw.polygon((tuple(temp4Pos)), fill=(255, 255, 255), outline=(255, 255, 0))
    im.save('/Users\ZBOOK\Desktop\Heuristics\image\PSO\PSOresults' + str(ite) + '.jpg', quality=95)

f = open('/Users\ZBOOK\Desktop\heuristics\csvRESULTSOFPSO.csv', 'w')

writer = csv.writer(f, lineterminator='\n')
writer.writerow(bestResultArea)

f.close()
