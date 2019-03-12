#Maximisation of the area
#par l'algorithme GA

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
import matplotlib.pyplot as plt
import csv

fig  = plt.figure()
canv = fig.add_subplot(1, 1, 1)
canv.set_xlim(0, 500)
canv.set_ylim(0, 500)


Nb_Cycles = 100
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
def modifCoord(pos, minMaxSets, flag):
    temp4Pos =  sol2rect(pos)    
    for i in range(len(temp4Pos)):
        if temp4Pos[i][0] < minMaxSets[0]:
            return flag 
        if minMaxSets[1] < temp4Pos[i][0]:
            return flag
        if temp4Pos[i][1] < minMaxSets[2]:
            return flag 
        if minMaxSets[3] < temp4Pos[i][1]:
            return flag 
    return -1
def sol2rect(pos):
    tempPos = []
    tempPos.append(pos[1])
    tempPos.append(rotate([pos[1][0] - pos[0][0], pos[1][1] - pos[0][1]], -pos[2], pos[0]))
    tempPos.append(rotate([pos[1][0] - pos[0][0], pos[1][1] - pos[0][1]], math.pi, pos[0]))
    tempPos.append(rotate([tempPos[1][0] - pos[0][0], tempPos[1 ][1] - pos[0][1]], math.pi, pos[0]))
    return tempPos

def makingCentGrav(polygon):
    coor0 = [0, 0]
    coor1 = [0, 0]
    angl = 0
    for line in polygon:
        coor0[0] += (line[0][0])/len(polygon)
        coor0[1] += (line[0][1])/len(polygon)
        coor1[0] += (line[1][0])/len(polygon)
        coor1[1] += (line[1][1])/len(polygon)
        angl += line[2]/len(polygon)
    return [coor0, coor1, angl]


def calcArea(pop):
    temp4Pos =  sol2rect(pop)
    return calcDist(temp4Pos[0], temp4Pos[1]) * calcDist(temp4Pos[0], temp4Pos[3])

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
def getBest(pop):
    bestArea  = 0
    for linePop in pop:
        if calcArea(linePop) > bestArea:
            bestCoord = linePop
            bestArea  = calcArea(linePop)
    return bestCoord

def getWorst(pop):
    worstArea  = 0
    index = 0
    worstCoord = pop[0]
    for i, linePop in enumerate(pop):
        if calcArea(linePop) < worstArea:
            worstCoord = linePop
            worstArea  = calcArea(linePop)
            index = i
    return [worstCoord, index]

def makingInd(parents, minMaxSets, polygon):
    flag = 1

    centGrav = makingCentGrav(parents)
    child = []
    child = [np.array([centGrav[0][0], centGrav[0][1]]), np.array([centGrav[1][0], centGrav[1][1]]), centGrav[2]]
    while flag == 1:
        for i in range (3):
            if i == 0 or i == 1:
                dVec = [parents[1][i][0] - parents[0][i][0], parents[1][i][1] - parents[0][i][1]]
                a = math.sqrt(pow(dVec[1], 2)/(pow(dVec[0], 2) + pow(dVec[1], 2)))
                b = math.sqrt(pow(dVec[0], 2)/(pow(dVec[0], 2) + pow(dVec[1], 2)))
                child[i][0] = (parents[0][i][0] + parents[1][i][0])/2 + randnp.randn()*(parents[1][i][0] - parents[0][i][0]) + calcDist(dVec, parents[2][0])*randnp.randn()*a
                child[i][1] = (parents[0][i][1] + parents[1][i][1])/2 + randnp.randn()*(parents[1][i][1] - parents[0][i][1]) + calcDist(dVec, parents[2][0])*randnp.randn()*b
            else:
                child[i] = (parents[0][i] + parents[1][i])/2 + randnp.randn()*(parents[1][i] - parents[0][i])
                child[i] = abs(child[i]) % (math.pi*2)
        tempArray = sol2rect(child)
        flag = inpolygon(tempArray, polygon)
    return child

bestResultArea = []
bestResultCmbi = []
for iteration in range(30):
    minMaxSets         = getBounds(polygon)
    centerGrav         = makingCenterGravCircums(polygon)
    pop                = initPop(Nb_Indiv, polygon)

    child = []
    best = getBest(pop)
    worst = getWorst(pop)
    child = []


    for i in range(Nb_Cycles):
        for j in range(Nb_Indiv):
            flag = 0
            num1 = num2 = num3 = 0
            while num1 == num2 == num3 or calcDist(pop[num1][0], pop[num2][0]) == 0 or calcDist(pop[num1][0], pop[num2][0]) == 0:
                num1 = random.randint(0, Nb_Indiv - 1)
                num2 = random.randint(0, Nb_Indiv - 1)
                num3 = random.randint(0, Nb_Indiv - 1)
                flag += 1
                if flag == 10:
                    for i in range(3):
                        pop[i] = muration(pop[i], minMaxSets, polygon)
            parents = [pop[num1], pop[num2], pop[num3]]
            child.append(makingInd(parents, minMaxSets, polygon))
        if calcArea(best) < calcArea(getBest(child)):
            best = getBest(child)
            pop[num1] = best
        if calcArea(worst[0]) < calcArea(getBest(child)):
            pop[worst[1]] = getBest(child)
            worst = getWorst(pop)
    print(iteration)

    
    bestResultArea.append(calcArea(best))
    bestResultCmbi.append(best)


print("Best Result Area")
for i in range(30):
    print(bestResultArea[i])
print(bestResultCmbi)
for i in range(30):
    print(bestResultCmbi[i])
	
	
	
print("The Best Result Area are : ")
for i in range(30):
    print(bestResultArea[i])
print(bestResultCmbi)
for i in range(30):
    print(bestResultCmbi[i])


for ite in range(30):
    temp4Pos =  sol2rect(bestResultCmbi[ite])
    for i in range(len(temp4Pos)):
        temp4Pos[i] = tuple(temp4Pos[i])
    im = Image.new('RGB', (600, 600), "black")
    draw = ImageDraw.Draw(im)
    draw.polygon((polygon), fill=200, outline=(255, 0, 255))
    draw.polygon((tuple(temp4Pos)), fill=(255, 255, 255), outline=(255, 0, 255))
    im.save('/Users\ZBOOK\Desktop\heuristics\image\GA\GAresults' + str(ite) + '.jpg', quality=95)

with open('/Users\ZBOOK\Desktop\heuristics\csvRESULTSOFGA.csv', 'w') as f:
    writer = csv.writer(f, lineterminator='\n') 
    writer.writerow(bestResultArea)     
