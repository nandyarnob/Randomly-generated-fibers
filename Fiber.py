import math
import random
import numpy as np
import matplotlib.pyplot as plt
import time
import os
import shutil
from scipy.optimize import curve_fit
from scipy.stats import poisson, norm

class line():
    def __init__(self,startPoint,endPoint):
        '''Initiate the line class with required properties of the line.'''
        self.startPoint = startPoint
        self.endPoint = endPoint
        self.m = (self.startPoint[1] - self.endPoint[1])/(self.startPoint[0] - self.endPoint[0])
        self.c = self.startPoint[1] - self.m*self.startPoint[0]
        self.maxLimit = [max(p,q) for p,q in zip(self.startPoint,self.endPoint)]    #This sets the higher limits of x and y co-ordinates
        self.minLimit = [min(p,q) for p,q in zip(self.startPoint,self.endPoint)]    #This sets the lower limits of x and y co-ordinates
        self.intersections = []
        self.nodeSerial = []
        
class circle():
    def __init__(self,center_c,radii):
        '''Initiate the circle class with required properties of the line.'''
        self.center_c = center_c
        self.radii = radii
        self.intersections = []
        
        #Set the limits of x and y co-ordinates around the circle
        self.xLowerLimit = self.center_c[0] - self.radii
        self.xUpperLimit = self.center_c[0] + self.radii
        self.yLowerLimit = self.center_c[1] - self.radii
        self.yUpperLimit = self.center_c[1] + self.radii
        
        
        
def lengthOfLine(startPoint,endPoint):
    '''This is a program to calculate the distance between 2 points. Input to be given as (x1,y1) and (x2,y2)'''
    length = math.sqrt(sum([(x-y)**2 for x, y in zip(startPoint, endPoint)]))
    return length

def circleLineIntersections(circles,lines):
    h,k = circles.center_c
    r = circles.radii
    a = 1 + lines.m**2
    b = 2*(lines.m*lines.c - lines.m*k - h)
    c = h**2 + k**2 + lines.c**2 - r**2 - 2*lines.c*k
    d1 = lengthOfLine(circles.center_c,lines.startPoint)
    d2 = lengthOfLine(circles.center_c,lines.endPoint)
    
    if(d1<=r and d2<=r):
        return
    
    d = b**2 - 4*a*c
    
    if(d==0):
        x1 = -b/(2*a)
        x2 = x1
    elif(d>0):
        x1 = (-b + np.sqrt(d))/(2*a)
        x2 = (-b - np.sqrt(d))/(2*a)
    else:
        return
    
    if(x1>lines.maxLimit[0] or x1<lines.minLimit[0]):
        x1 = 0
    if(x2>lines.maxLimit[0] or x2<lines.minLimit[0]):
        x2 = 0
        
    y1 = lines.m*x1 + lines.c
    y2 = lines.m*x2 + lines.c

    if(d1 < r):
        if(x1 != 0):
            lines.startPoint = (x1,y1)
        elif(x2 != 0):
            lines.startPoint = (x2,y2)

    elif(d2 < r):
        if(x1 != 0):
            lines.endPoint = (x1,y1)
        elif(x2 != 0):
            lines.endPoint = (x2,y2)
    
    
def createLine(x,y,theta,lineLength,boundary):
    l = lineLength
    
    x1 = x + (l/2)*np.cos(theta)
    x2 = x - (l/2)*np.cos(theta)
    
    y1 = y + (l/2)*np.sin(theta)
    y2 = y - (l/2)*np.sin(theta)
    
    [x1,x2,y1,y2] = correctBoundary(x1,x2,y1,y2,boundary=boundary)
    return [x1,y1],[x2,y2]


def correctBoundary(*arg,boundary):
    temp = np.array(arg)
    temp[temp>boundary] = boundary
    temp[temp<0] = 0
    
    return temp

def solve2Dequations(line1,line2):
        D = np.array([[line1.m, line2.m],[-1, -1]])     #D is denoted as the co-efficent Matrix, consisting coefficient of x and y (a1,a2,b1,b2)
        C = np.array([-line1.c, -line2.c])              #C is the Matrix of the constants. (c1,c2)
        DX = np.array([C,D[1]])                               #Replace the coefficients of x with constants
        DY = np.array([D[0],C])                               #Replace the coefficients of y with constants
        
        #This section is to eliminate the lines which are away from each other. This will improve the calculation time
        if(line2.startPoint[0] > line1.maxLimit[0] and line2.endPoint[0] > line1.maxLimit[0]):
            return [None]
        elif(line2.startPoint[0] < line1.minLimit[0] and line2.endPoint[0] < line1.minLimit[0]):
            return [None] 
        elif(line2.startPoint[1] > line1.maxLimit[1] and line2.endPoint[1] > line1.maxLimit[1]):
            return [None]
        elif(line2.startPoint[1] < line1.minLimit[1] and line2.endPoint[1] < line1.minLimit[1]):
            return [None]
        
        if(np.linalg.det(D) == 0):
            return [None]
        
        else:
            X = np.linalg.det(DX)/np.linalg.det(D)
            Y = np.linalg.det(DY)/np.linalg.det(D)
        
        return [X,Y]
        
        
convertFactor = 10**(-3)
lineLength = 7*convertFactor
noOfCells = 1
radius = 7.6*convertFactor
U = 0.10*radius
totalLines = 36000
boundary = radius*30
lineDensity = (totalLines*lineLength)/(boundary**2)



for D in range(1):
    

    print(f'{totalLines} numbers of fibers and {noOfCells} cell will be created.')

    start = time.time()

    lineObjects = []
    percentRadius = 0.05
    minElementLength = 0.0001*convertFactor

    circleObject = []

    xCenter = np.linspace(0,boundary,noOfCells+2)[1:-1]    #Create x co-ordinate as uniformly distributed in the boundary
    yCenter = np.ones(noOfCells)*boundary/2                #Same for all the circle objects
    centers = [(x,y) for x,y in zip(xCenter,yCenter)]

    for i in range(noOfCells):
        circleObject.append(circle(center_c=centers[i],radii=radius))


    count = 0

    while True:
        t = 0
        x = np.random.uniform(0,boundary)
        y = np.random.uniform(0,boundary)
        for tempCircle in circleObject:       
            d = lengthOfLine((x,y),tempCircle.center_c)
            if(d<=radius):
                t = 1
                break
        if(t == 1):
            continue
        theta = np.random.uniform(0,np.pi)
        p1,p2 = createLine(x,y,theta,lineLength,boundary=boundary)
        tempLine = line(p1,p2)
        for tempCircle in circleObject:
            circleLineIntersections(tempCircle,tempLine)

        newLine = line(tempLine.startPoint,tempLine.endPoint)
        lineObjects.append(newLine)

        count += 1
        if count == totalLines:
            break

    NodeX = []
    NodeY = []
    NodeNumber = 0

    for i in range(len(lineObjects)-1):
        for j in range(i+1,len(lineObjects)):
            temp = solve2Dequations(lineObjects[i],lineObjects[j])
            if len(temp) < 2 :
                continue
            else:
                if(lineObjects[i].minLimit[0] <= temp[0] <= lineObjects[i].maxLimit[0] and lineObjects[j].minLimit[0] <= temp[0] <= lineObjects[j].maxLimit[0]):
                    if(lineObjects[i].minLimit[1] <= temp[1] <= lineObjects[i].maxLimit[1] and lineObjects[j].minLimit[1] <= temp[1] <= lineObjects[j].maxLimit[1]):
                        NodeNumber += 1
                        lineObjects[i].nodeSerial.append(NodeNumber)
                        lineObjects[j].nodeSerial.append(NodeNumber)
                        NodeX.append(temp[0])
                        NodeY.append(temp[1])

        if((i+1)%100 == 0):
            print(f'{i+1} numbers of fibers checked')

    for lines in lineObjects:
        for circle1 in circleObject:
            d4 = lengthOfLine(lines.startPoint,circle1.center_c)
            if((d4 <= (1+percentRadius)*radius) or (boundary in lines.startPoint) or (0 in lines.startPoint)):
                NodeNumber += 1
                lines.nodeSerial.append(NodeNumber)
                NodeX.append(lines.startPoint[0])
                NodeY.append(lines.startPoint[1])

            d5 = lengthOfLine(lines.endPoint,circle1.center_c)
            if((d5 <= (1+percentRadius)*radius) or (boundary in lines.endPoint) or (0 in lines.endPoint)):
                NodeNumber += 1
                lines.nodeSerial.append(NodeNumber)
                NodeX.append(lines.endPoint[0])
                NodeY.append(lines.endPoint[1])


    NodeCord = list(zip(NodeX,NodeY))


    elementNodeNumbers = []
    elementLength = []
    for lines in lineObjects:
        temp = []
        for j in lines.nodeSerial:
            temp.append(NodeCord[j-1])
        sorted_nodeList = [temp1 for temp1,temp2 in sorted(zip(lines.nodeSerial,temp), key = lambda n:n[1][0])]
        for k in range(len(sorted_nodeList) - 1):
            elementNodeNumbers.append([sorted_nodeList[k],sorted_nodeList[k+1]])
            elementLength.append(lengthOfLine(NodeCord[sorted_nodeList[k]-1],NodeCord[sorted_nodeList[k+1]-1]))

    elementNumber = list(range(1,len(elementNodeNumbers)+1))


    count = 1
    delElement = []
    for t in elementNodeNumbers:
        i = t[0] - 1
        j = t[1] - 1

        if (lengthOfLine(NodeCord[i],NodeCord[j]) <= minElementLength):
            delElement.append(count)
        count += 1

    for i in sorted(delElement, reverse=True):
        del elementNodeNumbers[i-1]
        del elementNumber[i-1]
        del elementLength[i-1]

    nodeList = list(range(1,len(NodeCord)+1))

    XYCircle = []
    X0 = []
    Y0 = []
    XL = []
    YL = []
    Ux = []
    Uy = []
    Yc = []
    Xc = []

    tempRadius = (1+percentRadius)*radius
    for i in range(len(NodeCord)):
        for tempCircle in circleObject:
            if lengthOfLine(NodeCord[i],tempCircle.center_c) <= tempRadius:
                XYCircle.append(i+1)
                if (tempCircle.center_c[0] == NodeCord[i][0]):
                    ux = 0
                    uy = U                
                else:
                    m = abs(tempCircle.center_c[1] - NodeCord[i][1])/abs(tempCircle.center_c[0] - NodeCord[i][0])
                    theta = np.arctan(float(m))
                    ux = U*np.cos(theta)
                    uy = U*np.sin(theta)


                if (tempCircle.center_c[1] - NodeCord[i][1]) < 0:
                    uy = -uy
                if (tempCircle.center_c[0] - NodeCord[i][0]) < 0:
                    ux = -ux

                Ux.append(ux)
                Uy.append(uy)    

        if(NodeX[i] == 0):
            X0.append(i+1)
        elif(NodeX[i] == boundary):
            XL.append(i+1)
        if(NodeY[i] == 0):
            Y0.append(i+1)
        elif(NodeY[i] == boundary):
            YL.append(i+1)
        
        lowL = boundary/2 - (boundary/2)*percentRadius*0.1
        highL = boundary/2 + (boundary/2)*percentRadius*0.1
        if(NodeX[i] >= lowL and NodeX[i] <= highL):
            Yc.append(i+1)
        if(NodeY[i] >= lowL and NodeY[i] <= highL):
            Xc.append(i+1)        

    end = time.time()

    print(f'It took {end - start} seconds to complete the program, with {len(NodeCord)} nodes and {len(elementNodeNumbers)} elements.')

    fileName = f'L{lineLength}_D{lineDensity}_B{boundary}_C{noOfCells}_R{radius}_U{U}_{D}'

    basePath = "D:/My Stuffs/P/New/20_06_22/"
    path = os.path.join(basePath,fileName)
    os.mkdir(path)

    destination = path + "//"

    with open(destination + 'element_Final.inp','w') as f:
        tempList = list(zip(elementNumber,elementNodeNumbers))
        for i in tempList:
            f.write(f'{i[0]}, {i[1][0]}, {i[1][1]}\n')

    with open(destination + 'nodes_Final.inp','w') as f:
        tempList = list(zip(nodeList,NodeCord))
        for i in tempList:
            f.write(f'{i[0]}, {i[1][0]}, {i[1][1]}\n')

    with open(destination + 'nodes_Final.txt','w') as f:
        tempList = list(zip(nodeList,NodeCord))
        for i in tempList:
            f.write(f'{i[0]}, {i[1][0]}, {i[1][1]}\n')

    with open(destination + 'Ux_Final.inp','w') as f:
        tempList = list(zip(XYCircle,Ux))
        for i in tempList:
            f.write(f'{i[0]}, 1, 1, {i[1]}\n')

    with open(destination + 'Uy_Final.inp','w') as f:
        tempList = list(zip(XYCircle,Uy))
        for i in tempList:
            f.write(f'{i[0]}, 2, 2, {i[1]}\n')

    with open(destination + 'Circle_Nodes_Final.inp','w') as f:
        for i in XYCircle:
            f.write(f'{i}\n')


    with open(destination + 'X0_Final.inp','w') as f:
        for i in X0:
            f.write(f'{i}\n')

    with open(destination + 'XL_Final.inp','w') as f:
        for i in XL:
            f.write(f'{i}\n')

    with open(destination + 'Y0_Final.inp','w') as f:
        for i in Y0:
            f.write(f'{i}\n')

    with open(destination + 'YL_Final.inp','w') as f:
        for i in YL:
            f.write(f'{i}\n')
            
    with open(destination + 'Xc.inp','w') as f:
        for i in Xc:
            f.write(f'{i}\n')
            
    with open(destination + 'Yc.inp','w') as f:
        for i in Yc:
            f.write(f'{i}\n')

    with open(destination + 'Parameter_Details.txt','w') as f:
            f.write(f'Length of the fiber: {lineLength}\nFiber density: {lineDensity}\nBoundary: {boundary}\nNo of cell: {noOfCells}\nRadius: {radius}\nU: {U}')

    copy_files = ['Draft.inp','DraftEnd.inp','create DT File.bat', 'run abaqus.bat', '123.py','new123.py','getU.bat','Plots.bat','DT_plot.py']
    for p in copy_files:
            src_path = "D:/My Stuffs/P/New/Base Files/" + p
            dst_path = destination + p
            shutil.copy(src_path, dst_path)

    files = ['Draft.inp','Ux_Final.inp','Uy_Final.inp','DraftEnd.inp']

    with open(destination+'First.inp','w') as f:
        for filename in files:
            with open(destination+filename,'r') as fnew:
                temp = fnew.readlines()
                f.writelines(temp)

    end = time.time()
    
    print(np.average(elementLength))
    timeCount = end - start
    print(f'It took {round(timeCount,3)} seconds or {round(timeCount/60,3)} minutes to complete the program, with {len(NodeCord)} nodes and {len(elementNodeNumbers)} elements.')
