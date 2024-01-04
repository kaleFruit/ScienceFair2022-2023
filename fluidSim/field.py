import numpy as np
import math 
import pandas as pd 

class DistanceField:
    def __init__(self, width, height, kernel):
        #nx=2, ny =2
        #[0,0]
        #[0,0]
        self.x = int(width/kernel)
        self.y = int(height/kernel)
        self.extension = 0
        self.kernel = kernel
        self.unitVectors = [[-1/(2**0.5), 1/(2**0.5)], [0, 1], [1/(2**0.5), 1/(2**0.5)],
                            [-1, 0], [0,0], [1, 0],
                            [-1/(2**0.5), -1/(2**0.5)], [0,-1], [1/(2**0.5), -1/(2**0.5)]]
        self.distanceGrid = self.initGenDistance()
        self.normalGrid = self.initGenNorm()

    
    def initGenDistance(self):
        grid = np.zeros((self.x, self.y))
        for i in range(int(self.x/2)+1):
            for j in range(self.y):
                grid[i][j] = i*self.kernel
                grid[self.x-1-i][j] = i*self.kernel
        for j in range(int(self.y/2)+1):
            for i in range(self.x):
                if j < grid[i][j]:
                    grid[i][j] = j*self.kernel
                    grid[i][self.y-1-j] = j*self.kernel
        return grid
    
    def initGenNorm(self):
        grid = np.empty((self.x, self.y, 2))
        for row in range(self.x):
            for column in range(self.y):
                grid[row][column] = self.calculateNormal(row, column)
        
        test = pd.DataFrame(grid[:,:,0])
        test.to_csv("bruh.csv")

        return grid
    
    def genNorm(self, x, y):
        grid = np.empty((x, y, 2))
        for row in range(x):
            for column in range(y):
                grid[row][column] = self.calculateNormal(row, column)
        return grid

    def getDistance(self, particle):
        coorX = (math.ceil(particle.pos[0]/self.kernel) if particle.pos[0]>0 else math.floor(particle.pos[0]/self.kernel))
        coorY = (math.ceil(particle.pos[1]/self.kernel) if particle.pos[1]>0 else math.floor(particle.pos[1]/self.kernel))
        coors = [coorX, coorY]
        self.expand(coors)
        gridCoors = self.convertedCoors(coors)
        try:
            return self.distanceGrid[gridCoors[0]][gridCoors[1]]    
        except:
            print(particle.pos, gridCoors, self.distanceGrid.shape)

    def expand(self, coors):
        ls = []
        #print(f"coor {coors} and {self.distanceGrid}")
        if coors[0] < 0:
            difference = self.extension+coors[0]
            if difference >= 0:
                ls.append(0)
            else:
                ls.append(abs(difference))
        elif coors[0] >= 0:
            difference = (coors[0] - (self.x+self.extension-1))
            if difference <= 0:
                ls.append(0)
            else:
                ls.append(difference)
        
        if coors[1] < 0:
            difference = self.extension+coors[1]
            if difference >= 0:
                ls.append(0)
            else:
                ls.append(abs(difference))
        elif coors[1] >= 0:
            difference = (coors[1] - (self.y+self.extension-1))
            if difference <= 0:
                ls.append(0)
            else:
                ls.append(difference)
        
        maxExtension = max(ls)

        if maxExtension > 0:
            for _ in range (maxExtension):
                self.extension+=1
                self.distanceGrid = np.pad(self.distanceGrid, 1, "constant", constant_values=-self.extension*self.kernel)
                self.normalGrid = np.pad(self.normalGrid, ((1,1), (1,1), (0,0)), "constant", constant_values=np.empty(2))
                for col in range(self.extension, self.extension+self.x):
                    self.normalGrid[0][col] = [0,-1]
                    self.normalGrid[len(self.normalGrid)-1][col] = [0,1]
                for y in range(self.extension, self.extension + self.y):
                    self.normalGrid[y][0] = [1, 0] 
                    self.normalGrid[y][len(self.normalGrid[0])-1] = [-1, 0]
            #3,4,2,1
                for x in range(1, self.extension+1):
                    for y in range(1, self.extension+1):
                        if x == self.extension or y == self.extension:
                            resultVector = [-x, y]
                            magnitude = (resultVector[0]**2 + resultVector[1]**2)**0.5
                            #magnitude = 1
                            self.normalGrid[-(self.extension-y+1)][-(self.extension-x+1)] = [resultVector[0]/magnitude, resultVector[1]/magnitude]
                
                for x in range(1, self.extension+1):
                    for y in range(1, self.extension+1):
                        if x == self.extension or y == self.extension:
                            resultVector = [-x, -y]
                            magnitude = (resultVector[0]**2 + resultVector[1]**2)**0.5
                            #magnitude = 1
                            self.normalGrid[(self.extension-y)][-(self.extension-x+1)] = [resultVector[0]/magnitude, resultVector[1]/magnitude]
                
                for x in range(1, self.extension+1):
                    for y in range(1, self.extension+1):
                        if x == self.extension or y == self.extension:
                            resultVector = [x, y]
                            magnitude = (resultVector[0]**2 + resultVector[1]**2)**0.5
                            #magnitude = 1
                            self.normalGrid[-(self.extension-y+1)][(self.extension-x)] = [resultVector[0]/magnitude, resultVector[1]/magnitude]
                
                for x in range(1, self.extension+1):
                    for y in range(1, self.extension+1):
                        if x == self.extension or y == self.extension:
                            resultVector = [x, -y]
                            magnitude = (resultVector[0]**2 + resultVector[1]**2)**0.5
                            #magnitude = 1
                            self.normalGrid[(self.extension-y)][(self.extension-x)] = [resultVector[0]/magnitude, resultVector[1]/magnitude]


    def convertedCoors(self, coors):
        #after kerneled 
        if coors[1]< 0:
            return [self.y-1+self.extension-coors[1],coors[0] +self.extension]
        return [self.x+self.extension-1-coors[1],coors[0] +self.extension] 

    def calculateNormal(self, row, column):
        resultVector = [0,0]
        counter = 0
        for x in range(-1, 2):
            for y in range(-1, 2):
                if (row+x) >= 0 and (row+x) < len(self.distanceGrid) and (column+y) > 0 and (column+y) < len(self.distanceGrid):
                    resultVector = np.add(resultVector, np.multiply(self.unitVectors[counter], (self.distanceGrid[row+x][column+y])))
                counter+=1
        magnitude = (resultVector[0]**2 + resultVector[1]**2)**0.5
        if magnitude != 0:
            return [resultVector[0]/magnitude, resultVector[1]/magnitude]
        else:
            return [0,0]

    def getNormal(self, particle):
        coorX = (math.ceil(particle.pos[0]/self.kernel) if particle.pos[0]>0 else math.floor(particle.pos[0]/self.kernel))
        coorY = (math.ceil(particle.pos[1]/self.kernel) if particle.pos[1]>0 else math.floor(particle.pos[1]/self.kernel))
        coors = [coorX, coorY]
        if abs(coors[0])>int((len(self.normalGrid)-1)/2) or abs(coors[1])>int((len(self.normalGrid[0])-1)/2):
            self.expand(coors)
        gridCoors = self.convertedCoors(coors)
        return self.normalGrid[gridCoors[0]][gridCoors[1]]


# grid = DistanceField(3,3,1)
# test = pd.DataFrame(grid.distanceGrid[:,:])
# test.to_csv("bruh1.csv")
# grid.expand([-1,-3])
# test = pd.DataFrame(grid.distanceGrid[:,:])
# test.to_csv("bruh.csv")
# p = Particle(0, [0,-2])
# print(grid.getNormal(p))
# print(grid.convertedCoors(p.pos))
# print(grid.normalGrid.shape)
# print(grid.distanceGrid.shape)