import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from particle import Particle
from field import DistanceField
from grid import Grid
import numpy as np
import pandas as pd
from numpy import linalg as LA
import time

#need to convert to np array for speeds

class particleManager:
    def __init__(self, width, height, numParticles, radius):
        self.RADIUS = radius
        self.COLLISIONRADIUS = 0
        self.DENSITYREST = 0.05
        self.SIGMA = 0
        self.BETA = 0.1

        #0-0.2
        self.YIELDRATIO = 0.1
        self.ALPHA = 0.3
        self.K = 0.7
        self.KSPRING = 0.03
        self.KNEAR = 5
        self.GRAVITY = [0, -9.8]
        self.FRICTION = 0.01
        self.COLLISIONSOFTNESS = 0.01

        #numParticles * width * height = actual number particles
        self.numParticles = numParticles
        self.width = width
        self.height = height
        self.distanceField = DistanceField(width, height, 0.1)  
        self.particles=self.generateParticles()
        self.grid = Grid(width, height, self.particles, self.RADIUS)
        self.neighbors = [[] for _ in range(len(self.particles))]
        self.updateNeighbors()
        self.counter = 0
        self.springs = [{} for _ in range(len(self.particles))]

    def generateParticles(self):
        particles = []
        id = 0
        for x in range(self.width*self.numParticles):
            for y in range(self.height * self.numParticles):
                particles.append(Particle(id, [x/self.numParticles, y/self.numParticles]))
                id+=1
        return particles

    def update(self, t):
        print(self.counter)
        self.animate(self.counter)
        self.applyExternalForces(t)
        self.applyViscosity(t)
        self.advanceParticles(t)
        self.updateNeighbors()
        self.springAdjustment(t)
        self.applySprings(t)
        self.doubleDensityRelaxation(t)
        self.updateVelocity(t)
        #self.resolveCollisionsWithWorld(t)
        self.resolveCollisions(t)
        #self.updateVelocity(t)
        self.counter += 1

    def applyExternalForces(self, t):
        for particle in self.particles:
            particle.vel = np.add(particle.vel, np.multiply(t, self.GRAVITY))
    
    def applyViscosity(self, t):
        for particle in self.particles:
            for nID in self.neighbors[particle.id]:
                n = self.particles[nID]
                r = np.subtract(n.pos, particle.pos)
                rMAG = LA.norm(r)
                q = rMAG/self.RADIUS
                if q<1:
                    rUNIT = np.divide(r, LA.norm(r))
                    u = np.dot(np.subtract(particle.vel, n.vel), rUNIT)
                    if u>0:
                        I = np.multiply(rUNIT, t*(1-q)*(self.SIGMA*u + self.BETA*u**2))
                        particle.vel = np.subtract(particle.vel, np.divide(I,2))
                        n.vel = np.add(n.vel, np.divide(I,2))
    
    def advanceParticles(self, t):
        for particle in self.particles:
            particle.posPrev = particle.pos
            particle.pos = np.add(particle.pos, np.multiply(particle.vel, t))
            self.grid.moveParticle(particle)
    
    def updateNeighbors(self):
        for particle in self.particles:
            self.neighbors[particle.id] = []
            possibleNeighbors = self.grid.possibleNeighbors(particle)
            for n in possibleNeighbors:
                #FIX THIS
                if LA.norm(np.subtract(self.particles[n].pos, particle.pos)) < self.RADIUS:
                    self.neighbors[particle.id].append(n)

    def springAdjustment(self, t):
        for particle in self.particles:
            for nID in self.neighbors[particle.id]:
                n = self.particles[nID]
                q = LA.norm(np.subtract(n.pos, particle.pos))/self.RADIUS
                if q <1:
                    if nID not in self.springs[particle.id]:
                        self.springs[particle.id][nID] = self.RADIUS
                    d  = self.YIELDRATIO * self.springs[particle.id][nID]
                    r = np.subtract(n.pos, particle.pos)
                    rMAG = LA.norm(r)
                    L = self.springs[particle.id][nID]
                    if rMAG > L+d:
                        self.springs[particle.id][nID] = L + t*self.ALPHA*(rMAG-L-d)
                    elif rMAG < L-d:
                        self.springs[particle.id][nID] = L - t*self.ALPHA*(rMAG-L-d) 
        for i in self.springs:
            toDrop = []
            for key, value in i.items():
                if value > self.RADIUS:
                    toDrop.append(key)
            for key in toDrop:
                i.pop(key)
                
    def applySprings(self, t):
        for particle in self.particles:
            for springID in self.springs[particle.id]:
                j = self.particles[springID]
                r = np.subtract(j.pos, particle.pos)
                rMAG = LA.norm(r)
                rUnit = np.divide(r, rMAG)
                L = self.springs[particle.id][springID]
                D = np.multiply(t**2*self.KSPRING * (1-L/self.RADIUS) * (L-rMAG), rUnit)
                particle.pos = np.subtract(particle.pos, np.divide(D, 2))
                j.pos = np.add(j.pos, np.divide(D, 2))
    
    def doubleDensityRelaxation(self, t):
        for particle in self.particles:
            p = 0
            pNear = 0
            for nID in self.neighbors[particle.id]:
                n = self.particles[nID]
                r = np.subtract(n.pos, particle.pos)
                rMAG = LA.norm(r)
                q = 1-rMAG/self.RADIUS
                if q<1:
                    p+=q**2
                    pNear += q**3
            P = self.K*(p-self.DENSITYREST)
            PNear = self.KNEAR * pNear
            dx = [0,0]
            print(f"id: {particle.id} and {self.neighbors[particle.id]}")
            for nID in self.neighbors[particle.id]:
                n = self.particles[nID]
                r = np.subtract(n.pos, particle.pos)
                rMAG = LA.norm(r)
                q = 1-rMAG/self.RADIUS
                rUNIT = np.divide(r, rMAG)
                if q<1:
                    D = np.multiply(rUNIT, 0.5*(P*q+PNear*(q**2))*t**2)
                    n.pos = np.add(n.pos, D)
                    dx = np.subtract(dx, D)
            particle.pos = np.add(particle.pos, dx)

    def checkBoundary(self, position):
        if position[1] > self.height or position[1] < 0 or position[0] < 0 or position[1] > self.width:
            return True
            #OUTSIDE
        else:
            return False

    def resolveCollisions(self, t):
        for particle in self.particles:
            distance = self.distanceField.getDistance(particle)
            if distance < -self.COLLISIONRADIUS:
                normal = self.distanceField.getNormal(particle) 
                # incrementor = 1.5
                # counter = 0
                # while self.checkBoundary(particle.pos):
                #     counter += 1
                #     reversed = np.multiply(particle.vel, -1 * incrementor*counter)
                #     particle.pos = np.add(particle.pos, np.multiply(reversed, t))
                #     print(reversed)
                initVel = np.multiply(particle.vel, -1)
                project = np.divide(np.dot(np.multiply(normal, abs(distance)), initVel),(LA.norm(initVel)**2))
                projection = np.multiply(initVel, project)
                particle.pos = np.add(particle.pos, np.multiply(projection, t)) 
                particle.vel = np.multiply([particle.vel[1], -1* particle.vel[0]], self.COLLISIONSOFTNESS*distance)
          
    def resolveCollisionsWithWorld(self, t):
        for particle in self.particles:
            distance = self.distanceField.getDistance(particle)
            if distance < -self.COLLISIONRADIUS:
                normal = self.distanceField.getNormal(particle)
                vNormal = np.multiply(normal, abs(np.dot(particle.vel, normal)))
                vTangent = np.subtract(particle.vel, vNormal)
                I = np.subtract(vNormal, np.multiply(self.FRICTION, vTangent))
                particle.vel = np.add(particle.vel, np.multiply(I,t))
                if particle.pos[0]>self.width:
                    particle.pos[0] = self.width
                elif particle.pos[0]< 0:
                    particle.pos[0] = 0
                elif particle.pos[1] > self.height:
                    particle.pos[1] = self.height
                elif particle.pos[1] < 0:
                    particle.pos[1] = 0
                particle.vel = np.multiply([particle.vel[1], -1* particle.vel[0]], 0.01*distance)

    # def resolveCollisionsWithWorld(self, t):
    #     for particle in self.particles:
    #         distance = self.distanceField.getDistance(particle)
    #         if distance < self.COLLISIONRADIUS:
    #             vMAG = LA.norm(particle.vel)
    #             vUNIT = np.divide(particle.vel, vMAG)
    #             normal = self.distanceField.getNormal(particle)
    #             tangent = self.perpendicularCCW(normal)
    #             normal = np.multiply(normal, t*self.COLLISIONSOFTNESS*(abs(distance)+self.RADIUS))
    #             tangent = np.subtract(particle.vel, normal)
    #             tangent = np.multiply(tangent, t*self.FRICTION)
    #             particle.vel = np.subtract(particle.vel, normal)
    #             particle.vel = np.subtract(particle.vel, tangent)
    #             particle.pos = np.add(particle.pos, np.multiply(particle.vel, t))
    #         if particle.id==10:
    #             print(f"pos: {particle.pos} vel: {particle.vel}")

    def updateVelocity(self, t):
        for particle in self.particles:
            particle.vel = np.divide(np.subtract(particle.pos, particle.posPrev), t)

    def perpendicularCCW(self, vector):
        return [vector[1], vector[0]]

    def animate(self, counter):
        if counter % 1 == 0:
            plotter = [[],[]]
            for particle in self.particles:
                plotter[0].append(particle.pos[0])
                plotter[1].append(particle.pos[1])
            plt.scatter(plotter[0], plotter[1], s= 100, c=plotter[1], cmap = "Wistia_r")
            ids = [x.id for x in self.particles]
            # for i in range(len(ids)):
            #     plt.annotate(ids[i], (plotter[0][i], plotter[1][i]))
            plt.xlim((-1, self.width))
            plt.ylim((-1, self.height))
            plt.xticks([], [])
            plt.yticks([], [])
            plt.tight_layout()
            plt.draw()
            plt.pause(0.001)
            plt.clf()

bruh = particleManager(20, 20, 1, 1)
for i in range(500): 
    bruh.update(1/60)
