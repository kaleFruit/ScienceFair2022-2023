#Credit: Månsson, Daniel. “Interactive 2D Particle-based Fluid Simulationfor Mobile Devices.” (2013).

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from fluidSim.particle import Particle
from field import DistanceField
from fluidSim.grid import Grid
import numpy as np
import pandas as pd
from numpy import linalg as LA


class particleManager:
    def __init__(self, width, height, numParticles, radius):
        #width = 2, height = 2
        #[0,0]
        #[0,0]
        self.radius = radius
        self.collisionRadius = 1
        self.densityRest = 1
        self.sigma = 1
        self.beta = 0.01
        self.k = 0.1
        self.kNear = 1
        self.GRAVITY = [0, -9.8]
        self.friction = 0.2
        self.collisionSoftness = 0.8
        #numParticles * width * height = actual number particles
        self.numParticles = numParticles
        self.width = width
        self.height = height

        self.distanceField = DistanceField(width, height, 1)  

        self.particles=self.generateParticles()
        self.grid = Grid(width, height, self.particles, self.radius)
        self.neighbors = [[] for _ in range(len(self.particles))]
        self.updateNeighbors()

    def generateParticles(self):
        particles = []
        counter = 0
        for x in range(self.width*self.numParticles):
            for y in range(self.height * self.numParticles):
                particles.append(Particle(counter, [x, y]))
                counter+=1
        return particles

    def update(self, t):
        self.applyExternalForces(t)
        self.applyViscosity(t)
        self.advanceParticles(t)
        self.updateNeighbors()
        self.doubleDensityRelaxation(t)
        print("2")
        self.resolveCollisions(t)
        print("3")
        self.updateVelocity(t)

    def applyExternalForces(self, t):
        forces = self.inputForce()
        for particle in self.particles:
            particle.vel = np.add(particle.vel, self.GRAVITY)
            particle.vel = np.add(particle.vel, forces[particle.id])
    
    def applyViscosity(self, t):
        for particle in self.particles:
            for n in self.neighbors[particle.id]:
                n = self.particles[n]
                vPN = np.subtract(n.pos, particle.pos)
                velInward = np.dot(np.subtract(n.vel, particle.vel), vPN)
                if velInward > 0:
                    length = LA.norm(vPN)
                    velInward = velInward/length
                    vPN = np.multiply(vPN, 1/length)
                    q=length/self.radius
                    I = np.multiply(vPN, 0.5*t*(1-q)*(self.sigma*velInward+self.beta*velInward**2))
                    particle.vel = np.subtract(particle.vel, I)
    
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
                #CHECK LESS THAN AND EQUALS or just LESS THAN
                if LA.norm(np.subtract(particle.pos, self.particles[n].pos)) <= self.radius:
                    self.neighbors[particle.id].append(n)

    def doubleDensityRelaxation(self, t):
        for particle in self.particles:
            p = 0
            pNear = 0
            for n in self.neighbors[particle.id]:
                n = self.particles[n]
                tempN = LA.norm(np.subtract(particle.pos, n.pos))
                q = 1 - tempN/self.radius
                p+= q**2
                pNear += q**3
            P = self.k * (p - self.densityRest)
            PNear = self.kNear * pNear
            dx = [0,0]
            for n in self.neighbors[particle.id]:
                n = self.particles[n]
                tempN = LA.norm(np.subtract(particle.pos, n.pos))
                q = 1-tempN/self.radius
                if tempN == 0:
                    tempN = 0.001
                vPN = np.multiply(np.subtract(particle.pos,n.pos), 1/tempN)
                D = np.multiply(vPN, 0.5*(t**2)*(P*q+PNear*(q**2)))
                n.pos = np.add(n.pos, D)
                dx = np.subtract(dx, D)
            particle.pos = np.add(particle.pos, dx)

    # def resolveCollisions(self, t):
    #     for particle in self.particles:
    #         distance = self.distanceField.getDistance(particle)
    #         print(distance)
    #         if distance < -self.collisionRadius:
    #             test= particle.pos
    #             for n in self.neighbors[particle.id]:
    #                 n = self.particles[n]
    #                 tempN = LA.norm(np.subtract(particle.pos, n.pos))
    #                 #print(tempN)
    #                 #FIX DIVISION BY ZERO
    #                 if tempN==0:
    #                     tempN = 0.001
    #                 vPN = np.multiply(np.subtract(particle.pos, n.pos), 1/tempN)
    #                 normal = np.multiply(self.distanceField.getNormal(particle), 1)
    #                 tangent = np.multiply(self.perpendicularCCW(normal), 1)
    #                 tangent = np.multiply(tangent, t*self.friction*np.dot(vPN, tangent))
    #                 particle.pos = np.subtract(particle.pos, tangent)
    #                 tempVector = np.multiply(normal, self.collisionSoftness*(abs(distance)+self.radius))
    #                 #particle.posPrev = particle.pos
    #                 particle.pos = np.add(particle.pos, tempVector)
    #             print(f"previous: {test} tempV: {particle.pos}")
    #     test = pd.DataFrame(self.distanceField.normalGrid[:,:, 1])
    #     test.to_csv("bruh.csv")

    def resolveCollisions(self, t):
        for particle in self.particles:
            distance = self.distanceField.getDistance(particle)
            if distance < 0:
                test = particle.pos
                v = LA.norm(particle.vel)
                if v==0:
                    v = 0.001
                normal = self.distanceField.getNormal(particle)
                tangent = self.perpendicularCCW(normal)
                tangent = np.multiply(tangent, t*self.friction*np.dot(v, tangent))
                particle.pos = np.subtract(particle.pos, tangent)
                normal = np.multiply(normal, self.collisionSoftness*(abs(distance)+self.radius))
                particle.pos = np.add(particle.pos, normal)
                particle.vel = np.subtract(np.multiply(particle.vel, -1), normal)
                particle.vel = np.add(particle.vel, tangent)
                print(f"previous: {test} particle: {particle}")
            print(self.distanceField.distanceGrid.shape)
        test = pd.DataFrame(self.distanceField.normalGrid[:,:, 1])
        test.to_csv("bruh.csv")

    # def resolveCollisions(self, t):
    #     for particle in self.particles:
    #         distance = self.distanceField.getDistance(particle)
    #         if distance < 0:
    #             bounce = [particle.vel[0], particle.vel[1]*-1]
    #             particle.vel = np.multiply(bounce, 0.3)
    #             tangent = self.perpendicularCCW(normal)
    #             particle.pos = particle.pos
    #             print(f"previous: {test} particle: {particle}")
    #         print(self.distanceField.distanceGrid.shape)
    #     test = pd.DataFrame(self.distanceField.normalGrid[:,:, 1])
    #     test.to_csv("bruh.csv")

    
    def updateVelocity(self, t):
        for particle in self.particles:
            particle.vel = np.multiply(np.subtract(particle.pos, particle.posPrev), 1/t)

    def inputForce(self):
        #FIX
        forces = [[0,0] for x in range(len(self.particles))]
        # if input("Force? ") == "y":
        #     for particle in self.grid.matrix[0][self.height]:
        #         forces[particle] = [5, 5]
        #     return forces
        # else:
        #     return forces
        return forces

    def perpendicularCCW(self, vector):
        return [vector[1], vector[0]]

    def animate(self, time):
        bruh.update(time)
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
        plt.pause(0.0001)
        plt.clf()


bruh = particleManager(5, 5, 1, 3)
for i in range(500): 
    bruh.animate(1/60)

