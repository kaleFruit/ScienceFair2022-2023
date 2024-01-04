#Credit: Månsson, Daniel. “Interactive 2D Particle-based Fluid Simulationfor Mobile Devices.” (2013).
#Credit:  S. Clavet, P. Beaudoin, and P. Poulin / Particle-based Viscoelastic Fluid Simulation

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from fluidSim.particle import Particle
from field import DistanceField
from fluidSim.grid import Grid
import numpy as np
import pandas as pd
from numpy import linalg as LA
import time

#need to convert to np array for speeds

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
        self.collisionSoftness = 0.4
        self.numParticles = numParticles
        self.width = width
        self.height = height
        self.distanceField = DistanceField(width, height, 1)  
        self.particles=self.generateParticles()
        self.grid = Grid(width, height, self.particles, self.radius)
        self.neighbors = [[] for _ in range(len(self.particles))]
        self.updateNeighbors()
        self.counter = 0

    def generateParticles(self):
        particles = []
        counter = 0
        for x in range(self.width*self.numParticles):
            for y in range(self.height * self.numParticles):
                particles.append(Particle(counter, [x, y]))
                counter+=1
        return particles

    def update(self, t):
        print(self.counter)
        print(self.particles[89])
        self.animate(self.counter)
        self.applyExternalForces(t)
        self.applyViscosity(t)
        self.advanceParticles(t)
        self.updateNeighbors()
        self.doubleDensityRelaxation(t)
        self.resolveCollisions(t)
        self.updateVelocity(t)
        self.counter += 1

    def applyExternalForces(self, t):
        forces = self.inputForce()
        for particle in self.particles:
            particle.updateVelocity(self.GRAVITY, t)
            particle.updateVelocity(forces[particle.id], t)
    
    def applyViscosity(self, t):
        for particle in self.particles:
            for n in self.neighbors[particle.id]:
                n = self.particles[n]
                vPN = np.subtract(n.pos, particle.pos)
                velInward = np.dot(np.subtract(particle.vel, n.vel), vPN)
                if velInward > 0:
                    length = LA.norm(vPN)
                    velInward = velInward/length
                    vPN = np.multiply(vPN, 1/length)
                    q=length/self.radius
                    I = np.multiply(vPN, -1*0.5*t*(1-q)*(self.sigma*velInward+self.beta*velInward**2))
                    particle.updateVelocity(I, t)
    
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
                q = LA.norm(np.subtract(particle.pos, n.pos))
                q /= self.radius
                if q<1:
                    p+=(1-q)**2
                    pNear += (1-q)**3
            P = self.k*(p-self.densityRest)
            PNear = self.k * pNear
            dx = [0,0]
            for n in self.neighbors[particle.id]:
                n = self.particles[n]
                q = LA.norm(np.subtract(particle.pos, n.pos))
                q /= self.radius
                if q<1:
                    r = np.subtract(particle.pos,n.pos)
                    r = np.divide(r, LA.norm(r))
                    D = np.multiply(r, 0.5*(t**2)*(P*(1-q)+PNear*((1-q)**2)))
                    n.pos = np.add(n.pos, D)
                    dx = np.subtract(dx, D)
            particle.pos = np.add(particle.pos, dx)

    def resolveCollisions(self, t):
        for particle in self.particles:
            distance = self.distanceField.getDistance(particle)
            if distance < 0:
                v = LA.norm(particle.vel)
                if v==0:
                    v = 0.001
                normal = self.distanceField.getNormal(particle)
                tangent = self.perpendicularCCW(normal)
                tangent = np.multiply(tangent, t*self.friction*np.dot(v, tangent))
                particle.pos = np.subtract(particle.pos, tangent)
                normal = np.multiply(normal, self.collisionSoftness*(abs(distance)*self.radius))
                particle.pos = np.add(particle.pos, normal)
                #particle.updateVelocity(normal, t)
                #particle.updateVelocity(tangent, t)
                #particle.updateVelocity(np.multiply(particle.vel, -1.1), t)
                # particle.vel = np.add(particle.vel, tangent)

    
    def updateVelocity(self, t):
        for particle in self.particles:
            #CHECK THIS
            particle.vel = np.multiply(np.subtract(particle.pos, particle.posPrev), 1)
        #print(self.particles[100].vel)
        pass

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

    def animate(self, counter):
        if counter % 1 == 0:
            plotter = [[],[]]
            for particle in self.particles:
                plotter[0].append(particle.pos[0])
                plotter[1].append(particle.pos[1])
            #plt.scatter(plotter[0], plotter[1], s= 100, c=plotter[1], cmap = "Wistia_r")
            ids = [x.id for x in self.particles]
            plt.scatter(plotter[0], plotter[1], s= 100, c=ids, cmap = "Wistia_r")
            # plt.ylim = self.height
            # plt.xlim = self.width
            plt.xticks([], [])
            plt.yticks([], [])
            plt.tight_layout()
            plt.draw()
            plt.pause(0.0001)
            plt.clf()


bruh = particleManager(10, 10, 1, 5)
for i in range(1000):
    bruh.update(0.1)
