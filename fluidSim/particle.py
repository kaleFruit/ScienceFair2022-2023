import numpy as np

class Particle:
    def __init__(self, id, pos):
        self.pos = pos
        self.posPrev = self.pos
        self.vel = [0,0]
        self.id = id 
    
    def __repr__(self):
        return f"ID: {self.id}; Position: {self.pos}; Prev: {self.posPrev}; Velocity: {self.vel}\n"

    def updateVelocity(self, force, t):
        self.vel = np.add(self.vel, np.multiply(force, 0.5*t**2))
