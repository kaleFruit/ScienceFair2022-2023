class Grid:
    def __init__(self, width, height, particleList, radius):
        self.cellDict = {}
        width = width/radius
        height = height/radius
        self.radius = radius
        self.particles = [[] for _ in range(len(particleList))]
        for i in range(len(particleList)):
            cell = (int(particleList[i].pos[0]/self.radius), int(particleList[1].pos[1]/self.radius))
            self.particles[i] = cell
            try: 
                self.cellDict[cell].append(particleList[i].id)
            except:
                self.cellDict[cell] = []
                self.cellDict[cell].append(particleList[i].id)
        
    def possibleNeighbors(self, particle):
        currentCell = self.particles[particle.id]
        neighbors = []
        for i in range(-1, 2):
            for j in range(-1, 2):
                # if (i == 0 and j == 0):
                #     pass
                # else:
                try:
                    neighbors.append(self.cellDict[(currentCell[0]+i, currentCell[1]+j)])
                except:
                    pass
        return [j for sub in neighbors for j in sub if j!=particle.id]
    
    def moveParticle(self, particle):
        prevCell = self.particles[particle.id]
        newCell = (int(particle.pos[0]/self.radius), int(particle.pos[1]/self.radius))
        self.particles[particle.id] = newCell
        self.cellDict[prevCell].remove(particle.id)
        try:
            self.cellDict[newCell].append(particle.id)
        except:
            self.cellDict[newCell] = []
            self.cellDict[newCell].append(particle.id)
    
    def __repr__(self):
        return f"particles: {self.particles} and cells: {self.cellDict}"