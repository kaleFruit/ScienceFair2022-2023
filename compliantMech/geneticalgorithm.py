from PRBMstatic import Mesh
import pygad
import numpy as np
from numpy import *
import random 
from sklearn.preprocessing import StandardScaler

class GeneticAlgorithm:
    def __init__(self):
        self.nx = 5
        self.ny = 5 
        self.nNum = self.nx*self.ny
        self.chrLength = int(0.5*self.nNum*(self.nNum-1))
        self.scale = 300 #mm
        
        self.force = 9.8*1000 #1kg mass *9.8M/s/s ==> [N]
        self.k = 573 #N.mm/rad
        self.angleLimit = 0.3
        self.N = 10
        
        self.population_size = 5
        self.mutation_probability = 0.7
        self.generations = 5
        self.chromosomes = []
        #CREATE CROSSOVER PROBABILIITIES TOO

    def main(self):
        population = self.generate_population()
        for gen in range(self.generations):
            print(f"GEN: {gen}")
            parent1, parent2 = self.select_chromosomes(population)
            child1, child2 = self.crossover(parent1, parent2)
            if random.uniform(0, 1) < self.mutation_probability:
                child1 = self.mutate(child1)
            if random.uniform(0, 1) < self.mutation_probability:
                child2 = self.mutate(child2)
            population = [child1, child2] + population[2:]
        best = self.get_best(population)
        candidate = self.generateModel(best)
        print("FINAL")
        candidate.draw()
        xOPT = candidate.staticAnalysis()
        newCoors = xOPT.reshape(np.array(candidate.nodePos).shape)
        candidate.draw(newCoors)
        print(best)

    def generateModel(self, chromosome):
        return Mesh(self.nx, self.ny, self.scale, chromosome, outputNodes=[20], 
                         constrainedNodes = [9], inputNodes=[0], 
                         freedomOutput=[False, True],  freedomInput = [False, True], 
                         k=self.k, force=self.force)
    
    def generate_population(self):
        population = []
        for _ in range(self.population_size):
            chromosome = [0 for _ in range(self.chrLength)]
            numConnections = random.randint(3, int(np.log2(self.chrLength)*3))
            for _ in range(numConnections):
                chromosome[random.randint(0, self.chrLength-1)] = 1
            candidate = self.generateModel(chromosome)
            candidate.looseLink()
            while sum(candidate.createChr()) <= 3:
                for _ in range(numConnections):
                    chromosome[random.randint(0, self.chrLength-1)] = 1
                candidate = self.generateModel(chromosome)
                candidate.looseLink()
            candidate.detectDisconnect(0)
            chromosome = candidate.createChr()
            self.chromosomes.append(chromosome)
            population.append(chromosome)
        return population

    def calculate_fitness(self, chromosome):
        p1 = 2
        p2 = 0.5
        posConstant = 50
        candidate = self.generateModel(chromosome)
        candidate.looseLink()
        candidate.redudantLink()
        candidate.detectDisconnect(0)
        volumeRatio = candidate.volumeRatio()
        if candidate.links > self.N:
            print(f"Fitness: {p1*volumeRatio + posConstant}")
            return p1*volumeRatio + posConstant
        else:
            xOPT = candidate.staticAnalysis()
            if isinstance(xOPT, list):
                print(f"xopt: {xOPT}")
                displacement, maxDeform = candidate.results(xOPT)
                penaltyZeta = maxDeform - self.angleLimit
                print(f"displacement: {displacement} maxDeform {maxDeform}")
                print(f"Fitness: {-displacement + p1*volumeRatio + p2*max(penaltyZeta, 0) + posConstant}")
                return -displacement + p1*volumeRatio + p2*max(penaltyZeta, 0) + posConstant
            else:
                return p1*volumeRatio + posConstant
        
    def select_chromosomes(self, population):
        fitness_values = []
        for chromosome in population:
            fitness_values.append(self.calculate_fitness(chromosome))
        fitness_values = [float(i) for i in fitness_values]
        newWeights = []
        print(f"prev: {fitness_values}")
        for i, v in enumerate(fitness_values):
            newWeights.append((v, i))
        for currIndex in range(len(newWeights)):
            smallest = currIndex
            for index2 in range(currIndex+1, len(newWeights)):
                if newWeights[index2][0] < newWeights[smallest][0]:
                    smallest = index2
            (newWeights[currIndex], newWeights[smallest]) = (newWeights[smallest], newWeights[currIndex])
        newWeights.reverse()
        for i, val in enumerate(newWeights):
            fitness_values[int(val[1])] = i*2
        parent1, parent2 = random.choices(population, weights=fitness_values, k=2)
        return parent1, parent2

    def crossover(self, parent1, parent2):
        crossover_point = random.randint(0, self.chrLength-1)
        child1 = parent1[0:crossover_point] + parent2[crossover_point:]
        child2 = parent2[0:crossover_point] + parent1[crossover_point:]
        return child1, child2
    
    def mutate(self, chromosome):
        mutation_point = random.randint(0, self.chrLength-1)
        if chromosome[mutation_point] == 0:
            chromosome[mutation_point] = 1
        else:
            chromosome[mutation_point] = 0
        return chromosome
    
    def get_best(self, population):
        fitness_values = []
        for chromosome in population:
            fitness_values.append(self.calculate_fitness(chromosome))
        min_index = min(fitness_values)
        min_index = fitness_values.index(min_index)
        return population[min_index]
    
new = GeneticAlgorithm()
new.main()