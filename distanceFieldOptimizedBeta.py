class MinHeap:
    def __init__(self, maxsize):
        self.maxsize = maxsize
        self.size = 0
        self.Heap = [0]*(self.maxsize + 1)
        self.FRONT = 1
  
    def parent(self, pos):
        return pos//2
  
    def leftChild(self, pos):
        return 2 * pos
  
    def rightChild(self, pos):
        return (2 * pos) + 1
  
    def isLeaf(self, pos):
        return pos*2 > self.size
  
    def swap(self, fpos, spos):
        self.Heap[fpos], self.Heap[spos] = self.Heap[spos], self.Heap[fpos]
  
    def minHeapify(self, pos):
        if not self.isLeaf(pos):
            if (self.Heap[pos] > self.Heap[self.leftChild(pos)] or 
               self.Heap[pos] > self.Heap[self.rightChild(pos)]):
                if self.Heap[self.leftChild(pos)] < self.Heap[self.rightChild(pos)]:
                    self.swap(pos, self.leftChild(pos))
                    self.minHeapify(self.leftChild(pos))
                else:
                    self.swap(pos, self.rightChild(pos))
                    self.minHeapify(self.rightChild(pos))
  
    def insert(self, element):
        if self.size >= self.maxsize :
            return
        self.size+= 1
        self.Heap[self.size] = element
        current = self.size
        pos = current
        while self.Heap[current] < self.Heap[self.parent(current)]:
            self.swap(current, self.parent(current))
            current = self.parent(current)
        return current
  
    def Print(self):
        for i in range(1, (self.size//2)+1):
            print(" PARENT : "+ str(self.Heap[i])+" LEFT CHILD : "+ 
                                str(self.Heap[2 * i])+" RIGHT CHILD : "+
                                str(self.Heap[2 * i + 1]))
  
    def minHeap(self):
        for pos in range(self.size//2, 0, -1):
            self.minHeapify(pos)
    
    def remove(self):
        popped = self.Heap[self.FRONT]
        self.Heap[self.FRONT] = self.Heap[self.size]
        self.size-= 1
        self.minHeapify(self.FRONT)
        return popped
  

minHeap = MinHeap(15)

minHeap.Print()
print("The Min val is " + str(minHeap.remove()))

class Voxel:
    def __init__(self, posInLattice, state):
        self.state = state
        self.neighbors = []
        self.posInLattice = posInLattice
        self.posInHeap = 0 

    def setPosInHeap(self, posInHeap):
        self.posInHeap = posInHeap

    def setNeighbors(self, neighbors):
        self.neighbors = neighbors
    
    def setValue(self, state):
        self.state = state

class DistanceField:
    def __init__(self, nx, ny, numParticles, resolution):
        self.FROZEN = 0
        self.NARROWBAND = 1
        self.COMPUTED = 2

        self.lattice = [[Voxel([i,j], self.NARROWBAND) for j in range(int(ny/resolution))] for i in range(int(nx/resolution))]
        self.I = [v in self.lattice[i] for i in range(len(self.lattice))]
        self.I[nx/resolution][ny/resolution].setValue(self.FROZEN)
        self.voxels = []
        self.neighbors = {}

        self.heap = MinHeap(numParticles)
        for v in self.I:
            self.findNeighbors(v)
            for vn in self.neighbors[v]:
                d = self.computeDistance()
                if vn.state != self.NARROWBAND:
                    vn.changeState = self.NARROWBAND
                    posInHeap = self.heap.insert(d)
                    v.setPosInHeap(posInHeap)
                else:
                    self.heap.Heap[]

    def findNeighbors(self, voxel):
        neighbors = []
        if voxel.pos[0]+1 < len(self.lattice[0]):
            neighbors.append([voxel.pos[0]+1, voxel.pos[1]])
        if voxel.pos[0]-1 >= 0:
            neighbors.append([voxel.pos[0]-1, voxel.pos[1]])
        if voxel.pos[1]+1 < len(self.lattice):
            neighbors.append([voxel.pos[0], voxel.pos[1]+1])
        if voxel.pos[1]-1 >= 0:
            neighbors.append([voxel.pos[1]-1, voxel.pos[1]])
        return neighbors

    def computeDistance(self):
        pass

