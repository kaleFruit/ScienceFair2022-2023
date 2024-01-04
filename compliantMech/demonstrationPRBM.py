import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from sympy import *
from sympy.vector import CoordSys3D
from numpy import *
from functools import partial
from scipy import optimize
from scipy.optimize import NonlinearConstraint
from scipy.optimize import LinearConstraint


class Mesh:
    def __init__(
        self,
        nX,
        nY,
        scale,
        chr,
        outputNodes,
        constrainedNodes,
        inputNodes,
        freedomOutput,
        freedomInput,
        k,
        force,
    ):
        self.N = CoordSys3D("N")
        self.k = k
        self.forceInputNodes = inputNodes
        self.force = force

        self.freedomOutput = freedomOutput
        self.freedomInput = freedomInput

        self.scale = scale
        self.nNum = nX * nY
        self.nX = nX
        self.nY = nY
        self.constrainedNodes = []
        self.outputNodes = outputNodes
        self.M = np.zeros((self.nNum, self.nNum))
        self.chr = chr
        if not chr:
            self.createChr()
        else:
            self.createFromChr()
        self.links = sum(self.chr)
        for i in constrainedNodes:
            self.setNodeType(i, 1)
        xs = [x for x in symbols(f"x0:{self.nNum}")]
        ys = [y for y in symbols(f"y0:{self.nNum}")]
        self.nodePos = list(zip(xs, ys))
        self.initNodePos = [
            (x * scale, y * scale) for y in range(self.nY) for x in range(self.nX)
        ]
        self.edgeLabels = {}

    def draw(self, newCoors=[]):
        G = nx.Graph()
        color_map = []
        node_sizes = []
        for nodeID, nPos in enumerate(self.initNodePos):
            G.add_node(nodeID, pos=nPos)
            color_map.append("blue")
            node_sizes.append(1000)
        for nodeID, nPos in enumerate(newCoors):
            if sum(self.M[nodeID]) >= 1:
                G.add_node(float(nodeID) + 0.1, pos=nPos)
                color_map.append("green")
                node_sizes.append(500)
        for i in range(len(self.M)):
            for j in range(len(self.M[i])):
                if self.M[i][j] == 1:
                    G.add_edge(i, j)
        if len(newCoors) != 0:
            for i in range(len(self.M)):
                for j in range(len(self.M[i])):
                    if self.M[i][j] == 1:
                        G.add_edge(float(i) + 0.1, float(j) + 0.1)
        pos = nx.get_node_attributes(G, "pos")
        nx.draw(
            G,
            pos,
            with_labels=True,
            node_size=node_sizes,
            node_color=color_map,
            connectionstyle="arc3, rad=0.12",
            arrows=True,
        )
        nx.draw_networkx_edge_labels(
            G, pos, edge_labels=self.edgeLabels, font_color="red"
        ),
        plt.show()

    def newDraw(self, newCoors):
        G = nx.Graph()
        node_sizes = []
        color_map = []
        for nodeID, nPos in enumerate(newCoors):
            if sum(self.M[nodeID]) >= 1:
                G.add_node(float(nodeID) + 0.1, pos=nPos)
                color_map.append("green")
                node_sizes.append(500)
        for i in range(len(self.M)):
            for j in range(len(self.M[i])):
                if self.M[i][j] == 1:
                    G.add_edge(float(i) + 0.1, float(j) + 0.1)
        pos = nx.get_node_attributes(G, "pos")
        nx.draw(
            G,
            pos,
            with_labels=True,
            node_size=node_sizes,
            node_color=color_map,
            connectionstyle="arc3, rad=0.12",
            arrows=True,
        )
        plt.show()

    def createChr(self):
        chr = []
        for i in range(self.nNum - 1):
            for j in range(1 + i, self.nNum):
                chr.append(self.M[i][j])
        self.chr = chr
        self.links = sum(chr)
        return chr

    def createFromChr(self):
        counter = 0
        for row in range(self.nNum - 1):
            for column in range(row + 1, self.nNum):
                self.M[row][column] = self.chr[counter]
                self.M[column][row] = self.chr[counter]
                counter += 1

    def setNodeType(self, node, type):
        if type == 1:
            self.M[node][node] = 1
            self.constrainedNodes.append(node)
        elif type == 0:
            self.M[node][node] = 0
            self.constrainedNodes.remove(node)

    def connect(self, node1, node2):
        self.M[node1][node2] = 1
        self.M[node2][node1] = 1

    def disconnect(self, node1, node2):
        self.M[node1][node2] = 0
        self.M[node2][node1] = 0

    def chrLength(self):
        return 0.5 * self.nNum * (self.nNum - 1)

    def looseLink(self):
        tisfine = self.constrainedNodes + self.forceInputNodes + self.outputNodes
        for _ in range(2):
            for i in range(len(self.M)):
                if sum(self.M[i]) == 1:
                    for j in range(len(self.M[i])):
                        if self.M[i][j] == 1 and (i not in tisfine):
                            self.disconnect(i, j)

    def redudantLink(self):
        tisfine = self.constrainedNodes
        for i in range(self.nNum):
            if sum(self.M[i]) == 2 and self.M[i][i] not in tisfine:
                connections = []
                cut = False
                for j in range(self.nNum):
                    if self.M[i][j] == 1:
                        connections.append(j)
                for j in connections:
                    if cut:
                        break
                    for _ in [
                        x
                        for x in range(self.nNum)
                        if x in connections and self.M[j][x] == 1
                    ]:
                        self.disconnect(i, connections[0])
                        self.disconnect(i, connections[1])
                        cut = True

    def detectDisconnect(self, start):
        isConnected = self.DFS(start, [False] * self.nNum, [], start)
        while not isConnected:
            currConnected = []
            for i in range(self.nNum - 1):
                for j in range(i + 1, self.nNum):
                    if self.M[i][j] == 1:
                        currConnected.append(i)
                        currConnected.append(j)
            if currConnected:
                for input in self.forceInputNodes:
                    self.connect(input, random.choice(currConnected))
                for output in self.outputNodes:
                    self.connect(output, random.choice(currConnected))
                for constr in self.constrainedNodes:
                    self.connect(constr, random.choice(currConnected))
                isConnected = self.DFS(start, [False] * self.nNum, [], start)

    def DFS(self, start, visited, constrained, OG):
        myVisited = 0
        visited[start] = True
        if (
            start in self.constrainedNodes
            or start in self.forceInputNodes
            or start in self.outputNodes
        ):
            constrained.append(start)
        for i in range(self.nNum):
            if self.M[start][i] == 1 and (not visited[i]):
                self.DFS(i, visited, constrained, OG)
                myVisited += 1
        if start == OG:
            return set(constrained) == set(
                self.constrainedNodes + self.forceInputNodes + self.outputNodes
            )

    def unconstrainedAngle(self, j, i, k):
        curr = self.lawOfCosines(self.nodePos[j], self.nodePos[i], self.nodePos[k])
        initial = self.lawOfCosines(
            self.initNodePos[j], self.initNodePos[i], self.initNodePos[k]
        )
        return curr - initial

    def lawOfCosines(self, jPoint, iPoint, kPoint):
        lIJ = ((jPoint[0] - iPoint[0]) ** 2 + (jPoint[1] - iPoint[1]) ** 2) ** 0.5
        lIK = ((kPoint[0] - iPoint[0]) ** 2 + (kPoint[1] - iPoint[1]) ** 2) ** 0.5
        lKJ = ((kPoint[0] - jPoint[0]) ** 2 + (kPoint[1] - jPoint[1]) ** 2) ** 0.5
        return acos((lIJ**2 + lIK**2 - lKJ**2) / (2 * lIJ * lIK))

    def constrainedAngle(self, i, j):
        initial = atan2(
            (self.initNodePos[j][1] - self.initNodePos[i][1]),
            (self.initNodePos[j][0] - self.initNodePos[i][0]),
        )
        curr = atan2(
            self.nodePos[j][1] - self.nodePos[i][1],
            (self.nodePos[j][0] - self.nodePos[i][0]),
        )
        return curr - initial

    def angularDeformation(self):
        dA = []
        angleLabels = []
        for i in range(self.nNum):
            constrained = False
            if i in self.constrainedNodes:
                constrained = True
            neighbors = []
            for j in range(self.nNum):
                if j == i:
                    pass
                elif self.M[i][j] == 1:
                    if constrained:
                        dA.append(self.constrainedAngle(i, j))
                        angleLabels.append((i, j))
                    else:
                        neighbors.append(j)
            if neighbors:
                j = neighbors[0]
                og = neighbors[0]
                passCheck = True
                if len(neighbors) >= 2:
                    passCheck = False
                while len(neighbors) >= 2:
                    neighbors.remove(j)
                    potentialAngles = {}
                    for k in neighbors:
                        potentialAngles[k] = self.lawOfCosines(
                            self.initNodePos[j],
                            self.initNodePos[i],
                            self.initNodePos[k],
                        )
                    smallestAngleMate = min(potentialAngles, key=potentialAngles.get)
                    dA.append(self.unconstrainedAngle(j, i, smallestAngleMate))
                    angleLabels.append((j, i, smallestAngleMate))
                    j = smallestAngleMate
                if not passCheck:
                    latestAngle = self.lawOfCosines(
                        self.initNodePos[angleLabels[-1][0]],
                        self.initNodePos[angleLabels[-1][1]],
                        self.initNodePos[angleLabels[-1][2]],
                    )

                    if (
                        ogCheck := self.lawOfCosines(
                            self.initNodePos[og],
                            self.initNodePos[i],
                            self.initNodePos[j],
                        )
                    ) < latestAngle:
                        dA[-1] = ogCheck
                        angleLabels[-1] = (og, i, j)
        return Matrix(dA)

    def length(self, iPoint, jPoint):
        return ((jPoint[1] - iPoint[1]) ** 2 + (jPoint[0] - iPoint[0]) ** 2) ** 0.5

    def rigidConstraints(self):
        cB = []
        for i in range(self.nNum):
            for j in range(i):
                if self.M[i][j] == 1:
                    cB.append(
                        self.length(self.nodePos[i], self.nodePos[j])
                        - self.length(self.initNodePos[i], self.initNodePos[j])
                    )
                    self.edgeLabels[(i, j)] = self.length(
                        self.nodePos[i], self.nodePos[j]
                    ) - self.length(self.initNodePos[i], self.initNodePos[j])
        print(f"cB: {cB}")
        return cB

    def constrainedPoints(self):
        cR = []
        for i in self.constrainedNodes:
            cR.append(self.nodePos[i][0] - self.initNodePos[i][0])
            cR.append(self.nodePos[i][1] - self.initNodePos[i][1])
        for i in self.outputNodes:
            if not self.freedomOutput[0]:
                cR.append(self.nodePos[i][0] - self.initNodePos[i][0])
            if not self.freedomOutput[1]:
                cR.append(self.nodePos[i][1] - self.initNodePos[i][1])
        for i in self.forceInputNodes:
            if not self.freedomInput[0]:
                cR.append(self.nodePos[i][0] - self.initNodePos[i][0])
            if not self.freedomInput[1]:
                cR.append(self.nodePos[i][1] - self.initNodePos[i][1])
        print(f"cR: {cR}")
        return cR

    def staticAnalysis(self):
        coors = np.array(self.nodePos).flatten()
        dA = self.angularDeformation()
        print(f"da: {dA}")

        def objective(x):
            values = dict(zip(coors, x))
            displacement = (
                self.nodePos[self.outputNodes[0]][0]
                - self.initNodePos[self.outputNodes[0]][0]
            )
            f = 0.5 * self.k * dA.T.dot(dA) - (self.force * (displacement))
            val = double(f.evalf(subs=values))
            return val

        def constraint_eq(x, c):
            values = dict(zip(coors, x))
            val = double(c.evalf(subs=values))
            return val

        def linearConstraintEQ():
            constraintDict = []
            for c in self.constrainedPoints():
                matrix = np.zeros(len(coors))
                c = str(c)
                base = 0
                if c[0] == "x":
                    matrix[int(c[1]) * 2] = 1
                    base = self.initNodePos[int(c[1])][0]
                else:
                    matrix[int(c[1]) * 2 + 1] = 1
                    base = self.initNodePos[int(c[1])][1]
                constraintDict.append(LinearConstraint(matrix, base, base))
            return constraintDict

        xInit = np.array(self.initNodePos).flatten().astype(double)
        xInit += 0
        constraintDict = []
        for c in self.rigidConstraints():
            constraintDict.append(
                NonlinearConstraint(partial(constraint_eq, c=c), 0, 0)
            )

        for linear in linearConstraintEQ():
            constraintDict.append(linear)
        res = optimize.minimize(
            objective,
            xInit,
            method="trust-constr",
            constraints=constraintDict,
            options={"xtol": 1e-3, "maxiter": 10, "verbose": 2},
        )
        print(f"Optimal solution; x = {res}")
        return res.x

    def results(self, xopt):
        newCoors = xopt.reshape(np.array(self.nodePos).shape)
        coors = np.array(self.nodePos).flatten()
        values = dict(zip(coors, xopt))
        replaced = self.angularDeformation().evalf(subs=values)
        for i, replace in enumerate(replaced):
            print(replace)
            if isinstance(replace, complex):
                print(self.angularDeformation()[i])
        maxDeformation = sum(self.angularDeformation().evalf(subs=values))
        displacement = newCoors[1][0] - self.initNodePos[1][0]
        return (displacement, maxDeformation)


candidate = Mesh(
    nX=4,
    nY=3,
    scale=10,
    chr=None,
    outputNodes=[3],
    constrainedNodes=[0],
    inputNodes=[8],
    freedomOutput=[True, False],
    freedomInput=[False, True],
    k=573,
    force=9.8 * 1000,
)
candidate.connect(0, 5)
candidate.connect(5, 8)
candidate.connect(8, 3)
candidate.createChr() 
candidate.draw()

analysis = candidate.staticAnalysis()
newCoors = analysis.reshape(np.array(candidate.nodePos).shape)
candidate.draw(newCoors)
candidate.newDraw(newCoors)