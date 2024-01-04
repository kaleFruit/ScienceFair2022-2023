from compliantMech.PRBMoptmizationnlopt import Mesh

new = Mesh(2,2)
new.connect(0,1)
new.connect(2,3)
new.connect(0,3)
new.connect(3,1)
new.setNodeType(0,1)
new.setNodeType(1,1)
# new.draw()
print(new.optimization())
# new.createChr()
# print(new.chr)
# print(new.chrLength())
# new.looseLink()
# print(new.detectDisconnect(0))
# new.redudantLink()˙
# new.draw()


