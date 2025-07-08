import matplotlib.pyplot as plt

X, Y = [], []
for line in open('points.txt', 'r'):
  values = [int(s) for s in line.split()]
  X.append(values[0])
  Y.append(values[1])
  
FX, FY = [], []
for line in open('filteredPoints.txt', 'r'):
  values = [int(s) for s in line.split()]
  FX.append(values[0])
  FY.append(values[1])

SX, SY = [], []
for line in open('skyline.txt', 'r'):
  values = [int(s) for s in line.split()]
  SX.append(values[0])
  SY.append(values[1])

plt.plot(X, Y, 'bo')
plt.plot(FX, FY, 'bo')
plt.plot(SX, SY, 'ro')
plt.plot(SX, SY, 'r--', label="Conventional Skyline")
plt.plot(SX, SY, 'r-', label="Linear Skyline")

# Area dominated by the first upper-bound
#plt.plot([SX[0],SX[0]], [SY[0],max(Y)], 'r--')
#plt.plot([SX[0],max(X)], [SY[0],SY[0]], 'r--')

# Area dominated by the first upper-bound
#plt.plot([SX[len(SX)-1],max(X)], [SY[len(SY)-1],SY[len(SY)-1]], 'r--')
#plt.plot([SX[len(SX)-1],SX[len(SX)-1]], [SY[len(SY)-1],max(Y)], 'r--')

axes = plt.gca()
axes.set_xlim([5,50])
axes.set_ylim([50,200])
#plt.autoscale(enable='True', axis='both')
plt.xlabel('Travel distance')
plt.ylabel('POI cost')
#plt.title('Linear Skyline')
plt.legend(loc=9, ncol=3)
plt.show()
