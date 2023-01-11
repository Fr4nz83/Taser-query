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

plt.plot(X, Y, 'bo', label='Dominated')
plt.plot(FX, FY, 'bo', label='Dominated')
plt.plot(SX, SY, 'ro', label='Skyline')
plt.plot(SX, SY, 'r-')

# Area dominated by the first upper-bound
plt.plot([SX[0],SX[0]], [SY[0],max(Y)], 'r--')
plt.plot([SX[0],max(X)], [SY[0],SY[0]], 'r--')

# Area dominated by the first upper-bound
plt.plot([SX[len(SX)-1],max(X)], [SY[len(SY)-1],SY[len(SY)-1]], 'r--')
plt.plot([SX[len(SX)-1],SX[len(SX)-1]], [SY[len(SY)-1],max(Y)], 'r--')

plt.autoscale(enable='True', axis='both')
plt.xlabel('Travel distance')
plt.ylabel('POI cost')
plt.title('Linear Skyline')
plt.legend(loc=9, ncol=3)
plt.show()
