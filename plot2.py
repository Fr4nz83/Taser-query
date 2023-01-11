import matplotlib.pyplot as plt

X, Y = [9,27,13,32,31,26,22,36], [160,100,158,160,98,140,142,158]
SX, SY = [9,13,22,26,27,31], [160,158,142,140,100,98]
LSX, LSY = [9,27,31], [160,100,98]
# LSX, LSY = [13,31], [158,98]

plt.plot(X, Y, 'bo')
#plt.plot(LSX, LSY, 'ro-')
plt.plot(SX, SY, 'r--', label="Conventional Skyline")
plt.plot(LSX, LSY, 'r-', label="Linear Skyline")

# Area dominated by the first upper-bound
#plt.plot([LSX[0],LSX[0]], [LSY[0],200], 'r--')
#plt.plot([LSX[0],60], [LSY[0],LSY[0]], 'r--')

# Area dominated by the second upper-bound
#plt.plot([SX[len(SX)-1],60], [SY[len(SY)-1],SY[len(SY)-1]], 'r--')
#plt.plot([SX[len(SX)-1],SX[len(SX)-1]], [SY[len(SY)-1],200], 'r--')

#axes = plt.gca()
#axes.set_xlim([0,50])
#axes.set_ylim([60,190])
plt.autoscale(enable='True', axis='both')
plt.xlabel('Travel distance')
plt.ylabel('Travel cost')
#plt.title('Linear Skyline')
plt.legend(loc=9, ncol=2)
plt.show()
