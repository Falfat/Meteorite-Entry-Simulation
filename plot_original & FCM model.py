import meteor
import meteor_FCM
import matplotlib.pyplot as plt
plt.style.use('ggplot')

#original model
mathilde1 = meteor.meteor(v=10E3, m=1200E3, theta=20, r=30, z0=100E3, g=9.81, Rp=6400E3, sigma0=64*1E3, pa_tab=False, planet="Earth")
print('The result of original model is:')
mathilde1.results()

print('\n')

#FCM model
mathilde2 = meteor_FCM.meteor(v=10E3, m=1200E3, theta=20, r=30, g0=9.81, P=64*1E3,Rp=6400E3)
print('The result of alternative model(Fragment-cloud Model) is:')
mathilde2.results()

#read the data from two models
data1 = mathilde1.fire()
energydata1 = mathilde1.energypatch()
data2 = mathilde2.fire()
energydata2 = mathilde2.energypatch()

#plot datas
plt.subplot(1, 1, 1)
plt.plot(data1[0], data1[1][3],'b-',data2[0],data2[1][3], 'r--')
plt.legend(('$the\ original\ model$','$Fragment-cloud\ model$'),loc='best')
plt.title('Plot of Altitude vs Time')
plt.ylabel('Altitude (m)')
plt.xlabel('Time (s)')
plt.xticks(size='small')
plt.yticks(size='small')

plt.tight_layout()
plt.show()

plt.subplot(1, 1, 1)
plt.plot(data1[1][3], data1[1][0], 'b-',data2[1][3],data2[1][0],'r--')
plt.legend(('$the\ original\ model$','$Fragment-cloud\ model$'),loc='best')
plt.title('Plot of Velocity vs Altitude')
plt.xlabel('Altitude(m)')
plt.ylabel('Velocity (m/s)')
plt.xticks(size='small')
plt.yticks(size='small')
plt.gca().invert_xaxis()

plt.tight_layout()
plt.show()

plt.subplot(1, 1, 1)
plt.plot(data1[1][3], data1[1][1], 'b-',data2[1][3],data2[1][1],'r--')
plt.legend(('$the\ original\ model$','$Fragment-cloud\ model$'),loc='best')
plt.title('Plot of Mass vs Altitude')
plt.xlabel('Altitude(m)')
plt.ylabel('Mass (kg)')
plt.xticks(size='small')
plt.yticks(size='small')
plt.gca().invert_xaxis()

plt.tight_layout()
plt.show()

plt.subplot(1, 1, 1)
plt.plot(data1[1][3], data1[1][5], 'b-',data2[1][3],data2[1][4],'r--')
plt.legend(('$the\ original\ model$','$Fragment-cloud\ model$'),loc='best')
plt.title('Plot of Radius vs Altitude')
plt.xlabel('Altitude(m)')
plt.ylabel('Radius (m)')
plt.xticks(size='small')
plt.yticks(size='small')
plt.gca().invert_xaxis()

plt.tight_layout()
plt.show()

plt.subplot(1, 1, 1)
plt.plot(energydata1[1], energydata1[2],'b-',energydata2[1],energydata2[2],'r--')
plt.legend(('$the\ original\ model$','$Fragment-cloud\ model$'),loc='best')
plt.title('Plot of Altitude vs Energy Release')
plt.xlabel('Energy Release (Kt/Km)')
plt.ylabel('Altitude(m)')

plt.tight_layout()
plt.show()
