import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

L = 1

pvd_files = ['spin_transport_soln/rho_{}.pvd'.format(i) for i in range(1,4)]

pvd_datasets = [ET.parse(pvd_file).getroot().find('Collection') for pvd_file in pvd_files]

def getVtuData(filename):
	vtu_data = ET.parse(filename).getroot()[0][0].find('PointData')[0].text.strip().split('  ')
	return [float(point) for point in vtu_data]

t = np.ndarray([len(pvd_datasets[0])])
rho = np.ndarray([len(pvd_datasets), len(pvd_datasets[0]), 1])

for i in range(len(pvd_datasets[0])):
	t[i] = float(pvd_datasets[0][i].attrib['timestep'])
	for j in range(len(pvd_datasets)):
		data = getVtuData('spin_transport_soln/' + pvd_datasets[j][i].attrib['file'])
		if len(data) != rho.shape[2]:
			rho.resize([len(pvd_datasets), len(pvd_datasets[0]), len(data)])
		rho[j,i,:] = data

x = range

def update_plot(num, lines, time):
	for i in range(len(lines)):
		lines[i][0].set_data(np.linspace(0, L, rho.shape[2]), rho[i,num,:])
	time.set_text('t={:.3f} s'.format(t[num]))
	return lines[0]

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata={}, bitrate=1800)

fig1 = plt.figure()

lines = [
	plt.plot([], [], 'r-', label = r'$\rho_1$'),
	plt.plot([], [], 'b-', label = r'$\rho_2$'),
	plt.plot([], [], 'g-', label = r'$\rho_3$')
]

plt.legend()

plt.xlim(0, 1)
plt.ylim(np.min(rho), np.max(rho))

plt.xlabel('$x$')

plt.title('Spin Transport')

time_label = plt.text(0.95, 0.05, 't=0 s',
	horizontalalignment='right',
	verticalalignment='baseline',
	transform=fig1.axes[0].transAxes)

line_ani = animation.FuncAnimation(fig1,
	update_plot,
	t.shape[0],
	fargs=(lines, time_label),
	interval = 100,
	blit = True)

line_ani.save('spin_transport_soln.mp4', writer = writer)
