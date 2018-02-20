import numpy as np
import os
import matplotlib
if 'DISPLAY' not in os.environ:
	matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def update_plot(num, lines, time, rho1, rho2, rho3, t, x):
	lines[0][0].set_data(x, rho1[num,:])
	lines[1][0].set_data(x, rho2[num,:])
	lines[2][0].set_data(x, rho3[num,:])
	time.set_text('t={:.3f} s'.format(t[num]))

def plot(
	rho1 = None,
	rho2 = None,
	rho3 = None,
	t = None,
	x = None,
	title = 'Spin Transport',
	interval = 100,
	filename = None):

	if rho1 is None:
		raise ValueError('rho1 not defined')
	if rho2 is None:
		raise ValueError('rho2 not defined')
	if rho3 is None:
		raise ValueError('rho3 not defined')
	if t is None:
		raise ValueError('t not defined')
	if x is None:
		raise ValueError('x not defined')

	fig1 = plt.figure()

	lines = [
		plt.plot([], [], 'r-', label = r'$\rho_1$'),
		plt.plot([], [], 'b-', label = r'$\rho_2$'),
		plt.plot([], [], 'g-', label = r'$\rho_3$')
	]

	plt.legend()

	plt.xlim(0, 1)
	plt.ylim(
		np.min(np.hstack([rho1, rho2, rho3])),
		np.max(np.hstack([rho1, rho2, rho3])))

	plt.xlabel('$x$')

	plt.title(title)

	time_label = plt.text(0.95, 0.05, 't=0 s',
		horizontalalignment='right',
		verticalalignment='baseline',
		transform=fig1.axes[0].transAxes)

	line_ani = animation.FuncAnimation(fig1,
		update_plot,
		t.shape[0],
		fargs=(lines, time_label, rho1, rho2, rho3, t, x),
		interval = interval)

	if filename is not None:
		extension = filename.split('.')[-1]
		if extension == 'mp4':
			Writer = animation.writers['ffmpeg']
			writer = Writer(fps=15, metadata={}, bitrate=1800)
		elif extension == 'gif':
			Writer = animation.writers['imagemagick']
			writer = Writer(fps=15)
		elif extension == 'html':
			from JSAnimation import HTMLWriter
			writer = HTMLWriter()
		else:
			raise ValueError('"{}" is not a known file extension'.format(filename.split('.')[-1]))
		line_ani.save(filename, writer = writer)
	else:
		plt.show()

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(
		description = 'Spin Transport Solution Plotter',
		add_help = True)
	parser.add_argument('-f',
		type = str,
		help = 'The file to plot.',
		default = 'spin_transport_soln/soln.npz')
	parser.add_argument('-s',
		type = str,
		help = 'Filename to save an animation of the data to.')
	parser.add_argument('-i',
		type = int,
		help = 'How long to wait between frames in ms.',
		default = 100)
	parser.add_argument('-t',
		type = str,
		help = 'The plot title to use.',
		default = 'Spin Transport')
	args = parser.parse_args()
	data = np.load(args.f)
	plot(
		rho1 = data['rho1'],
		rho2 = data['rho2'],
		rho3 = data['rho3'],
		t = data['t'],
		x = data['x'],
		filename = args.s,
		interval = args.i,
		title = args.t)
