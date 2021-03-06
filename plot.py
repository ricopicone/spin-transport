import numpy as np
import os
import matplotlib
if 'DISPLAY' not in os.environ:
	matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def update_plot(num, lines, time, rho1, rho2, rho3, t, x, n):
	lines[0][0].set_data(x, rho1[num,:])
	lines[1][0].set_data(x, rho2[num,:])
	lines[2][0].set_data(x, rho3[num,:])
	time.set_text(('t={:.' + str(n) + 'f} s').format(t[num]))

def plot(
	rho1 = None,
	rho2 = None,
	rho3 = None,
	t = None,
	x = None,
	title = 'Spin Transport',
	interval = 100,
	filename = None,
	frames = [0, -1],
	time_p = 3,
	subplots = False):

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

	rho1 = rho1[frames[0]:frames[1],:]
	rho2 = rho2[frames[0]:frames[1],:]
	rho3 = rho3[frames[0]:frames[1],:]
	t = t[frames[0]:frames[1]]

	fig1 = plt.figure()

	axes = [plt] * 3
	if subplots:
		axes = [
			plt.subplot(3, 1, 1),
			plt.subplot(3, 1, 2),
			plt.subplot(3, 1, 3)]

	lines = [
		axes[0].plot([], [], 'r-', label = r'$\rho_1$'),
		axes[1].plot([], [], 'b-', label = r'$\rho_2$'),
		axes[2].plot([], [], 'g-', label = r'$\rho_3$')
	]

	for axis in axes:
		axis.legend()

	xlim = [np.min(x), np.max(x)]

	if subplots:
		axes[0].set_title(title)
		for (axis, flow) in zip(axes, [rho1, rho2, rho3]):
			axis.set_xlim(*xlim)
			ylim = [np.min(flow), np.max(flow)]
			if ylim[0] - ylim[1] == 0:
				ylim = [ylim[0] - 1e-9, ylim[1] + 1e-9]
			axis.set_ylim(*ylim)
			axis.set_xlabel('$x$')
	else:
		plt.title(title)
		plt.xlim(*xlim)
		plt.ylim(
			np.min(np.hstack([rho1, rho2, rho3])),
			np.max(np.hstack([rho1, rho2, rho3])))
		plt.xlabel('$x$')

	text = [0.95, 0.05, 't=0 s']
	if subplots:
		text[1] = -2.35

	time_label = axes[2].text(*text,
		horizontalalignment='right',
		verticalalignment='baseline',
		transform=fig1.axes[0].transAxes)

	line_ani = animation.FuncAnimation(fig1,
		update_plot,
		t.shape[0],
		fargs=(lines, time_label, rho1, rho2, rho3, t, x, time_p),
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
	parser.add_argument('--start',
		type = int,
		default = 0,
		help = 'The frame to start the animation at.')
	parser.add_argument('--end',
		type = int,
		default = -1,
		help = 'The ending frame of the animation.')
	parser.add_argument('-n',
		type = int,
		default = 3,
		help = 'The time precision to use.')
	parser.add_argument('-p',
		action = 'store_true',
		help = 'Enable subplots')
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
		title = args.t,
		frames = [args.start, args.end],
		time_p = args.n,
		subplots = args.p)
