import sys
from paraview.simple import *

reader = PVDReader(FileName = sys.argv[1])

UpdatePipeline()

PlotOverLine(reader)

Show()
Render()

outfile = None
try:
	outfile = sys.argv[2]
except:
	pass
AnimateReader(reader, filename = '{}.avi'.format(outfile))

if outfile is not None:
	import subprocess, os
	subprocess.call('mencoder -speed 45 -o {0}.mp4 -ovc lavc {0}.avi'.format(outfile).split(' '))
	os.remove('{}.avi'.format(outfile))
