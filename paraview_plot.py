import sys
from paraview.simple import *

reader = PVDReader(FileName = sys.argv[1])

UpdatePipeline()

PlotOverLine(reader)

Show()
Render()

avifile = None
mp4file = None
try:
	avifile = sys.argv[2] + '.avi'
	mp4file = sys.argv[2] + '.mp4'
except:
	pass

AnimateReader(reader, filename = avifile)

if avifile is not None:
	import subprocess, os
	subprocess.call('mencoder -speed 45 -o {} -ovc lavc {}'.format(mp4file, avifile).split(' '))
	os.remove(avifile)
