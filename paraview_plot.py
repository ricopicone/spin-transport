import sys
from paraview.simple import *

reader = PVDReader(FileName = sys.argv[1])

UpdatePipeline()

PlotOverLine(reader)

Show()
Render()
AnimateReader(reader)
