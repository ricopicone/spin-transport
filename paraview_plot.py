from paraview.simple import *

files = [
	'spin_transport_soln/rho_1.pvd',
	'spin_transport_soln/rho_2.pvd',
	'spin_transport_soln/rho_3.pvd']

readers = [PVDReader(FileName = f) for f in files]

UpdatePipeline()

lines = [PlotOverLine(reader) for reader in readers]

[Show(line) for line in lines]

Render()

input('press enter to exit')
