import pymysql as sql
import numpy as np
import json

def upload(
	host = 'sql.camerondevine.me',
	username = 'SpinTransSim',
	password = '',
	database = 'MRFM',
	table = 'spin_transport',
	port = 3306,
	hash = None,
	rho1 = None,
	rho2 = None,
	rho3 = None,
	t = None,
	x = None):

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

	if hash is None:
		import git
		hash = git.Repo().head.object.hexsha
	elif len(hash) != 40:
		raise ValueError('hash is not of correct length (40)')

	if type(rho1) == np.ndarray:
		rho1 = json.dumps(rho1.tolist())
	if type(rho2) == np.ndarray:
		rho2 = json.dumps(rho2.tolist())
	if type(rho3) == np.ndarray:
		rho3 = json.dumps(rho3.tolist())
	if type(t) == np.ndarray:
		t = json.dumps(t.tolist())
	if type(x) == np.ndarray:
		x = json.dumps(x.tolist())

	connection = sql.connect(
		host = host,
		user = username,
		password = password,
		port = port,
		database = database)
	cursor = connection.cursor()
	cursor.execute(
		"INSERT INTO `{}` (`hash`, `rho1` , `rho2`, `rho3`, `t`, `x`) VALUES (%s, %s, %s, %s, %s, %s);".format(table.replace(';', '')),
		(
			hash,
			rho1,
			rho2,
			rho3,
			t,
			x))
	cursor.close()
	connection.commit()
	connection.close()

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(
		description = 'Spin Transport SQL Uploader',
		add_help = False)
	parser.add_argument('-f',
		type = str,
		help = 'The file to upload.',
		default = 'spin_transport_soln/soln.npz')
	parser.add_argument('-h',
		type = str,
		help = 'The SQL server to upload results to.',
		default = 'sql.camerondevine.me')
	parser.add_argument('-u',
		type = str,
		help = 'The SQL server username.',
		default = 'SpinTransSim')
	parser.add_argument('-p',
		type = str,
		help = 'The SQL server password.',
		default = '')
	parser.add_argument('--hash',
		type = str,
		help = 'The commit hash to from which the results were generated')
	parser.add_argument('--port',
		type = int,
		help = 'The SQL server port',
		default = 3306)
	parser.add_argument('-t',
		type = str,
		help = 'The SQL table to upload to',
		default = 'spin_transport')
	parser.add_argument('-d',
		type = str,
		help = 'The SQL database to uplaod to',
		default = 'MRFM')
	args = parser.parse_args()
	data = np.load(args.f)
	upload(
		password = args.p,
		rho1 = data['rho1'],
		rho2 = data['rho2'],
		rho3 = data['rho3'],
		t = data['t'],
		x = data['x'],
		username = args.u,
		table = args.t,
		database = args.d,
		hash = args.hash,
		port = args.port)
