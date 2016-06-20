from numpy import f2py
with open("rungeadd.f") as sourcefile:
	sourcecode=sourcefile.read()
	f2py.compile(sourcecode,modulename='coefs')

