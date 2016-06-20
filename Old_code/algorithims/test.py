#import add
#print add.zadd([1,2,3],[4,5,6])


from numpy import f2py
with open("addd.f") as sourcefile:
	sourcecode=sourcefile.read()
	f2py.compile(sourcecode,modulename='hello')

import hello
print hello.zadd([1,2,3],[4,5,6])
