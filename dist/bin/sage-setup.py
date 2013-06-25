#!/usr/bin/env python

# this program is intended to be called instead of "python setup.py"
# see #14804

# FIXME: too ugly

import distutils.file_util, distutils.dir_util, os
import distutils.command.install
#import distutils.util
distcopy = distutils.file_util.copy_file
distmkdir = distutils.dir_util.mkpath
distbytecompile = distutils.util.byte_compile

destdir = ""

try:
	destdir = os.environ["DESTDIR"]
except:
	pass

def ownmkdir(*args, **kwds):
	if args[0][0] == "/":
		args = ( destdir + args[0], )
	return distmkdir(*args, **kwds)

def ownbytecompile(*args, **kwds):
	newargs=[]
	for i in args[0]:
		if i[0] == "/":
			i = destdir + i
		newargs.append(i)
	return distbytecompile(newargs, **kwds)

def myowncopy(*args, **kwds):
	# preserve_mode=1, preserve_times=1, update=0, link=None, verbose=1, dry_run=0
	update = kwds.get("update",0)
	if args[0][0] == "/":
		args = ( destdir + args[0], args[1] )
	if args[1][0] == "/":
		args = ( args[0], destdir + args[1] )

	a = distcopy(*args, **kwds)
	return a

distutils.file_util.copy_file = myowncopy
distutils.dir_util.mkpath = ownmkdir
distutils.util.byte_compile = ownbytecompile

import sys
sys.argv[0] = "setup.py"
c = os.getcwd()
s = os.path.join(c, "setup.py")
sys.path.append(os.getcwd())
__file__ = c + "/setup.py"
execfile("setup.py")
