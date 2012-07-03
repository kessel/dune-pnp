#!/usr/bin/python

import sys, os

name="dune_pnp"
solvers=["BCGS_SSORk", "BCGS_NOPREC", "CG_NOPREC", "CG_Jacobi", "CG_AMG_SSOR" ]

from optparse import OptionParser

usage = """usage: %prog [options] [ configfile ]
If no configuration file is given, only the binary corresponding to the given options is created."""
parser = OptionParser(usage=usage)

parser.add_option("-s", "--linearsolver", dest="linearsolver", default="BCGS_SSORk", 
                          help="choose linear solver", metavar="LINEARSOLVER")
parser.add_option("-p", "--pdegree", dest="pdegree", default=1,
                          help="choose degree of ansatz polynomials", metavar="PDEGREE")

parser.add_option("-n", "--np", dest="np", default=1,
                          help="number of processors", metavar="NP")

(options, args) = parser.parse_args()


for i in range(len(solvers)):
    if options.linearsolver==solvers[i]:
        solverid=i+1

progname=name+"_"+solvers[solverid-1]+"_"+str(options.pdegree)

print "If necessary, now we compile", progname

pathname = os.path.dirname(sys.argv[0])        
oldpath=os.path.abspath(os.path.curdir)

os.chdir(pathname+"/../src/")
os.system("make "+progname+" 1>out 2>err")
print "Done."
os.chdir(oldpath)
if len(args) > 0:
    os.system("mpirun -np "+str(options.np)+" " + pathname+"/../src/"+progname+" "+args[0])


