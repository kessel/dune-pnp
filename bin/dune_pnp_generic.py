
import sys, os

name="dune_pnp"
solvers=["BCGS_SSORk", "BCGS_NOPREC", "CG_NOPREC", "CG_Jacobi", "CG_AMG_SSOR" ]

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-s", "--linearsolver", dest="linearsolver", default="BCGS_SSORk", 
                          help="choose linear solver", metavar="LINEARSOLVER")
parser.add_option("-p", "--pdegree", dest="pdegree", default=1,
                          help="choose degree of ansatz polynomials", metavar="PDEGREE")

(options, args) = parser.parse_args()


for i in range(len(solvers)):
    if options.linearsolver==solvers[i]:
        solverid=i+1

progname=name+"_"+solvers[solverid-1]+"_"+str(options.pdegree)

print "you want me to make ", progname

pathname = os.path.dirname(sys.argv[0])        
print 'full path =', os.path.abspath(pathname)
oldpath=os.path.abspath(os.path.curdir)

os.chdir(pathname+"/../src/")
os.system("make "+progname)
os.chdir(oldpath)
os.system(pathname+"/../src/"+progname+" "+args[0])


