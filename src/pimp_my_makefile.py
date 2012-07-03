
import sys
solvers=["BCGS_SSORk", "BCGS_NOPREC", "CG_NOPREC", "CG_Jacobi", "CG_AMG_SSOR" ]

sys.stdout.write("EXTRA_PROGRAMS= ")
for i in range(len(solvers)):
    for j in range(3):
        sys.stdout.write("dune_pnp_"+solvers[i]+"_"+str(j)+ " ")

sys.stdout.write("\n\n\n")

for i in range(len(solvers)):
    for j in range(1,4):
        sys.stdout.write("""
dune_pnp_"""+solvers[i]+"""_"""+str(j)+"""_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_"""+solvers[i]+"""_"""+str(j)+"""_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE="""+str(j)+""" -DLINEARSOLVER="""+str(i+1)+"""
dune_pnp_"""+solvers[i]+"""_"""+str(j)+"""_LDFLAGS =$(dune_pnp_LDFLAGS)
        """)



