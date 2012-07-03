
SUBDIRS =

noinst_PROGRAMS = dune_pnp

dune_pnp_SOURCES = sysparams.cc dune_pnp.cc pnp_solver_main.cc GridLearning.cc

dune_pnp_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(UG_CPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	$(ALBERTA_CPPFLAGS) \
	$(ALUGRID_CPPFLAGS)
# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
dune_pnp_LDADD = \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) \
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	$(UG_LDFLAGS) $(UG_LIBS) \
	$(DUNEMPILIBS)	\
	$(LDADD)
dune_pnp_LDFLAGS = $(AM_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(UG_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	$(ALBERTA_LDFLAGS) \
	$(ALUGRID_LDFLAGS) \
	$(DUNE_LDFLAGS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules

#define BCGS_SSORk    1
#define BCGS_NOPREC   2
#define CG_NOPREC     3
#define CG_Jacobi     4
#define CG_AMG_SSOR   5
        
EXTRA_PROGRAMS= dune_pnp_BCGS_SSORk_1 dune_pnp_BCGS_SSORk_2 dune_pnp_BCGS_SSORk_3 dune_pnp_BCGS_NOPREC_1 dune_pnp_BCGS_NOPREC_2 dune_pnp_BCGS_NOPREC_3 dune_pnp_CG_NOPREC_1 dune_pnp_CG_NOPREC_2 dune_pnp_CG_NOPREC_3 dune_pnp_CG_Jacobi_1 dune_pnp_CG_Jacobi_2 dune_pnp_CG_Jacobi_3 dune_pnp_CG_AMG_SSOR_1 dune_pnp_CG_AMG_SSOR_2 dune_pnp_CG_AMG_SSOR_3 



dune_pnp_BCGS_SSORk_1_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_BCGS_SSORk_1_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=1 -DLINEARSOLVER=1
dune_pnp_BCGS_SSORk_1_LDFLAGS =$(dune_pnp_LDFLAGS)
        
dune_pnp_BCGS_SSORk_2_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_BCGS_SSORk_2_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=2 -DLINEARSOLVER=1
dune_pnp_BCGS_SSORk_2_LDFLAGS =$(dune_pnp_LDFLAGS)
        
dune_pnp_BCGS_SSORk_3_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_BCGS_SSORk_3_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=3 -DLINEARSOLVER=1
dune_pnp_BCGS_SSORk_3_LDFLAGS =$(dune_pnp_LDFLAGS)
        
dune_pnp_BCGS_NOPREC_1_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_BCGS_NOPREC_1_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=1 -DLINEARSOLVER=2
dune_pnp_BCGS_NOPREC_1_LDFLAGS =$(dune_pnp_LDFLAGS)
        
dune_pnp_BCGS_NOPREC_2_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_BCGS_NOPREC_2_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=2 -DLINEARSOLVER=2
dune_pnp_BCGS_NOPREC_2_LDFLAGS =$(dune_pnp_LDFLAGS)
        
dune_pnp_BCGS_NOPREC_3_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_BCGS_NOPREC_3_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=3 -DLINEARSOLVER=2
dune_pnp_BCGS_NOPREC_3_LDFLAGS =$(dune_pnp_LDFLAGS)
        
dune_pnp_CG_NOPREC_1_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_CG_NOPREC_1_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=1 -DLINEARSOLVER=3
dune_pnp_CG_NOPREC_1_LDFLAGS =$(dune_pnp_LDFLAGS)
        
dune_pnp_CG_NOPREC_2_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_CG_NOPREC_2_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=2 -DLINEARSOLVER=3
dune_pnp_CG_NOPREC_2_LDFLAGS =$(dune_pnp_LDFLAGS)
        
dune_pnp_CG_NOPREC_3_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_CG_NOPREC_3_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=3 -DLINEARSOLVER=3
dune_pnp_CG_NOPREC_3_LDFLAGS =$(dune_pnp_LDFLAGS)
        
dune_pnp_CG_Jacobi_1_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_CG_Jacobi_1_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=1 -DLINEARSOLVER=4
dune_pnp_CG_Jacobi_1_LDFLAGS =$(dune_pnp_LDFLAGS)
        
dune_pnp_CG_Jacobi_2_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_CG_Jacobi_2_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=2 -DLINEARSOLVER=4
dune_pnp_CG_Jacobi_2_LDFLAGS =$(dune_pnp_LDFLAGS)
        
dune_pnp_CG_Jacobi_3_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_CG_Jacobi_3_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=3 -DLINEARSOLVER=4
dune_pnp_CG_Jacobi_3_LDFLAGS =$(dune_pnp_LDFLAGS)
        
dune_pnp_CG_AMG_SSOR_1_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_CG_AMG_SSOR_1_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=1 -DLINEARSOLVER=5
dune_pnp_CG_AMG_SSOR_1_LDFLAGS =$(dune_pnp_LDFLAGS)
        
dune_pnp_CG_AMG_SSOR_2_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_CG_AMG_SSOR_2_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=2 -DLINEARSOLVER=5
dune_pnp_CG_AMG_SSOR_2_LDFLAGS =$(dune_pnp_LDFLAGS)
        
dune_pnp_CG_AMG_SSOR_3_SOURCES =$(dune_pnp_SOURCES)
dune_pnp_CG_AMG_SSOR_3_CPPFLAGS=$(dune_pnp_CPPFLAGS) -DPDEGREE=3 -DLINEARSOLVER=5
dune_pnp_CG_AMG_SSOR_3_LDFLAGS =$(dune_pnp_LDFLAGS)
        