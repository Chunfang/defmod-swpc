# Makefile for Defmod

FFLAGS =
FPPFLAGS =
LDLIBS =

-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/conf/rules
-include ${PETSC_DIR}/lib/petsc/conf/variables
-include ${PETSC_DIR}/lib/petsc/conf/rules

ifeq ($(debug),true)
	FFLAGS += -O0 -g
else
	ifeq ($(CC),icc)
		FFLAGS += -ipo -diag-disable ipo -static_intel
	else
		FFLAGS += -Wno-maybe-uninitialized -Wno-tabs -Wno-unused-function
	endif
endif

OBJS = main.o m_esh3d.o m_h5io.o m_fvfe.o m_galpha.o m_fefd.o m_global.o m_local.o m_elems.o m_utils.o

m_utils.o  : m_utils.F90
m_esh3d.o  : m_esh3d.F90 m_utils.o
m_elems.o  : m_elems.F90 m_utils.o
m_local.o  : m_local.F90 m_elems.o m_esh3d.o
m_global.o : m_global.F90 m_local.o
m_fefd.o   : m_fefd.F90 m_global.o
m_galpha.o : m_galpha.F90 m_global.o
m_fvfe.o   : m_fvfe.F90 m_global.o
m_h5io.o   : m_h5io.F90 m_global.o
main.o     : main.F90 m_global.o m_fefd.o m_galpha.o m_fvfe.o m_h5io.o

all: ${OBJS}
	-${FLINKER} ${OBJS} -o ../../bin/defmod ${PETSC_LIB} ${LDLIBS}
