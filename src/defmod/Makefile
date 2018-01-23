# Makefile for Defmod

FFLAGS =
FPPFLAGS =
LDLIBS =

-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/conf/rules
-include ${PETSC_DIR}/lib/petsc/conf/variables
-include ${PETSC_DIR}/lib/petsc/conf/rules

OBJS = main.o m_galpha.o m_fefd.o m_global.o m_local.o m_elems.o m_utils.o

m_utils.o  : m_utils.F90
m_elems.o  : m_elems.F90 m_utils.o
m_local.o  : m_local.F90 m_elems.o
m_global.o : m_global.F90 m_local.o 
m_fefd.o   : m_fefd.F90 m_global.o
m_galpha.o : m_galpha.F90 m_global.o
main.o     : main.F90 m_global.o m_fefd.o m_galpha.o

all: ${OBJS}
	-${FLINKER} ${OBJS} -o ../../bin/defmod ${PETSC_LIB} ${LDLIBS}
