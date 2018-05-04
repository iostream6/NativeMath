#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=GCC_OPENBLAS_DEBUG
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/e7bc4a80/arpack.o \
	${OBJECTDIR}/_ext/1e2cce5f/matvec.o \
	${OBJECTDIR}/_ext/1e2cce5f/xblas.o \
	${OBJECTDIR}/_ext/2abcb726/mmio.o \
	${OBJECTDIR}/_ext/2abcb726/xalloc.o \
	${OBJECTDIR}/_ext/2abcb726/xmmio.o \
	${OBJECTDIR}/_ext/2abcb726/zrng.o \
	${OBJECTDIR}/test_mat.o \
	${OBJECTDIR}/test_xblas.o


# C Compiler Flags
CFLAGS=-fopenmp

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L${MUNOOR_SMATH} -L${MUNOOR_SMATH}/openblas/lib

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/nativemathdev

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/nativemathdev: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.c} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/nativemathdev ${OBJECTFILES} ${LDLIBSOPTIONS} -lm -lopenblas -larpack64

${OBJECTDIR}/_ext/e7bc4a80/arpack.o: ../../src/drivers/arpack.c
	${MKDIR} -p ${OBJECTDIR}/_ext/e7bc4a80
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -D__DEBUG__ -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/e7bc4a80/arpack.o ../../src/drivers/arpack.c

${OBJECTDIR}/_ext/1e2cce5f/matvec.o: ../../src/matvec/matvec.c
	${MKDIR} -p ${OBJECTDIR}/_ext/1e2cce5f
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -D__DEBUG__ -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1e2cce5f/matvec.o ../../src/matvec/matvec.c

${OBJECTDIR}/_ext/1e2cce5f/xblas.o: ../../src/matvec/xblas.c
	${MKDIR} -p ${OBJECTDIR}/_ext/1e2cce5f
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -D__DEBUG__ -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1e2cce5f/xblas.o ../../src/matvec/xblas.c

${OBJECTDIR}/_ext/2abcb726/mmio.o: ../../src/utils/mmio.c
	${MKDIR} -p ${OBJECTDIR}/_ext/2abcb726
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -D__DEBUG__ -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2abcb726/mmio.o ../../src/utils/mmio.c

${OBJECTDIR}/_ext/2abcb726/xalloc.o: ../../src/utils/xalloc.c
	${MKDIR} -p ${OBJECTDIR}/_ext/2abcb726
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -D__DEBUG__ -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2abcb726/xalloc.o ../../src/utils/xalloc.c

${OBJECTDIR}/_ext/2abcb726/xmmio.o: ../../src/utils/xmmio.c
	${MKDIR} -p ${OBJECTDIR}/_ext/2abcb726
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -D__DEBUG__ -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2abcb726/xmmio.o ../../src/utils/xmmio.c

${OBJECTDIR}/_ext/2abcb726/zrng.o: ../../src/utils/zrng.c
	${MKDIR} -p ${OBJECTDIR}/_ext/2abcb726
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -D__DEBUG__ -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2abcb726/zrng.o ../../src/utils/zrng.c

${OBJECTDIR}/test_mat.o: test_mat.c
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -D__DEBUG__ -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/test_mat.o test_mat.c

${OBJECTDIR}/test_xblas.o: test_xblas.c
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -D__DEBUG__ -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/test_xblas.o test_xblas.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
