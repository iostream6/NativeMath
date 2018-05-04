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
CND_CONF=GCC_OPENBLAS
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/e7bc4a80/arpack.o \
	${OBJECTDIR}/_ext/e7bc4a80/jni.o \
	${OBJECTDIR}/_ext/1e2cce5f/matvec.o \
	${OBJECTDIR}/_ext/1e2cce5f/xblas.o \
	${OBJECTDIR}/_ext/2abcb726/mmio.o \
	${OBJECTDIR}/_ext/2abcb726/xalloc.o \
	${OBJECTDIR}/_ext/2abcb726/xmmio.o \
	${OBJECTDIR}/_ext/2abcb726/zrng.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libNativeMath.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libNativeMath.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.c} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libNativeMath.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -Wl,--version-script=exportmap.inc -lm -lopenblas -larpack64 -shared -fPIC

${OBJECTDIR}/_ext/e7bc4a80/arpack.o: ../../src/drivers/arpack.c
	${MKDIR} -p ${OBJECTDIR}/_ext/e7bc4a80
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Wall -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/e7bc4a80/arpack.o ../../src/drivers/arpack.c

${OBJECTDIR}/_ext/e7bc4a80/jni.o: ../../src/drivers/jni.c
	${MKDIR} -p ${OBJECTDIR}/_ext/e7bc4a80
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Wall -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/e7bc4a80/jni.o ../../src/drivers/jni.c

${OBJECTDIR}/_ext/1e2cce5f/matvec.o: ../../src/matvec/matvec.c
	${MKDIR} -p ${OBJECTDIR}/_ext/1e2cce5f
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Wall -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1e2cce5f/matvec.o ../../src/matvec/matvec.c

${OBJECTDIR}/_ext/1e2cce5f/xblas.o: ../../src/matvec/xblas.c
	${MKDIR} -p ${OBJECTDIR}/_ext/1e2cce5f
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Wall -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1e2cce5f/xblas.o ../../src/matvec/xblas.c

${OBJECTDIR}/_ext/2abcb726/mmio.o: ../../src/utils/mmio.c
	${MKDIR} -p ${OBJECTDIR}/_ext/2abcb726
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Wall -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2abcb726/mmio.o ../../src/utils/mmio.c

${OBJECTDIR}/_ext/2abcb726/xalloc.o: ../../src/utils/xalloc.c
	${MKDIR} -p ${OBJECTDIR}/_ext/2abcb726
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Wall -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2abcb726/xalloc.o ../../src/utils/xalloc.c

${OBJECTDIR}/_ext/2abcb726/xmmio.o: ../../src/utils/xmmio.c
	${MKDIR} -p ${OBJECTDIR}/_ext/2abcb726
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Wall -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2abcb726/xmmio.o ../../src/utils/xmmio.c

${OBJECTDIR}/_ext/2abcb726/zrng.o: ../../src/utils/zrng.c
	${MKDIR} -p ${OBJECTDIR}/_ext/2abcb726
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Wall -D__SYSMATH_OPENBLAS__ -I${MUNOOR_SMATH}/openblas/include -I${JDK_HOME}/include -I${JDK_HOME}/include/linux -std=c11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2abcb726/zrng.o ../../src/utils/zrng.c

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
