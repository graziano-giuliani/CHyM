#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([ICHyM], [1.0.1], [ggiulian@ictp.it])

AC_SUBST([PACKAGE_URL], ["https://github.com/graziano-giuliani/CHyM"])

LT_INIT

AM_INIT_AUTOMAKE([foreign])

AC_CONFIG_MACRO_DIRS([m4])

# Checks for programs.

AM_PROG_AR
AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_FC

AC_LANG([Fortran])

AC_ARG_WITH([netcdf],
  AC_HELP_STRING([--with-netcdf],
            [Path to NetCDF installation (default: NETCDF environment)]),
  [NC_PREFIX="$withval"; NF_PREFIX="$withval"; NFAUTO="no"; NCAUTO="no"],
  [FAKE="no"]
)

PARALLEL="1"
AC_ARG_ENABLE(mpiserial,
  AC_HELP_STRING([--enable-mpiserial],
                 [Use the included MPI replacement library for single processor]),
  [DONOTCHECKMPI="DONOTCHECKMPI"; PARALLEL="0"]
)

AC_ARG_ENABLE(nc4,
  AC_HELP_STRING([--enable-nc4],
                 [Supply this option if you plan on using NetCDF4/HDF5 option.]),
  [NETCDF4_HDF5="1"],[NETCDF4_HDF5="0"]
)
if eval "test x$NETCDF4_HDF5 = x1"; then
  AM_CPPFLAGS="\$(DEFINE)NETCDF4_HDF5 $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

NC_INCLUDES=""
NC_LDFLAGS=""
NC_LIBS=""

AX_PROG_NF_CONFIG([NF_PREFIX=`nf-config --prefix`; NFAUTO="yes"],
                  [NF_PREFIX="${NETCDF}"; NFAUTO="no"])

AX_PROG_NC_CONFIG([NC_PREFIX=`nc-config --prefix`; NCAUTO="yes"],
                  [NC_PREFIX="${NETCDF}"; NCAUTO="no"])

##old # Checks for libraries.
##old # FIXME: Replace `main' with a function in `-lnetcdf':
##old AC_CHECK_LIB([netcdf], [-lnetcdf])
##old # FIXME: Replace `main' with a function in `-lnetcdff':
##old AC_CHECK_LIB([netcdff], [-lnetcdff])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

if eval "test x$NFAUTO = xyes"; then
  if eval "test x$NCAUTO = xyes"; then
    NC_INCLUDES="`nc-config --cflags` `nf-config --fflags`"
    #NC_LDFLAGS="-L`nf-config --prefix`/lib -L`nc-config --prefix`/lib"
    NC_LIBS="`nc-config --libs` `nf-config --flibs`"
  else
    NC_INCLUDES="-I$NC_PREFIX/include `nf-config --fflags`"
    NC_LDFLAGS="-L$NC_PREFIX/lib -L`nc-config --prefix`/lib"
    NC_LIBS="-L$NC_PREFIX/lib `nf-config --flibs`"
  fi
else
  if eval "test x$NCAUTO = xno"; then
    NC_INCLUDES="-I$NC_PREFIX/include -I$NF_PREFIX/include"
    NC_LDFLAGS="-L$NC_PREFIX/lib -L$NF_PREFIX/lib"
    NC_LIBS="-lnetcdff -lnetcdf"
  else
    NCVER=`nc-config --version | cut -d " " -f 2`
    if eval "test x$NCVER = x4.0.1"; then
      NC_INCLUDES="`nc-config --fflags`"
    else
      # new netcdf library has this
      NC_INCLUDES="-I`nc-config --includedir`"
    fi
    NC_LDFLAGS="-L`nc-config --prefix`/lib"
    NC_LIBS="`nc-config --libs`"
  fi
fi

AC_SUBST(NC_PREFIX)
AC_SUBST(NF_PREFIX)
AC_SUBST(NC_INCLUDES)
AC_SUBST(NC_LDFLAGS)
AC_SUBST(NC_LIBS)

AC_LANG_PUSH([C])
AC_C_BIGENDIAN
AC_LANG_POP([C])
__AC_FC_NAME_MANGLING

case $ac_cv_fc_mangling in
  "lower case, underscore, no extra underscore")
    name_mangling="FORTRAN_UNDERSCORE_";;
  "lower case, no underscore")
    name_mangling="FORTRAN_SAME";;
  "upper case, no underscore")
    name_mangling="FORTRAN_CAPS_";;
  "lower case, double underscore")
    name_mangling="FORTRAN_GNUF2C";;
  "lower case, no underscore, no extra underscore")
    name_mangling="FORTRAN_SAME";;
  *)
    name_mangling="FORTRAN_MANGLING_UNKNOWN"
    AC_MSG_WARN([UNKNOWN FORTRAN NAME-MANGLING SCHEME]);;
esac

RR_PATH_NETCDF
RR_PATH_NETCDF_F90

test "$osname"=NONE && osname=`uname -s`
OSNAME=$osname
AC_SUBST(OSNAME)

test "$machinename"=NONE && machinename=`uname -m`
fullhostname=`hostname -f`

# DEFINE VARIABLES ACCORDING TO OS AND COMPILER

if echo $osname | grep -i aix >/dev/null 2>&1; then
   SYSDEF="AIX"
fi
if echo $osname | grep -i darwin >/dev/null 2>&1; then
   SYSDEF="DARWIN"
fi
if echo $osname | grep -i unix_system_v >/dev/null 2>&1; then
   SYSDEF="UNIXSYSTEMV"
fi
if echo $osname | grep -i irix >/dev/null 2>&1; then
   SYSDEF="IRIX"
fi
if echo $osname | grep -i irix64 >/dev/null 2>&1; then
   SYSDEF="IRIX64"
fi
if echo $osname | grep -i linux >/dev/null 2>&1; then
   SYSDEF="LINUX"
fi
if echo $osname | grep -i osf1 >/dev/null 2>&1; then
   SYSDEF="OSF1"
fi
if echo $osname | grep -i super >/dev/null 2>&1; then
   SYSDEF="SUPERUX"
fi
if echo $osname | grep -i sun >/dev/null 2>&1; then
   SYSDEF="SUNOS"
fi
if echo $osname | grep -i t3e >/dev/null 2>&1; then
   SYSDEF="T3E"
fi
if echo $osname | grep -i unicos >/dev/null 2>&1; then
   SYSDEF="UNICOS"
fi
if test -z "$SYSDEF"; then
   AC_MSG_WARN([OPERATING SYSTEM UNKNOWN])
   SYSDEF="UNKNOWNOS"
fi

AM_CONDITIONAL(OSNAME_AIX, [ test x"$osname" = xAIX ] )
AM_CONDITIONAL(OSNAME_LINUX, [ test x"$osname" = xLinux ] )

FCPROG=`basename $(echo $FC | cut -d " " -f 1)`
AC_SUBST(FCPROG)

AM_CONDITIONAL(COMPILER_IBM, [ test x"$FCPROG" = xxlf2003 -o \
                                    x"$FCPROG" = xxlf2003_r -o \
                                    x"$FCPROG" = xxlf2008 -o \
                                    x"$FCPROG" = xxlf2008_r -o \
                                    x"$FCPROG" = xbgxlf2003 -o \
                                    x"$FCPROG" = xbgxlf2003_r -o \
                                    x"$FCPROG" = xbgxlf2008 -o \
                                    x"$FCPROG" = xbgxlf2008_r ] )
AM_CONDITIONAL(COMPILER_GNU, [ test x"$FCPROG" = xgfortran ] )
AM_CONDITIONAL(COMPILER_G95, [ test x"$FCPROG" = xg95 ] )
AM_CONDITIONAL(COMPILER_SUN, [ test x"$FCPROG" = xsunf95 ] )
AM_CONDITIONAL(COMPILER_PGI, [ test x"$FCPROG" = xpgf95  ] )
AM_CONDITIONAL(COMPILER_INTEL, [ test x"$FCPROG" = xifort -o \
                                      x"$FCPROG" = xifx -o \
                                      x"$FCPROG" = xmpiifx -o \
                                      x"$FCPROG" = ximpif90 ])
AM_CONDITIONAL(CHYM_PARALLEL, [ test x"$PARALLEL" = "x1"])

# Force Big endian I/O for known compilers
if test -z "$COMPILER_GNU_TRUE"; then
#for debugging  FCFLAGS="-fconvert=big-endian -fno-range-check -finline-functions -g -fcheck=all -Wall"
  #FCFLAGS="-fconvert=big-endian -fno-range-check -finline-functions "
  FCFLAGS="-fconvert=big-endian -g -fcheck=all -Wall"
  CPRDEF="GNU"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_INTEL_TRUE"; then
  FCFLAGS="-convert big_endian -assume byterecl \
           -heap-arrays -finline-functions"
  CPRDEF="INTEL"
  AC_SUBST(FCFLAGS)
  LIBS=`echo $LIBS | sed -e 's/-Wl,--as-needed//g'`
  AC_SUBST(LIBS)
fi
if test -z "$COMPILER_G95_TRUE"; then
  FCFLAGS="-fendian=big -fmultiple-save"
  CPRDEF="GNU"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_PGI_TRUE"; then
  FCFLAGS="-byteswapio"
  CPRDEF="PGI"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_IBM_TRUE"; then
  FCFLAGS="-WF,-DIBM"
  CPRDEF="XLF"
  AC_SUBST(FCFLAGS)
fi
if test -x "$COMPILER_SUN_TRUE"; then
  CPRDEF="WORKSHOP"
fi

AC_SUBST(CPRDEF)

MCTPATH="${PWD}/external"
if test -z "$COMPILER_G95_TRUE"; then
  FCFLAGS="-O2 $FCFLAGS"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_GNU_TRUE"; then
  #FCFLAGS="-O0 -g -pg -ffpe-trap=zero,invalid,overflow,underflow -fbacktrace \
  #         -fdump-core -Wall -pedantic -Wconversion -Wsurprising -Waliasing \
  #         -Wampersand -Warray-bounds -Wcharacter-truncation -Wuninitialized \
  #         -Wline-truncation -Wintrinsics-std -Wno-tabs -Wunderflow \
  #         -Wunused-parameter -Wno-align-commons $FCFLAGS \
  #         -fno-range-check -Wimplicit-interface -Wimplicit-procedure \
  #         -Wintrinsic-shadow -fcheck=array-temps,bounds,do,mem,recursion"
  #FCFLAGS="-O2 -mtune=native $FCFLAGS"
  FCFLAGS="-O0 -fbacktrace $FCFLAGS"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_INTEL_TRUE"; then
  #FCFLAGS="-O0 -fpe0 -g -check all -traceback -warn all -ftz \
  #         -fimf-arch-consistency=true -fltconsistency -mieee-fp $FCFLAGS"
  FCFLAGS="-O2 -fp-model fast=2 $FCFLAGS"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_IBM_TRUE"; then
  FCFLAGS="-qsuppress=1501-245 -O2 -qstrict -qxlf2003=polymorphic \
           -qmaxmem=-1 -qzerosize $FCFLAGS"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_SUN_TRUE"; then
  FCFLAGS="-O2 -xfilebyteorder=big16:%all \
           -D__DATE__="'`date +%D`'" -D__TIME__="'`date +%T`'" $FCFLAGS"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_PGI_TRUE"; then
  FCFLAGS="-fastsse $FCFLAGS"
  AC_SUBST(FCFLAGS)
fi

CHYM_ROOT=${PWD}
if eval "test x"$PARALLEL" = x1"; then
  ACX_MPI
else
  FORT_SIZE="real4double8"
  MPIFC=$FC
  LIBS="$LIBS -L$CHYM_ROOT/external -lmpiser"
  AC_SUBST(MPIFC)
  AC_SUBST(FORT_SIZE)
  AM_CPPFLAGS="\$(DEFINE)MPI_SERIAL $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi
MPEUCFLAGS="\$(DEFINE)SYS$SYSDEF \$(DEFINE)CPR$CPRDEF \$(DEFINE)$name_mangling"
MPEUFCFLAGS="\$(DEFINE)SYS$SYSDEF \$(DEFINE)CPR$CPRDEF"
if test -n "$SRKDEF"; then
  MPEUCFLAGS="$MPEUCFLAGS \$(DEFINE)$SRKDEF"
  MPEUFCFLAGS="$MPEUFCFLAGS \$(DEFINE)$SRKDEF"
fi
AC_SUBST(MPEUCFLAGS)
AC_SUBST(MPEUFCFLAGS)

F90FLAGS="$FCFLAGS"
AC_SUBST(F90FLAGS)

AC_SUBST(CHYM_ROOT)
AC_PREFIX_DEFAULT(${PWD})

if test -z "$COMPILER_SUN_TRUE"; then
  CPPFLAGS="$AM_CPPFLAGS -M$CHYM_ROOT/share"
  AC_SUBST(CPPFLAGS)
else
  CPPFLAGS="$AM_CPPFLAGS -I$CHYM_ROOT/share"
  AC_SUBST(CPPFLAGS)
fi
LIBS="-L$CHYM_ROOT/share -lrchlib $LIBS"
echo $LIBS
SYSLIBS="$LIBS"
AC_SUBST(SYSLIBS)
##

##wait AC_SUBST(NC_INCLUDES)
##wait AC_SUBST(NC_LIBS)
##wait AC_SUBST(NC_LDFLAGS)

AC_CONFIG_HEADERS(config.h)
AC_CONFIG_FILES([Makefile external/Makefile share/Makefile preproc/Makefile utility/Makefile utility/dewtimeseries/Makefile utility/dewvalidation/Makefile utility/hourlydb/Makefile utility/dranetwrite/Makefile main/Makefile])
AC_OUTPUT
