# Process this file with autoconf to produce a configure script.
# SRE, Wed Jul 19 11:57:36 2000
# CVS $Id: configure.ac,v 1.2 2003/10/07 17:00:17 jason Exp $
#
# Hint: See autoconf macro archive at 
#     http://www.gnu.org/software/ac-archive/


# AC_INIT information.
#
# The four AC_INIT args set the following output variables and preprocessor symbols:
#     PACKAGE_NAME      <package>     e.g. "SQUID"
#     PACKAGE_VERSION   <version>     e.g. "1.9a"
#     PACKAGE_BUGREPORT <bug-report>  e.g. "eddy@genetics.wustl.edu"
#     PACKAGE_TARNAME   <tarname>     e.g. "squid"
# From them, it derives one more:
#     PACKAGE_STRING    <package> <version>, e.g. "SQUID 1.9a"
# We define some of our own:
#     SQUID_RELCODE     unique package/version CVS tag without .'s: e.g. "squid1_9a"
#     SQUID_DATE        release date: e.g. "January 2003"
#     SQUID_COPYRIGHT   one-line copyright string
#     SQUID_LICENSE     one-line license string
#     SQUID_LICENSETAG  which license to bundle from Licenses/ subdirectory.
#
# Because SQUID is designed to be a library and must coexist with
# packages that include it as a subdirectory (HMMER, INFERNAL...),
# we don't want to use AC_INIT's PACKAGE_ variables as preprocessor
# symbols to put version info into executables; we'll get name clashes,
# plus we might want to have both squid version info and main package
# version info. So, we use the PACKAGE_ stuff to make a
# new preprocessor symbol of our own:
#     SQUID_VERSION    e.g. "1.9a"
#
# The preprocessor symbols (version info) goes into squidconf.h.
# Output variables are used in Makefiles.
#
AC_PREREQ(2.57)
AC_INIT(SQUID, 1.9g, eddy@genetics.wustl.edu, squid)
AC_MSG_NOTICE([configuring the SQUID library for your system.])

SQUID_RELCODE="squid1_9g"
SQUID_DATE="January 2003"
SQUID_COPYRIGHT="Copyright (C) 1992-2003 HHMI/Washington University School of Medicine"
SQUID_LICENSE="Freely distributed under the GNU General Public License (GPL)"
SQUID_LICENSETAG=gnu
SQUID_VERSION=$PACKAGE_VERSION

# Make output variables.
AC_SUBST(SQUID_RELCODE)
AC_SUBST(SQUID_DATE)
AC_SUBST(SQUID_COPYRIGHT)
AC_SUBST(SQUID_LICENSE)
AC_SUBST(SQUID_LICENSETAG)
AC_SUBST(SQUID_VERSION)

# Make preprocessor symbols.
AC_DEFINE_UNQUOTED(SQUID_DATE, "$SQUID_DATE")
AC_DEFINE_UNQUOTED(SQUID_COPYRIGHT, "$SQUID_COPYRIGHT")
AC_DEFINE_UNQUOTED(SQUID_LICENSE, "$SQUID_LICENSE")
AC_DEFINE_UNQUOTED(SQUID_VERSION, "$SQUID_VERSION")


# Checks for programs.
#
AC_PROG_CC
AC_PROG_RANLIB
AC_PATH_PROG([AR], [ar], [:], [$PATH:/usr/ccs/bin:/usr/xpg4/bin])

# GNU make check; format of makefile dependency lines for executables.
# original from John Darrington <j.darrington@elvis.murdoch.edu.au> 
# w/ heavy modifications.
#
# We need this because GNU make and SYSV make use different systems
# specifying variables for dependencies: $$@ in sysv, %: %.o in GNU.
# Would love to hear a better way of doing this.
# 
# I use two different conventions in my Makefiles. Sometimes 
# executable "foo" has a file "foo.c" - this is the HMMER convention.
# Sometimes executable "foo" has a file "foo_main.c" - this is
# the SQUID convention. The configure script sets the
# EXEC_DEPENDENCY appropriately: here, HMMER style.
#
# This creates a function CHECK_GNU_MAKE, which we immediately call.
# It sets an output variable EXEC_DEPENDENCY. 
# This is used in the src/Makefile.in.
#
AC_DEFUN(CHECK_GNU_MAKE,[ 
  AC_MSG_CHECKING(whether your make is GNU make)
  foundGNUmake='nope, assuming sysv make.' ;
  EXEC_DEPENDENCY=[\$\$\@_main.o] ;
  if ( make --version nothing 2> /dev/null | grep GNU > /dev/null ) ;  then
     foundGNUmake='yes, it is.' ;
     EXEC_DEPENDENCY='%: %_main.o' ;
  fi
  AC_MSG_RESULT($foundGNUmake)
  AC_SUBST(EXEC_DEPENDENCY)
])
CHECK_GNU_MAKE

	
# ================================================================
# Provide for unsigned integers of known size
# Sets SQD_UINT16, 32, and 64.
# Will substitute for @SQD_UINT16@, etc. in an output file (config.h)
# Substitutes "FIXME" if no appropriate typedef is found.
# This info is substituted in squid.h.
# SRE, Thu Dec 28 13:58:51 2000
#
AC_C_BIGENDIAN()
AC_CHECK_FUNCS(ntohs, , AC_CHECK_LIB(socket, ntohs))
AC_CHECK_FUNCS(ntohl, , AC_CHECK_LIB(socket, ntohl))
AC_CHECK_FUNCS(htons, , AC_CHECK_LIB(socket, htons))
AC_CHECK_FUNCS(htonl, , AC_CHECK_LIB(socket, htonl))

AC_CHECK_SIZEOF(unsigned short, 2)
AC_CHECK_SIZEOF(unsigned int, 4)
AC_CHECK_SIZEOF(unsigned long, 4)
AC_CHECK_SIZEOF(unsigned long long, 8)
AC_CHECK_FUNCS(strtoul strtoull)
if test "$ac_cv_sizeof_unsigned_short" = "2"; then
  SQD_UINT16="unsigned short    "
else
  SQD_UINT16="FIXME"
  AC_MSG_WARN("No 16-bit int? Manually edit config file to typedef sqd_uint16.")
fi
if test "$ac_cv_sizeof_unsigned_int" = "4"; then
  SQD_UINT32="unsigned int      "
elif test "$ac_cv_sizeof_unsigned_long" = "4"; then
  SQD_UINT32="unsigned long     "
else
  SQD_UINT32="FIXME"
  AC_MSG_WARN("No 32-bit int? Manually edit config file to typedef sqd_uint32.")
fi
if test "$ac_cv_sizeof_unsigned_long" = "8"; then
  SQD_UINT64="unsigned long     "
elif test "$ac_cv_sizeof_unsigned_long_long" = "8"; then
  SQD_UINT64="unsigned long long"
else
  SQD_UINT64="FIXME"
  AC_MSG_WARN("No 64-bit int? Manually edit config file to typedef sqd_uint64.")
fi
AC_SUBST(SQD_UINT16)
AC_SUBST(SQD_UINT32)
AC_SUBST(SQD_UINT64)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# ================================================================
# Test for whether we can cheat and treat fpos_t
# (used by fgetpos() and fsetpos()) as an arithmetic datatype.
# This is essential if we're to be able to save an fpos_t to a file
# in an architecture-neutral format. We need this on FreeBSD
# systems, which don't provide a 64-bit ftell64() or ftello(),
# so we have to use fgetpos().
#
# if true, squidconf.h will contain
#   #define ARITHMETIC_FPOS_T 1
# else if false, 
#   #undef ARITHMETIC_FPOS_T
AC_DEFUN(SQ_ARITHMETIC_FPOS_T, [
   AC_MSG_CHECKING(whether fpos_t is an arithmetic datatype)
   fpos_arithmetic="no."
   AC_TRY_COMPILE([#include <stdio.h>],
      [int main(void) { fpos_t f1, f2; if (f1 == f2) f1 = 0;}],
      [AC_DEFINE(ARITHMETIC_FPOS_T) 
       fpos_arithmetic="yes."])
   AC_MSG_RESULT($fpos_arithmetic)
])
SQ_ARITHMETIC_FPOS_T()
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ================================================================
# Tests for 64-bit file offset functions
# Note: only AC_CHECK_FUNCS seems to properly define HAVE_FOO?
#
AC_CHECK_FUNCS(strtoull)
AC_CHECK_FUNCS(ftello fseeko)
AC_CHECK_FUNCS(ftello64 fseeko64)
AC_CHECK_FUNCS(ftell64 fseek64)
AC_CHECK_FUNC(stat64, 
   [AC_MSG_CHECKING(for struct stat64)
    stat64_struct="no!"	    
    AC_TRY_COMPILE([#include <sys/types.h>
                    #include <sys/stat.h>
                    #include <unistd.h>],
                   [int main(void) { struct stat64 s1;} ],
	           [AC_DEFINE(HAVE_STAT64)
                    stat64_struct="yes."])
    AC_MSG_RESULT($stat64_struct)]
)
AC_CHECK_SIZEOF(off_t)
AC_CHECK_SIZEOF(off64_t)
AC_CHECK_SIZEOF(fpos_t)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#################################################################################
# Now we're into our "optional features"
#################################################################################


# --enable-ccmalloc      - configure for compiling in ccmalloc memory checks
#
# The CC output variable is modified.
# Requires ccmalloc version >= 0.4.0, because it uses --no-wrapper.
#
AC_ARG_ENABLE(ccmalloc,
[  --enable-ccmalloc       turn ccmalloc memory debugging on ],
[case $enable_ccmalloc in
   yes)  AC_MSG_NOTICE([enabled ccmalloc memory debugging])
         CC="ccmalloc --no-wrapper $CC" 
         ;;
    no)  AC_MSG_NOTICE([ccmalloc memory debug/tracing disabled])
         ;;
     *)  echo "Ignoring unknown argument to --enable-ccmalloc: $enable_ccmalloc"
	 ;;
esac])



# --enable-debugging=x    - set debugging level to <x> (1-3)
# 
# At all levels, including 0, replaces CFLAGS w/ "-g -Wall" (so it assumes gcc).
# Sets the DEBUGLEVEL preprocessor symbol to <x>
#
AC_ARG_ENABLE(debugging,
[  --enable-debugging      set CFLAGS for debugging
  --enable-debugging=x    also set diagnostics level to <x> (1-3) (3 = most verbose)],
[ case $enable_debugging in
   yes)  AC_MSG_NOTICE([enabled debugging diagnostics level 0 (CFLAGS only, no verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(DEBUGLEVEL, 0)
         ;;
     1)  AC_MSG_NOTICE([enabled debugging diagnostics level 1 (low verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(DEBUGLEVEL, 1)
         ;;
     2)  AC_MSG_NOTICE([enabled debugging diagnostics level 2 (moderate verbosity)])
         CFLAGS="-g -Wall"
	 AC_DEFINE(DEBUGLEVEL, 2)
         ;;
     3)  AC_MSG_NOTICE([enabled debugging diagnostics level 3 (high verbosity)])
         CFLAGS="-g -Wall"
	 AC_DEFINE(DEBUGLEVEL, 3)
         ;;
    no)  AC_MSG_NOTICE([debugging diagnostics disabled])
         AC_DEFINE(DEBUGLEVEL, 0)
         ;;
     *)  echo "Ignoring unknown argument to --enable-debugging: $enable_debugging"
	 ;;
esac])



# --enable-lfs           Large File Summit (LFS) support for >2GB files
# See: http://ftp.sas.com/standards/large.file/x_open.20Mar96.html
#
AC_ARG_ENABLE(lfs, 
[  --enable-lfs            enable LFS, Large File Support],
[case $enable_lfs in
   yes) AC_MSG_NOTICE([configured for optional LFS, large file support])
	AC_DEFINE(_LARGEFILE_SOURCE)
	AC_DEFINE(_LARGEFILE64_SOURCE)
	AC_DEFINE(_FILE_OFFSET_BITS, 64)
        ;;
   no)  ;;
    *)  echo "Ignoring unknown argument to --enable-lfs: $enable_lfs"
        ;;
esac])



# --enable-pvm       Enable Parallel Virtual Machine (PVM) support
# 
# Sets PVMLIBDIR, PVMINCDIR, PCMPROGS, PVMLIBS output variables
# Sets SRE_ENABLE_PVM preprocessor variable.
#
AC_ARG_ENABLE(pvm, 
[  --enable-pvm            enable PVM, Parallel Virtual Machine],
[case $enable_pvm in
   yes)  AC_MSG_NOTICE([enabled optional PVM (Parallel Virtual Machine) support])
 	 PVMLIBDIR="-L${PVM_ROOT}/lib/${PVM_ARCH}"
	 PVMINCDIR="-I${PVM_ROOT}/include"
	 PVMLIBS="-lpvm3"
	 AC_DEFINE(SRE_ENABLE_PVM)
	 ;;
   no)   AC_MSG_NOTICE([PVM (Parallel Virtual machine) support disabled])
         ;;
   *)    echo "Ignoring unknown argument to --enable-pvm: $enable_pvm"
	 ;;
esac])
AC_SUBST(PVMLIBDIR)
AC_SUBST(PVMINCDIR)
AC_SUBST(PVMLIBS)


# Write out squidconf.h header
AC_CONFIG_HEADER(squidconf.h)

# Write out Makefiles, and squid.h also has some output variable substitution.
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(Testsuite/Makefile)
AC_CONFIG_FILES(squid.h)

AC_OUTPUT


