#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
m4_include([m4/ax_check_zlib.m4])
m4_include([m4/ax_cxx_compile_stdcxx.m4])
m4_include([m4/ax_compiler_version.m4])
m4_include([m4/ax_gcc_x86_cpu_supports.m4])
AC_PREREQ([2.67])
AC_INIT([kmc], [3.7.5], [gtg240@gmail.com])
AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([1.11 subdir-objects foreign])
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
# Checks for programs.
OLD_CXXFLAGS=${CXXFLAGS}
AC_PROG_CXXCPP
AC_PROG_CXX
CXXFLAGS="${OLD_CXXFLAGS} -O3"
AC_PROG_RANLIB
AC_PROG_LN_S
AC_PROG_SED
AC_PROG_AWK
[
if test "x`${CXX} --version | $SED -n -e 's/.*\Apple \(LLVM\) version.*clang.*/\1/p'`" = "xLLVM" ;
then
    gcc_version=`${CC} --version | $SED -n -e 's/.*LLVM version \([0-9.]*\).*$/\1/p'`
elif test "x`${CXX} --version | $SED -n -e 's/.*\(clang\) version.*/\1/p'`" = "xclang" ;
then
    # Don't use -dumpversion with clang: it will always produce "4.2.1".
    gcc_version=`${CXX} --version | $SED -n -e 's/.*clang version \([0-9.]*\).*$/\1/p'`
elif test "x`${CXX} --version | $SED -n -e 's/icc.*\(ICC\).*/\1/p'`" = "xICC" ; 
then
    gcc_version=`${CXX} -dumpversion 2>/dev/null`
else
    gcc_version=`${CXX} -dumpversion 2>/dev/null`
    if test "x$gcc_version" = x; then
	gcc_version=`${CXX} --version | $SED -n -e 's/[^ ]*gcc[^ ]* ([^)]*) \([0-9.]*\).*$/\1/p'`
    fi
fi
]
AX_CHECK_ZLIB
AC_CHECK_PROG([KMC_HAS_GCC], [gcc], [yes], [no])
AC_CHECK_PROG([KMC_HAS_CLANG], [clang], [yes], [no])
# Checks for header files.
AC_CHECK_HEADERS([zlib.h stdint.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
# Checks for available std
AX_CXX_COMPILE_STDCXX([14], [ext], [mandatory])
# Checks for available instructions
have_avx2=false
have_avx=false
have_sse41=false
have_neon=false
AH_TEMPLATE([HAVE_NEON_INSTRUCTIONS], [Define to 1 if NEON SIMD instructions are available on aarch64])
if test ${target_cpu} = "aarch64";
then
have_neon=true
AC_DEFINE([HAVE_NEON_INSTRUCTIONS], [1])
KMC_CXXFLAGS="-march=armv8-a+simd"
else
AX_GCC_X86_CPU_SUPPORTS([avx2], 
			[have_avx2=true
			 KMC_CXXFLAGS="-mavx2"],
			[AX_GCC_X86_CPU_SUPPORTS([avx],
						 [have_avx=true
						  KMC_CXXFLAGS="-mavx"],
						 [AX_GCC_X86_CPU_SUPPORTS([sse4.1],
									  [have_sse41=true
									   KMC_CXXFLAGS="-msse4.1"],
									  [AC_MSG_ERROR([*** Requires AVX2, AVX, or SSE4.1 instructions available on x86_64])])])])
fi
AM_CONDITIONAL([HAVE_NEON], [$have_neon])
AM_CONDITIONAL([HAVE_AVX2], [$have_avx2])
AM_CONDITIONAL([HAVE_AVX], [$have_avx])
AM_CONDITIONAL([HAVE_SSE41], [$have_sse41])
# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_MKTIME
AC_CHECK_FUNCS([gettimeofday memmove memset modf strchr strdup strtoull])
KMC_CXXFLAGS="$KMC_CXXFLAGS -fPIC -Wall -fsigned-char"
KMC_LDFLAGS="$KMC_LDFLAGS -Wl,--whole-archive -lpthread -Wl,--no-whole-archive"
AC_SUBST([KMC_CXXFLAGS])
AC_SUBST([KMC_LDFLAGS])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

echo \
"----------------------------------------------------------------

${PACKAGE_NAME} Version ${PACKAGE_VERSION}

Current compiler options:
   Target arch:                         ${target_cpu}: ${target_os}
   CXX:                                 ${CXX}
   CXX version:                         ${gcc_version}
   ranlib:                              ${RANLIB}
   Default CXXFLAGS:                    ${CXXFLAGS}
   KMC CXXFLAGS:                        ${KMC_CXXFLAGS}
   KMC LDFLAGS:                         ${KMC_LDFLAGS}
   KMC LIB:                             ${LIBS}

----------------------------------------------------------------"
