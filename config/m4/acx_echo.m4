# AC_PROG_ECHO
#
# Inspired from https://invisible-island.net/autoconf/portability-echo.html
# -------------
# Check whether to use -e or not to intepret \backslash commands
#
AC_DEFUN([AC_PROG_ECHO],
[
ECHO=" "
case `sh -c 'echo "\tx"'` in
*t*)
  ECHO='echo -e';;
*)
  ECHO='echo';;
esac
AC_SUBST(ECHO)
])
