#
# License-Identifier: GPL
#
# Copyright (C) 2020 The Yambo Team
#
# Authors (see AUTHORS file for details): AM
#
# Variable definitions
#
I_PRECMP=
ifneq (,$(findstring p2y,$(MAKECMDGOALS)))
 I_PRECMP=$(p2ycpp)
endif
#
GOALS=a2y c2y
ifeq ($(do_p2y),yes)
 GOALS+=p2y
endif
ifeq ($(do_e2y),yes)
 GOALS+=e2y
endif
#
# Compilation
#
$(GOALS):
	@rm -f ${compdir}/log/"compile_"$@".log"
	@rm -f ${compdir}/config/stamps_and_lists/compilation_stop_$@.stamp
	@touch ${compdir}/config/stamps_and_lists/compiling_$@.stamp
	@$(MAKE) $(MAKEFLAGS) dependencies
	@$(MAKE) $(MAKEFLAGS) ext-libs
	@$(MAKE) $(MAKEFLAGS) int-libs
	@+LIBS="$(2Y_LIBS)";LAB="_Y_";BASE="src";ADF="$(I_PRECMP) -D_$@";$(todo_lib);$(mk_lib)
	@+LIBS="int_modules";BASE="interfaces";ADF="$(I_PRECMP)";$(todo_lib);$(mk_lib)
	@+X2DO="$@";BASE="interfaces/$@";XLIBS="$(2Y_LIBS_LD)";ADF="$(I_PRECMP)";$(todo_driver)
	@+X2DO="$@";BASE="interfaces/$@";XLIBS="$(2Y_LIBS_LD)";ADF="$(I_PRECMP)";$(mk_exe)
