a2y: ext-libs
	@+LIBS2DO="$(INT_LIBS)"; DIR2GO="lib" ; VPATH="$(topdir)/lib" ; $(mk_internal_lib)
	@+LIBS2DO="$(YLIBDRIVER)"; NAME="a2y_driver_"; DIR2GO="lib/yambo/driver/src"; VPATH="$(topdir)/lib/yambo/driver/src"; ADF="-D_a2y"; $(mk_external_yambo_lib)
	@+LIBS2DO="$(2YLIBS)"; XPATH="src" ; VPATH="$(topdir)/src" ; $(mk_src)
	@+LIBS2DO="int_modules"; DIR2GO="interfaces" ; VPATH="$(topdir)/interfaces" ; $(mk_internal_lib)
	@+X2DO="a2y"; XPATH="interfaces/a2y"; VPATH="$(topdir)/interfaces/a2y"; XLIBS="$(2YLIBS_LD)"; DRILIBS="$(YLIBDRIVER_LD)"; $(mk_yambo)
c2y: ext-libs
	@+LIBS2DO="$(INT_LIBS)"; DIR2GO="lib" ; VPATH="$(topdir)/lib" ; $(mk_internal_lib)
	@+LIBS2DO="$(YLIBDRIVER)"; NAME="c2y_driver_"; DIR2GO="lib/yambo/driver/src"; VPATH="$(topdir)/lib/yambo/driver/src"; ADF="-D_c2y"; $(mk_external_yambo_lib)
	@+LIBS2DO="$(2YLIBS)"; XPATH="src" ; VPATH="$(topdir)/src" ; $(mk_src)
	@+LIBS2DO="int_modules"; DIR2GO="interfaces" ; VPATH="$(topdir)/interfaces" ; $(mk_internal_lib)
	@+X2DO="c2y"; XPATH="interfaces/c2y"; VPATH="$(topdir)/interfaces/c2y"; XLIBS="$(2YLIBS_LD)"; DRILIBS="$(YLIBDRIVER_LD)"; $(mk_yambo)
ifeq ($(do_p2y),yes)
p2y: ext-libs
	@+LIBS2DO="$(INT_LIBS)"; DIR2GO="lib" ; VPATH="$(topdir)/lib" ; $(mk_internal_lib)
	@+LIBS2DO="$(YLIBDRIVER)"; NAME="p2y_driver_"; DIR2GO="lib/yambo/driver/src"; VPATH="$(topdir)/lib/yambo/driver/src"; ADF="-D_p2y"; $(mk_external_yambo_lib)
	@+LIBS2DO="$(2YLIBS)"; XPATH="src" ; VPATH="$(topdir)/src" ; $(mk_src)
	@+LIBS2DO="int_modules"; DIR2GO="interfaces" ; VPATH="$(topdir)/interfaces" ; $(mk_internal_lib)
	@+X2DO="p2y" ; XPATH="interfaces/p2y"; VPATH="$(topdir)/interfaces/p2y"; XLIBS="$(2YLIBS_LD)"; DRILIBS="$(YLIBDRIVER_LD)"; ADF="$(p2ycpp)"; $(mk_yambo) ;
endif
ifeq ($(do_e2y),yes)
e2y: ext-libs
	@+LIBS2DO="$(INT_LIBS)"; DIR2GO="lib" ; VPATH="$(topdir)/lib" ; $(mk_internal_lib)
	@+LIBS2DO="$(YLIBDRIVER)"; NAME="e2y_driver_"; DIR2GO="lib/yambo/driver/src"; VPATH="$(topdir)/lib/yambo/driver/src"; ADF="-D_e2y"; $(mk_external_yambo_lib)
	@+LIBS2DO="$(2YLIBS)"; XPATH="src" ; VPATH="$(topdir)/src" ; $(mk_src)
	@+LIBS2DO="int_modules"; DIR2GO="interfaces" ; VPATH="$(topdir)/interfaces" ; $(mk_internal_lib)
	@+X2DO="e2y" ; XPATH="interfaces/e2y"; VPATH="$(topdir)/interfaces/e2y"; XLIBS="$(2YLIBS_LD)"; DRILIBS="$(YLIBDRIVER_LD)"; $(mk_yambo) ;
endif
