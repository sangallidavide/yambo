# Download external libs #
define get_external_libraries
 mkdir -p lib/archive; \
 cd lib/archive; \
 if [ "$$LIB2DO" = "all" ] ; then \
  $(MAKE) -f Makefile.loc all $(MAKEFLAGS) ; \
 else \
  $(MAKE) -f Makefile.loc $$LIB2DO $(MAKEFLAGS) ; \
 fi
endef
