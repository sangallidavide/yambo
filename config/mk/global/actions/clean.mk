clean_all: dependencies
	@$(call clean_driver,all)
clean: dependencies
	@$(call clean_driver,$(what))
reset: dependencies
	@$(call clean_driver,"")
	@$(call clean_driver,"int-libs")
check: 
	@FILES=`git ls-files --others|grep -v .tar.gz|grep -v .sw`; for target in $$FILES ; do echo $$target; done
	@find . -empty -type d
