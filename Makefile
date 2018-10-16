# Main Makefile for Project Wannier Tools

include make.inc

default : wt.x

wt.x :
	@echo "Begin compiling $@"
	@if test -d src ; then \
	( cd src ; $(MAKE) || exit 1) ; fi

.PHONY : clean
clean :
	@for dir in src ; do \
		if test -d $$dir ; then \
			( cd $$dir ; $(MAKE) clean ) \
		fi \
	done
	@- rm -rf bin/*.x
