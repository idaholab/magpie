# Check the existence of the contrib submodules and build accordingly
SPPARKS_SUBMODULE    := $(CURDIR)/spparks
ifneq ($(wildcard $(SPPARKS_SUBMODULE)/src/Makefile),)
  include spparks.mk
endif
