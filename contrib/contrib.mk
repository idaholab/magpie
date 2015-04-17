# Check the existence of the contrib submodules and build accordingly
SPPARKS_DIR    ?= $(CURDIR)/spparks
ifneq ($(wildcard $(SPPARKS_DIR)/src/Makefile),)
  include spparks.mk
endif
