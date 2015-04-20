# Check the existence of the contrib submodules and build accordingly
SPPARKS_DIR    ?= $(CURDIR)/spparks
ifneq ($(wildcard $(SPPARKS_DIR)/src/Makefile),)
	libmesh_CXXFLAGS += -DSPPARKS_ENABLED
  include spparks.mk
endif
