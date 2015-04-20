# add contrib to include path
app_INCLUDES   += -I $(CURDIR)/contrib

# Check the existence of the contrib submodules and build accordingly
SPPARKS_DIR    ?= $(CURDIR)/contrib/spparks
ifneq ($(wildcard $(SPPARKS_DIR)/src/Makefile),)
	libmesh_CXXFLAGS += -DSPPARKS_ENABLED
  include contrib/spparks.mk
endif
