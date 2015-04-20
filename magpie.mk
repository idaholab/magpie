# add contrib to include path
app_INCLUDES   += -I $(APPLICATION_DIR)contrib

# Check the existence of the contrib submodules and build accordingly
SPPARKS_DIR    ?= $(APPLICATION_DIR)contrib/spparks
ifneq ($(wildcard $(SPPARKS_DIR)/src/Makefile),)
  EXTRA_FLAGS += -DSPPARKS_ENABLED
  include $(APPLICATION_DIR)contrib/spparks.mk
endif

