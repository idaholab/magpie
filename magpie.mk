# Check the existence of the contrib submodules and build accordingly
SPPARKS_DIR    ?= $(APPLICATION_DIR)/contrib/spparks
ifneq ($(wildcard $(SPPARKS_DIR)/src/Makefile),)
  ADDITIONAL_CPPFLAGS += -DSPPARKS_ENABLED
  app_INCLUDES   += -I $(SPPARKS_DIR)/..
  include $(APPLICATION_DIR)/contrib/spparks.mk
endif

