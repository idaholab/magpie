#
# Check the existence of the contrib submodules and build accordingly
#

SPPARKS_DIR    ?= $(MAGPIE_DIR)/contrib/spparks
ifneq ($(wildcard $(SPPARKS_DIR)/src/Makefile),)
  ADDITIONAL_CPPFLAGS += -DSPPARKS_ENABLED
  app_INCLUDES   += -I $(SPPARKS_DIR)/..
  include $(MAGPIE_DIR)/contrib/spparks.mk
endif

MYTRIM_DIR    ?= $(MAGPIE_DIR)/contrib/mytrim
ifneq ($(wildcard $(MYTRIM_DIR)/trim.h),)
  ADDITIONAL_CPPFLAGS += -DMYTRIM_ENABLED
  app_INCLUDES   += -I $(MYTRIM_DIR)/..
  include $(MAGPIE_DIR)/contrib/mytrim.mk
endif

include $(MAGPIE_DIR)/contrib/gsl.mk
include $(MAGPIE_DIR)/contrib/fftw3.mk
