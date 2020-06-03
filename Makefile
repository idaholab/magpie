###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
#
###############################################################################

MOOSE_SUBMODULE    := $(CURDIR)/moose
MOOSE_PARENT := $(shell dirname `pwd`)/moose

# MOOSE_DIR is empty or unset
ifeq ($(wildcard $(MOOSE_DIR)/framework/Makefile),)
  # submodule contains valid moose
  ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
    MOOSE_DIR        := $(MOOSE_SUBMODULE)
  # valid moose next to the magpie directory
  else ifneq ($(wildcard $(MOOSE_PARENT)/framework/Makefile),)
    MOOSE_DIR        := $(MOOSE_PARENT)
  else
    $(error MOOSE framework does not seem to be available. Make sure that either the submodule is checked out or that your MOOSE_DIR points to the correct location)
  endif
endif
$(info Using moose framework at $(MOOSE_DIR))

# framework
FRAMEWORK_DIR      := $(MOOSE_DIR)/framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
PHASE_FIELD      := yes
TENSOR_MECHANICS := yes
HEAT_CONDUCTION  := yes
include $(MOOSE_DIR)/modules/modules.mk
###############################################################################

# dep apps
APPLICATION_DIR    := $(CURDIR)
MAGPIE_DIR         := $(CURDIR)
APPLICATION_NAME   := magpie
BUILD_EXEC         := yes
DEP_APPS           := $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include            $(FRAMEWORK_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here
