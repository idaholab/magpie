###############################################################################
######################## GNU Scientific Library (GSL) #########################
###############################################################################
GSL_DIR ?= $(realpath $(APPLICATION_DIR)/contrib/gsl)

# check if gsl submodule is checked out
ifeq ($(wildcard $(GSL_DIR)/Makefile.am),)
  $(error Please checkout the gsl submodule or point GSL_DIR to a valid gsl checkout)
endif

#
# Set gsl as a dependency for magpie
#

$(APPLICATION_DIR)/lib/libmagpie-$(METHOD).la: $(GSL_DIR)/libgsl.la

#
# compile and build gsl first! (if no command line target was specified)
#

ifeq ($(MAKECMDGOALS),)

# configure GSL
ifeq ($(shell [ ! -s $(GSL_DIR)/Makefile -o $(GSL_DIR)/configure -nt $(GSL_DIR)/Makefile ] && echo go),go)
$(info Configuring GSL...)
$(info $(shell cd $(GSL_DIR) && ./configure))
endif

# make GSL
ifeq ($(shell [ ! -s $(GSL_DIR)/libgsl.la -o $(GSL_DIR)/Makefile -nt $(GSL_DIR)/libgsl.la ] && echo go),go)
$(info Building GSL...)
$(info $(shell $(MAKE) -C $(GSL_DIR)))
endif
endif

ADDITIONAL_INCLUDES += -I$(GSL_DIR)
ADDITIONAL_CPPFLAGS += -DGSL_ENABLED
ADDITIONAL_LIBS += -L$(GSL_DIR) -lgsl -L$(GSL_DIR)/cblas -lgslcblas -I$(GSL_DIR)
ADDITIONAL_APP_DEPS += $(APPLICATION_DIR)/contrib/gsl.d
