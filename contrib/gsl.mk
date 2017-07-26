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
# compile and build gsl on demand
#

$(GSL_DIR)/libgsl.la: $(GSL_DIR)/Makefile
	@echo ===========================
	@echo ====== Building GSL =======
	@echo ===========================
	$(MAKE) -C $(shell dirname $@)

$(GSL_DIR)/Makefile: $(GSL_DIR)/configure
	@echo ===========================
	@echo ===== Configuring GSL =====
	@echo ===========================
	cd $(shell dirname $@) && ./configure

ADDITIONAL_INCLUDES += -I$(GSL_DIR)
ADDITIONAL_CPPFLAGS += -DGSL_ENABLED
ADDITIONAL_LIBS += -L$(GSL_DIR) -lgsl -L$(GSL_DIR)/cblas -lgslcblas -I$(GSL_DIR)
ADDITIONAL_APP_DEPS += $(APPLICATION_DIR)/contrib/gsl.d
