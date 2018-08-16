###############################################################################
######################## GNU Scientific Library (GSL) #########################
###############################################################################
GSL_DIR ?= $(realpath $(MAGPIE_DIR)/contrib/gsl)

# check if gsl submodule is checked out
ifeq ($(wildcard $(GSL_DIR)/Makefile.am),)
  $(error Please checkout the gsl submodule or point GSL_DIR to a valid gsl checkout)
endif

#
# Symlink GSL includes
#

ADDITIONAL_INCLUDES += -I$(MAGPIE_DIR)/contrib/gsl
ADDITIONAL_CPPFLAGS += -DGSL_ENABLED
contrib/gsl/gsl/.linked : $(MAGPIE_DIR)/.git/modules/contrib/gsl/HEAD
	@ln -sf `find $(GSL_DIR) -name gsl_\*.h -not -path $(GSL_DIR)/gsl` $(GSL_DIR)/gsl
	@touch $(MAGPIE_DIR)/contrib/gsl/gsl/.linked
all:: contrib/gsl/gsl/.linked

# configure GSL
ifeq ($(shell [ ! -s $(GSL_DIR)/Makefile -o $(GSL_DIR)/configure -nt $(GSL_DIR)/Makefile ] && echo go),go)
$(info Configuring GSL...)
$(info $(shell cd $(GSL_DIR) && ./configure))
endif

# make gsl modules
define GSL_MODULE_RULE
ADDITIONAL_APP_OBJECTS += $(GSL_DIR)/$(strip $(1))/libgsl$(strip $(1)).la
$(GSL_DIR)/$(strip $(1))/libgsl$(strip $(1)).la: $(shell find $(GSL_DIR)/$(strip $(1)) -name "*.c" -or -name "*.h")
	$(MAKE) -C $(GSL_DIR)/$(strip $(1))
endef
$(foreach module,sys specfunc integration err eigen complex vector, $(eval $(call GSL_MODULE_RULE, $(module))))
