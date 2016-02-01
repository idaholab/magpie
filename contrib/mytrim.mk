###############################################################################
################################ MYTRIM #######################################
###############################################################################

# source files
mytrim_srcfiles := $(shell find $(MYTRIM_DIR) -maxdepth 1 -name "*.C")

app_HEADERS := $(shell find $(MYTRIM_DIR) -maxdepth 1 -name "*.h")

# object files
ADDITIONAL_APP_OBJECTS     += $(patsubst %.C, %.$(obj-suffix), $(mytrim_srcfiles))

# dependencies (C, C++ files only)
ADDITIONAL_APP_DEPS += $(patsubst %.C, %.$(obj-suffix).d, $(mytrim_srcfiles))

# set MYTRIM_DATA_DIR to point to the crossection data
ADDITIONAL_CPPFLAGS += -DMYTRIM_DATA_DIR=\""$(MYTRIM_DIR)/data"\"
