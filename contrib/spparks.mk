###############################################################################
############################### SPPARKS #######################################
###############################################################################

# run Make.sh to recreate the stytle_*.h files
DUMMY := $(shell $(APPLICATION_DIR)/contrib/src/Make.sh style)

# source files
SPPARKS_cppsrcfiles := $(shell find $(ELK_DIR)/contrib/SPPARKS/src -maxdepth 1 -name "*.cpp" | grep -v "main.cpp")
SPPARKS_csrcfiles   := $(shell find $(ELK_DIR)/contrib/SPPARKS -name "*.c")

# object files
MAGPIE_objects     += $(patsubst %.cpp, %.$(obj-suffix), $(SPPARKS_cppsrcfiles))
MAGPIE_objects     += $(patsubst %.c, %.$(obj-suffix), $(SPPARKS_csrcfiles))

# dependencies (C, C++ files only)
MAGPIE_deps += $(patsubst %.cpp, %.$(obj-suffix).d, $(SPPARKS_cppsrcfiles)) \
               $(patsubst %.c, %.$(obj-suffix).d, $(SPPARKS_csrcfiles))

# header files
include_dirs   += $(APPLICATION_DIR)/contrib
app_INCLUDE    := $(foreach i, $(include_dirs), -I$(i)) $(ADDITIONAL_INCLUDES)
app_INCLUDES   += $(app_INCLUDE)
