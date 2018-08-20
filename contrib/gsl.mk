###############################################################################
######################## GNU Scientific Library (GSL) #########################
###############################################################################

# check if GSL is installed
ifeq ($(shell pkg-config gsl && echo go),go)
ADDITIONAL_INCLUDES += $(shell pkg-config gsl --cflags)
ADDITIONAL_CPPFLAGS += -DGSL_ENABLED
ADDITIONAL_LIBS += $(shell pkg-config gsl --libs)
endif
