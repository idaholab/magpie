###############################################################################
################### FFTW3 - fast fourier transform library ####################
###############################################################################

# check if FFTW3 is installed
ifeq ($(shell pkg-config fftw3 && echo go),go)
ADDITIONAL_INCLUDES += $(shell pkg-config fftw3 --cflags)
ADDITIONAL_CPPFLAGS += -DFFTW3_ENABLED
ADDITIONAL_LIBS += $(shell pkg-config fftw3 --libs)
endif
