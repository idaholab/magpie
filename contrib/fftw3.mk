###############################################################################
################### FFTW3 - fast fourier transform library ####################
###############################################################################

# check if FFTW3 is installed
ifeq ($(shell pkg-config fftw3 && echo go),go)
ADDITIONAL_INCLUDES += $(shell pkg-config fftw3 --cflags)
ADDITIONAL_CPPFLAGS += -DFFTW3_ENABLED
FFTW3LIBS := $(shell pkg-config fftw3 --libs)
ADDITIONAL_LIBS += $(FFTW3LIBS)
LDFLAGS += $(FFTW3LIBS)
endif
