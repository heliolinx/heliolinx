#
# Makefile for heliolinx
#
# To build all programs (in parallel), run:      make -j4
# To install to ../bin, run:                     sudo make install
# To install to (say) /usr/local, run:           PREFIX=/usr/local make install
#

PREFIX ?= ..
CXXFLAGS ?= -O3

###

PROGRAMS_PATH = $(PREFIX)/bin
PROGRAMS = make_tracklets heliolinc heliolinc_omp heliolinc_lowmem heliovane link_purify link_planarity parse_clust2det_MPC80 parse_clust2det modsplit_hlfile merge_tracklet_files make_trailed_tracklets parse_trk2det

LIB = libheliolinx.a
LIB_SOURCES = solarsyst_dyn_geo01.cpp

SOURCES = $(LIB_SOURCES) $(PROGRAMS:%=%.cpp)

.PHONY: all
all: $(PROGRAMS)

.PHONY: install
install: $(PROGRAMS)
	mkdir -p $(PREFIX)/bin
	cp -p $(PROGRAMS) $(PREFIX)/bin

$(PROGRAMS): %: %.o $(LIB)
	$(CXX) $< -L. -lheliolinx -fopenmp $(OUTPUT_OPTION)

$(LIB): $(LIB_SOURCES:%.cpp=%.o)
	$(AR) -rc $@ $^

clean:
	rm -f $(PROGRAMS) $(LIB) $(SOURCES:%.cpp=%.o)
	rm -rf .deps

###

DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

%.o : %.cpp
%.o : %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	$(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -std=c++11 -fopenmp -I../include -c $(OUTPUT_OPTION) $<
	@ $(POSTCOMPILE)

$(DEPDIR): ; @mkdir -p $@

DEPFILES := $(SOURCES:%.cpp=$(DEPDIR)/%.d)
$(DEPFILES):

include $(wildcard $(DEPFILES))
