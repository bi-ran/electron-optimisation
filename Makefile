CXX = clang++
CXXFLAGS += -O2 -Wall -Werror -Wextra
ROOTFLAGS := `root-config --cflags --libs` -lEG

BUILDDIR = ./build

CSRCS   = $(wildcard *.C)
CEXES   = $(patsubst %.C,%,$(CSRCS))
CDEPS   = $(patsubst %.C,$(BUILDDIR)/%.d,$(CSRCS))

CPPSRCS = $(wildcard *.cpp)
CPPEXES = $(patsubst %.cpp,%,$(CPPSRCS))
CPPDEPS = $(patsubst %.cpp,$(BUILDDIR)/%.d,$(CPPSRCS))

EXES = $(CEXES) $(CPPEXES)
DEPS = $(CDEPS) $(CPPDEPS)

.PHONY: all clean

all: $(CEXES) $(CPPEXES)

$(CEXES) : % : %.C
	@mkdir -p $(BUILDDIR)/$(@D)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -MMD -MF $(BUILDDIR)/$(@D)/$(*F).d $< -o $@

$(CPPEXES) : % : %.cpp
	@mkdir -p $(BUILDDIR)/$(@D)
	$(CXX) $(CXXFLAGS) -MMD -MF $(BUILDDIR)/$(@D)/$(*F).d $< -o $@

clean:
	@$(RM) $(EXES) $(DEPS)
	@rm -rf $(BUILDDIR)/*

-include $(DEPS)
