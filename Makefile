CXX = g++
CXXFLAGS  += -O2 -Wall -Werror -Wextra
RCXXFLAGS := `root-config --cflags --libs`
LDFLAGS   += -lconf -L./git/config/lib
RLDFLAGS  := -lEG -lTMVA -lTMVAGui

BUILDDIR = ./build

RPPSRCS = $(wildcard *.C)
RPPEXES = $(patsubst %.C,%,$(RPPSRCS))
RPPDEPS = $(patsubst %.C,$(BUILDDIR)/%.d,$(RPPSRCS))

CPPSRCS = $(wildcard *.cpp)
CPPEXES = $(patsubst %.cpp,%,$(CPPSRCS))
CPPDEPS = $(patsubst %.cpp,$(BUILDDIR)/%.d,$(CPPSRCS))

EXES = $(RPPEXES) $(CPPEXES)
DEPS = $(RPPDEPS) $(CPPDEPS)

.PHONY: all clean

all: $(RPPEXES) $(CPPEXES)

$(RPPEXES) : % : %.C
	@mkdir -p $(BUILDDIR)/$(@D)
	$(CXX) $(CXXFLAGS) $(RCXXFLAGS) -MMD -MF $(BUILDDIR)/$(@D)/$(*F).d $< -o $@ \
		$(LDFLAGS) $(RLDFLAGS)

$(CPPEXES) : % : %.cpp
	@mkdir -p $(BUILDDIR)/$(@D)
	$(CXX) $(CXXFLAGS) -MMD -MF $(BUILDDIR)/$(@D)/$(*F).d $< -o $@ \
		$(LDFLAGS)

clean:
	@$(RM) $(EXES) $(DEPS)
	@rm -rf $(BUILDDIR)/*

-include $(DEPS)
