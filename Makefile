CXX = clang++
CXXFLAGS += -O2 -Wall -Werror -Wextra --std=c++14
ROOTFLAGS := `root-config --cflags --libs` -lEG

BUILDDIR = ./build

SRCS  = $(wildcard *.C)
EXES  = $(patsubst %.C,%,$(SRCS))
DEPS  = $(patsubst %.C,$(BUILDDIR)/%.d,$(SRCS))

.PHONY: all clean

all: $(EXES)

%: %.C
	@mkdir -p $(BUILDDIR)/$(@D)
	$(CXX) $(ROOTFLAGS) $(CXXFLAGS) -MMD -MF $(BUILDDIR)/$(@D)/$(*F).d $< -o $@

clean:
	@$(RM) $(EXES) $(DEPS)
	@rm -rf $(BUILDDIR)/*

-include $(DEPS)
