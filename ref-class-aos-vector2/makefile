CXX				= g++
CFLAGS 			+= -std=c++11

ifeq ($(config),debug)
	CFLAGS 		+= -g -D__DEBUG -O0# -shared-intel
	LDFLAGS 	+= -g 
else
	CFLAGS 		+= -O3#  -xHost -ipo
endif


####################################################
#add subdirs
DIRS 		= $(shell find * -type d)
#add include directives for subdirs
CFLAGS 		+= $(DIRS:%=-I%)


#collect all cpp files
SRCS := $(shell find . -name \*.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(SRCS:.cpp=.d)


###############################################################
EXEC		= liz
default: $(EXEC)

help:
	@echo 'Usage                 : make [clean] [openmp={gnu|intel}] [config=debug]'
	@echo 'Environment variables : CXX (c++ compiler) HDF5_ROOT_DIR (HDF5 dir with /include and /lib)'

release:
	make

debug:
	make config=debug

clean:
	rm -f $(OBJS) $(DEPS) $(PYHEADERS) $(EXEC) $(PYCS)

distclean: clean
	rm -f $(EXEC)

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif

# this generates a .h file containing a char[] with the python script in binary then
#you can just include this file to get the contents
%.pyh: %.py
	@ echo "Creating binary char for $< : $@"
	@ xxd -i $< > $@

%.d: %.cpp $(PYHEADERS)
	@ echo "Checking dependencies for $<"
# create and modify dependecy file .d to take into account the location subdir
	@ $(CXX) $(CFLAGS) -MM $< 2>/dev/null | sed -e "s@\(^.*\)\.o:@$(shell  dirname $<)/\1.d $(shell  dirname $<)/\1.o:@" > $@

%.o : %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

$(EXEC): $(PYHEADERS) $(OBJS)
	$(CXX) $(CFLAGS) -o $(EXEC) $(OBJS) $(LDFLAGS)

# DO NOT DELETE
