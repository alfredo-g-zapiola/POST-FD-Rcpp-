PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)
export CMD_INSTALL_FLAGS = --preclean --no-multiarch --with-keep.source

############################
# C++ compiling and linking
###########################
mkCxxCompiler?=g++
CXX = ${mkCxxCompiler}
MAKEFILEH_DIR=${PACS_ROOT}
#
# include $(MAKEFILEH_DIR)/Makefile.inc
ifeq ($(DEBUG),yes)

else
 
endif
ifeq ($(DEBUG),no)
 # Since aa 18-19 I have decided to indicate the libraries even if DEBUG=no.
  # Remember however that normally release code does not have hardwired in the
  # executables the directories where the loader looks for libraries.
  # You take the responsibility of installing them in the right place.
  # But here I want to avoid students messing around with LD_LIBRARY_PATH
  # or ldconfig.
  LDFLAGS=-Wl,-rpath=. -Wl,-rpath=$(PACS_ROOT)/lib#
  OPTFLAGS=-O3 -funroll-loops #maybe you want -O3
  DEFINES+=-DNDEBUG -D MYNDEBUG
 else
   OPTFLAGS=-g
# If debugging we use the local dynamic libraries and avoid ldconfig -d
# or setting LD_LIBRARY_PATH
  LDFLAGS=-Wl,-rpath=. -Wl,-rpath=$(PACS_ROOT)/lib# 
endif
CXXFLAGS+=$(OPTFLAGS) # -fopenmp -lpthread

# Includes 



export INCLS += -I/usr/local/include -I./src  -I./inst/include -I/usr/include/R/ -I/usr/lib64/R/library/Rcpp/include -I/usr/lib64/R/library/RcppArmadillo/include -I$(PACS_ROOT)/include -I/usr/lib64/R/library/splines2/include/ 
# INCLS += -I/usr/lib/R/library/Rcpp/include/ -I/usr/lib/R/library/Rcpp/include/
INCLS += -I./inst/dirichlet-cpp/



export DEFINES
	
  
	

# linker flags
# linking libraries
export LDLIBS
# cppad
LDLIBS += -L/usr/local/lib/ -lcppad_lib -lcppad_ipopt # -L/usr/lib64/R/lib
# ipopt 
LDLIBS += -L/usr/lib -lipopt
# splines3
# LDLIBS += -L/usr/lib/R/library/splines2/libs/ -Wl,-rpath,/usr/lib/R/library/splines2/libs/
# Rcpp
LDLIBS += -L/usr/lib64/R/library/Rcpp/libs/ -Wl,-rpath,/usr/lib64/R/library/Rcpp/libs/
# R botth normal and 64 (to change?TODO)
# LDLIBS += -L/usr/lib64/R/lib/ -L/usr/lib/R/lib/  -Wl,-rpath,/usr/lib64/R/lib/ -Wl,-rpath,/usr/lib/R/lib/
#Rcpparma
LDLIBS += -L/usr/lib64/R/library/RcppArmadillo/libs/ -Wl,-rpath,/usr/lib64/R/library/RcppArmadillo/libs/ 
# splines 2
LDLIBS += -L/usr/lib64/R/library/splines2/libs/ -Wl,-rpath,/usr/lib64/R/library/splines2/libs/
# R
LDLIBS += -L/usr/lib64/R/lib -lR -L/usr/lib/R/lib/ -lR
# PACS
LDLIBS += -L${PACS_ROOT}/lib -Wl,-rpath,${PACS_ROOT}/lib -lpacs -lquadrature -lMesh1D  -lNewton -lMesh -lquadrules -lblas	


# openmp
#DEFINES += -fopenmp  -openmp


# Loader flags
export LDFLAGS = -Wl,-rpath,/usr/local/lib/ 
# -Wl,rpath,. # just need to look for cppad, since it is .so it loads ipopt by itself.


# Suggestion by Prof. Formaggia: since we are only uning c++ I set the linker LINK.o for object files
# so that it uses the c++ compiler. Make has no separate macro for the linker! 
# note this way the std lib is loaded automatically for eg.
LINK.o = $(CXX) $(LDFLAGS) $(LDLIBS)

# Warnings
export WARNFLAGS # Note these should not be changed, they are in /etc/R/Makeconf= -Wall -Wextra -Wsuggest-override -Wnon-virtual-dtor


export STDFLAGS= -std=c++17


export CPPFLAGS=$(INCLS) $(DEFINES)
export CXXFLAGS=$(OPTFLAGS) $(STDFLAGS) $(WARNFLAGS)

# Documentation
# Common file for Doxygen documentation
DOXYFILE=$(PACS_ROOT)/DoxyfileCommon


###################################
## tests
##################################
TEST_SRCS  = $(wildcard tests/*.cpp)
TEST_SRCS +=  $(filter-out Rcpp_fdpot.cpp RcppExports.cpp, $(wildcard src/*.cpp))
# get the corresponding object file
TEST_OBJS = $(TEST_SRCS:.cpp=.o)
TEST_EXEC=$(TEST_SRCS:.cpp=)

#PARALLEL_TEST_SRCS  = $(filter-out main_tests.cpp, $(wildcard tests/*.cpp))
#PARALLEL_TEST_SRCS +=  $(filter-out Rcpp_fdpot.cpp RcppExports.cpp, $(wildcard src/*.cpp))

#PARALLEL_TEST_OBJS = $(PARALLEL_TEST_SRCS:.cpp=.o)
#PARALLEL_TEST_EXEC=$(PARALLEL_TEST_SRCS:.cpp=)


all: check clean

.PHONY: all tests debugging
.DEFAULT_GOAL: install

install-cppad:
	# git clone https://github.com/coin-or/CppAD.git cppad.git
	
install:
	R -e "Rcpp::compileAttributes()"
	cd ..;\
	R CMD INSTALL FdPot $(CMD_INSTALL_FLAGS)
build:
	cd ..;\
	R CMD build --no-manual $(PKGSRC)
	
build-cran:
	$(shell "$(R_HOME).build/compile_attributes.R" )
	cd ..;\
	R CMD build $(PKGSRC) 
	
check: build-cran
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --as-cran

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/
	$(RM) $(TEST_EXEC) $(TEST_OBJS)
	-$(RM) make.dep

docs:
	cd ./src;\
	doxygen $(DOXYFILE)

tests: make.dep $(TEST_EXEC) 



debugging: tests
	gdb $(TEST_EXEC)
	
$(TEST_EXEC): $(TEST_OBJS)

$(TEST_OBJS): $(TEST_SRCS)


make.dep: $(TEST_SRCS);
	-\rm make.dep	
	for f in $(TEST_SRCS); do \
	$(CXX) $(STDFLAGS) $(CPPFLAGS) -MM $$f >> make.dep; \
	done
-include make.dep

	
