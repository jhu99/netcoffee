OS = $(shell uname -s)
CXX = g++-4.6
DEBUG = yes

# Default mode is "Release"
DEFAULT_MODE  = Debug
MODE         ?= $(DEFAULT_MODE)

# If mode is something other than "Debug" or "Release",throw a fit
ifneq ($(MODE),Debug)
ifneq ($(MODE),Release)
$(error MODE must be one of {Debug,Release})
endif
endif

LEMON = include/lemon-1.2.3

ifeq ($(MODE),Debug)
	CXXFLAGS = -Wall -g3 -DDEBUG -std=c++0x -DVERBOSE -I$(LEMON)/ -Isrc/ -Isrc/input/ -Isrc/algorithms/
else
	CXXFLAGS = -Wall -O3 -ffast-math -fcaller-saves -finline-functions -std=c++0x -DNDEBUG -I$(LEMON)/ -Isrc/ -Isrc/input/ -Isrc/algorithms/
endif

all: netcoffee move

netcoffee: src/main.cpp src/verbose.o $(LEMON)/lemon/arg_parser.o
	${CXX} ${CXXFLAGS} -o $@ $^ 

move:
	mv netcoffee ./bin

#lemon: lemon-config lemon-make

#lemon-config:
#$(LEMON)/configure

#lemon-make:
#$(LEMON)/make	

clean:
	rm ./bin/netcoffee
