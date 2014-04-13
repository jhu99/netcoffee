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
MYSQLCPPCONN = include/mysqlcppconn-1.1.3
BOOST = include/boost_1_55_0

ifeq ($(MODE),Debug)
	CXXFLAGS = -Wall -g3 -DDEBUG -std=c++0x -fopenmp -DVERBOSE -I$(LEMON)/ -I$(MYSQLCPPCONN)/include/ -I$(MYSQLCPPCONN)/include/cppconn/ -I$(BOOST)/ -Isrc/ -Isrc/input/ -Isrc/algorithms/  -Isrc/format/ -I-L$(MYSQLCPPCONN)/lib/
else
	CXXFLAGS = -Wall -O3 -ffast-math -fcaller-saves -finline-functions -std=c++0x -fopenmp -DNDEBUG -I$(LEMON)/ -I$(MYSQLCPPCONN)/include/ -Isrc/ -Isrc/input/ -Isrc/algorithms/ -Isrc/format/ -L$(MYSQLCPPCONN)/lib/
endif

all: netcoffee move

netcoffee: src/main.cpp src/verbose.o
	${CXX} ${CXXFLAGS} -o $@ $^ -lmysqlcppconn

move:
	sudo mv netcoffee /usr/lib/cgi-bin/netcoffee.cgi

#lemon: lemon-config lemon-make

#lemon-config:
#$(LEMON)/configure

#lemon-make:
#$(LEMON)/make	

clean:
	rm ./bin/netcoffee.cgi $(CGI_INTERFACE)/*.o
