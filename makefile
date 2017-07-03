# Change these to suit your system
CXX      ?= g++
IFLAGS   = -I/usr/include/
LFLAGS   = -L/usr/lib
CFLAGS   = -Wall -Wno-reorder -Wno-sign-compare

# Flags for optimized compilation
CFLAGS_O = -O3 -fopenmp $(CFLAGS)
LIBS_O   = -lnetcdf -lgsl -lcblas

# Flags for debug compilation
CFLAGS_D = -fPIC -g -pg -rdynamic -fprofile-arcs -ftest-coverage -coverage -DDEBUG $(CFLAGS)
LIBS_D   = -lnetcdf -lgsl -lcblas -L build/gtest -lgtest -lpthread

# Don't change below here
SRCDIR   = src/
BUILDDIR = build
BUILDDIR_O = build/opt
BUILDDIR_D = build/debug
EXE_O    = gridpp
EXE_D    = gridpp_debug

PREFIX = $(DESTDIR)/usr
BINDIR = $(PREFIX)/bin

CORESRC 	= $(wildcard src/*.cpp)
CALSRC  	= $(wildcard src/Calibrator/*.cpp)
FILESRC 	= $(wildcard src/File/*.cpp)
DOWNSRC 	= $(wildcard src/Downscaler/*.cpp)
PARSRC 	= $(wildcard src/ParameterFile/*.cpp)
DRVSRC  	= src/Driver/Gridpp.cpp
DRVOBJ_O	= $(BUILDDIR_O)/Driver/Gridpp.o
DRVOBJ_D	= $(BUILDDIR_D)/Driver/Gridpp.o
KFSRC  	= src/Driver/Kf.cpp
KFOBJ_O	= $(BUILDDIR_O)/Driver/Kf.o
KFOBJ_D	= $(BUILDDIR_D)/Driver/Kf.o
TRAINRC  = src/Driver/Train.cpp
TRAINOBJ_O= $(BUILDDIR_O)/Driver/Train.o
TRAINOBJ_D= $(BUILDDIR_D)/Driver/Train.o
SRC     	= $(CORESRC) $(CALSRC) $(FILESRC) $(DOWNSRC) $(PARSRC)
HEADERS 	= $(SRC:.cpp=.h)
OBJ0_O   = $(patsubst src/%,$(BUILDDIR_O)/%,$(SRC))
OBJ0_D  	= $(patsubst src/%,$(BUILDDIR_D)/%,$(SRC))
OBJ_O   	= $(OBJ0_O:.cpp=.o)
OBJ_D   	= $(OBJ0_D:.cpp=.o)
TESTSRC 	= $(wildcard src/Testing/*.cpp)
TESTS0  	= $(patsubst src/Testing%,testing%,$(TESTSRC))
TESTS   	= $(TESTS0:.cpp=.exe)
INCS    	= makefile $(HEADERS)

.PHONY: tags count coverage doxygen test

default: $(EXE_O)

all: gridpp gridpp_train gridpp_kf

debug: $(EXE_D)

$(BUILDDIR):
	@mkdir build build/Calibrator build/Downscaler build/File build/ParameterFile build/Driver build/Testing

$(BUILDDIR_O)/%.o : src/%.cpp $(INCS)
	$(CXX) $(CFLAGS_O) $(IFLAGS) -c $< -o $@

$(BUILDDIR_D)/%.o : src/%.cpp $(INCS)
	$(CXX) $(CFLAGS_D) $(IFLAGS) -c $< -o $@

$(BUILDDIR_D)/%.E : src/%.cpp $(INCS)
	$(CXX) $(CFLAGS_D) $(IFLAGS) -c $< -o $@ -E

$(EXE_O): $(OBJ_O) $(DRVOBJ_O) makefile
	$(CXX) $(CFLAGS_O) $(LFLAGS) $(OBJ_O) $(DRVOBJ_O) $(LIBS_O) -o $@

gridpp_kf: $(OBJ_O) $(KFOBJ_O) makefile
	$(CXX) $(CFLAGS_O) $(LFLAGS) $(OBJ_O) $(KFOBJ_O) $(LIBS_O) -o $@

gridpp_train: $(OBJ_O) $(TRAINOBJ_O) makefile
	$(CXX) $(CFLAGS_O) $(LFLAGS) $(OBJ_O) $(TRAINOBJ_O) $(LIBS_O) -o $@

$(EXE_D): $(OBJ_D) $(DRVOBJ_D) makefile gtest
	$(CXX) $(CFLAGS_D) $(LFLAGS) $(OBJ_D) $(DRVOBJ_D) $(LIBS_D) -o $@

gridpp_kf_debug: $(OBJ_D) $(KFOBJ_D) makefile gtest
	$(CXX) $(CFLAGS_D) $(LFLAGS) $(OBJ_D) $(KFOBJ_D) $(LIBS_D) -o $@

gridpp_train_debug: $(OBJ_D) $(TRAINOBJ_D) makefile gtest
	$(CXX) $(CFLAGS_D) $(LFLAGS) $(OBJ_D) $(TRAINOBJ_D) $(LIBS_D) -o $@

test: gtest $(TESTS)
	./runAllTests.sh

$(BUILDDIR_D)/libgridpp.a: $(OBJ_D)
	ar rvs $@ $(OBJ_D)

$(BUILDDIR_D)/libgridpp.so: $(OBJ_D)
	$(CXX) $(CFLAGS_D) -shared -Wl,-soname,libgridpp.so -o $@ $(OBJ_D)

testing/%.exe: $(BUILDDIR_D)/Testing/%.o $(INCS) $(BUILDDIR_D)/libgridpp.so gtest
	$(CXX) $(CFLAGS_D) $< $(LFLAGS) -L$(BUILDDIR_D) -Wl,-rpath,./$(BUILDDIR_D) -lgridpp $(LIBS_D) -o $@

#testing/%.exe: $(BUILDDIR_D)/Testing/%.o $(INCS) $(OBJ_D) gtest
#	$(CXX) $(CFLAGS_D) $(OBJ_D) $< $(LFLAGS) $(LIBS_D) -o $@

count:
	@wc src/*.h src/*.cpp src/*/*.h src/*/*.cpp -l | tail -1

clean:
	rm -rf build/*/*.o build/*/*/*.o build/*/*.E build/*/*/*.E gmon.out $(EXE) testing/*.exe\
		*.gcno build/*/*.gcda build/*/*.gcno build/*/*/*.gcda build/*/*/*.gcno\
		coverage/* coverage.* build/gtest gridpp gridpp_kf gridpp_debug gridpp_kf_debug build/*/libgridpp.a build/*/libgridpp.so

tags:
	ctags -R --c++-kinds=+pl --fields=+iaS --extra=+q -f tags ./*.h ./*.cpp */*.h */*.cpp

install:
	install -D gridpp $(BINDIR)/gridpp

doxygen:
	doxygen doxygen/config

coverage: test
	#rm -f build/*.gcno build/*.gcda build/*/*.gcno build/*/*.gcda
	lcov -b `pwd` -c -i -d `pwd`/ -o coverage.init
	./runAllTests.sh
	lcov -b `pwd` -c -d `pwd`/ -o coverage.run
	lcov -a coverage.init -a coverage.run -o coverage.total
	lcov -e coverage.total "`pwd`/*" -o coverage.total.filtered
	lcov -r coverage.total.filtered "`pwd`/*Testing*" -o coverage.info
	genhtml -o ./coverage/ coverage.info

depend:
	makedepend -- $(CFLAGS) -- $(SRC)

gtest: build/gtest/libgtest.a

build/gtest/libgtest.a: /usr/src/gtest
	mkdir -p build/gtest
	cd build/gtest && cmake /usr/src/gtest && cmake --build .
