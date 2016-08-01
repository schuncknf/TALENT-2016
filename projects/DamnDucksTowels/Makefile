TARGET = libhffs.a

CC =gcc-5
CXX=g++-5
AR=ar

DEBUG_CFLAGS  =-g -fno-inline-functions -O0 -D_GLIBCXX_DEBUG
RELEASE_CFLAGS =-O3 -march=native -DARMA_NO_DEBUG
COMMON_CFLAGS =-fPIC -Wall -Wextra -std=gnu++11 -Icpp/include

CFLAGS  =$(DEBUG_CFLAGS) $(COMMON_CFLAGS)
#CFLAGS  =$(RELEASE_CFLAGS) $(COMMON_CFLAGS)
LDFLAGS = -larmadillo -lgsl -lgslcblas

CPP_FILES=cpp/src/utils.cpp cpp/src/Basis.cpp cpp/src/HartreeFockBogo.cpp cpp/src/HartreeFock.cpp \
            cpp/src/quadrature.cpp cpp/src/ReducedSpBasis.cpp cpp/src/Solver.cpp cpp/src/SpBasis.cpp cpp/src/System.cpp cpp/src/NeutronDrop.cpp \
            cpp/src/Nucleus.cpp cpp/src/Interaction.cpp cpp/src/RawInteraction.cpp cpp/src/MinnesotaS0.cpp cpp/src/MinnesotaRaw.cpp cpp/src/FullSpBasis.cpp
OBJ_FILES=$(CPP_FILES:.cpp=.o)

DEPS=$(CPP_FILES:.cpp=.d)

all: $(TARGET) pybind

indent:
	astyle -A1 -s2 -S -f -p -H -xe -k3 cpp/include/*.h cpp/src/*.cpp

pybind: $(TARGET)
	cd ./python/bindings;./setup.py install --user;cd -

test: cpp/src/test.o $(TARGET)
	$(CXX) $^ -o $@ $(CFLAGS) $(LDFLAGS)

$(TARGET): $(OBJ_FILES)
	$(AR) rcs $@ $(OBJ_FILES)

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

.PHONY: clean mrproper depend

depend: $(DEPS)

clean:
	rm -rf $(OBJ_FILES)

mrproper: clean
	cd ./python/bindings;./setup.py clean;rm -Rf build;rm -f hffs.pyc hffs.py hffs_wrap.cpp hffs_wrap.cpp.bkp;cd -
	rm -rf $(TARGET) $(DEPS)


%.d: %.cpp
	$(CC) -c -MM $< -o $@ $(CFLAGS)

-include $(DEPS)
