TARGET = libhffs.a

CC =gcc-5
CXX=g++-5
AR=ar

DEBUG_CFLAGS  =-g -fno-inline-functions -O0 -D_GLIBCXX_DEBUG
RELEASE_CFLAGS =-O3 -march=native -DARMA_NO_DEBUG
COMMON_CFLAGS =-fPIC -Wall -Wextra -std=gnu++11 -Iinclude

#CFLAGS  =$(DEBUG_CFLAGS) $(COMMON_CFLAGS)
CFLAGS  =$(RELEASE_CFLAGS) $(COMMON_CFLAGS)
LDFLAGS = -larmadillo -lgsl -lgslcblas

CPP_FILES=src/Basis.cpp src/FieldMinnesota.cpp src/Fields.cpp src/FieldSph.cpp src/HartreeFockBogo.cpp src/HartreeFock.cpp \
            src/quadrature.cpp src/ReducedSpBasis.cpp src/Solver.cpp src/SpBasis.cpp src/FullSpBasis.cpp src/System.cpp src/NeutronDrop.cpp \
            src/Nucleus.cpp src/Interaction.cpp src/RawInteraction.cpp src/MinnesotaS0.cpp
OBJ_FILES=$(CPP_FILES:.cpp=.o)

DEPS=$(CPP_FILES:.cpp=.d)

all: $(TARGET) pybind

pybind: $(TARGET)
	cd ../python/bindings;./setup.py install --user;cd -

test: src/test.o $(TARGET)
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
	rm -rf $(TARGET) $(DEPS)


%.d: %.cpp
	$(CC) -c -MM $< -o $@ $(CFLAGS)

-include $(DEPS)
