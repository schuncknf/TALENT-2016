TARGET = libhffs.a

CC =gcc-5
CXX=g++-5
AR=ar

DEBUG_CFLAGS  =-g -fno-inline-functions -O0 -D_GLIBCXX_DEBUG
RELEASE_CFLAGS =-O3 -march=native -DARMA_NO_DEBUG
COMMON_CFLAGS =-fPIC -Wall -Wextra -std=gnu++11 -Icpp/include

#CFLAGS  =$(DEBUG_CFLAGS) $(COMMON_CFLAGS)
CFLAGS  =$(RELEASE_CFLAGS) $(COMMON_CFLAGS)
LDFLAGS = -larmadillo -lgsl -lgslcblas

CPP_FILES=cpp/src/Basis.cpp cpp/src/FieldMinnesota.cpp cpp/src/Fields.cpp cpp/src/FieldSph.cpp cpp/src/HartreeFockBogo.cpp cpp/src/HartreeFock.cpp \
            cpp/src/quadrature.cpp cpp/src/ReducedSpBasis.cpp cpp/src/Solver.cpp cpp/src/SpBasis.cpp cpp/src/System.cpp cpp/src/NeutronDrop.cpp \
            cpp/src/Nucleus.cpp cpp/src/Interaction.cpp cpp/src/RawInteraction.cpp cpp/src/MinnesotaS0.cpp
OBJ_FILES=$(CPP_FILES:.cpp=.o)

DEPS=$(CPP_FILES:.cpp=.d)

all: $(TARGET) pybind

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
	rm -rf $(TARGET) $(DEPS)


%.d: %.cpp
	$(CC) -c -MM $< -o $@ $(CFLAGS)

-include $(DEPS)
