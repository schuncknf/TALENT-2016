TARGET = hffs.a

CC =gcc
CXX=g++
AR=ar

DEBUG_CFLAGS  =-g -fno-inline-functions -O0 -D_GLIBCXX_DEBUG
RELEASE_FLAGS =-O3 -march=native -Wall -fPIC -Wextra -DARMA_NO_DEBUG
COMMON_CFLAGS =-Wall -Wextra -std=c++11 -Iinclude -fPIC

CFLAGS  =$(DEBUG_CFLAGS) $(COMMON_CFLAGS)
#CFLAGS  =$(RELEASE_CFLAGS) $(COMMON_CFLAGS)
LDFLAGS = -larmadillo

CPP_FILES=src/Basis.cpp src/FieldMinnesota.cpp src/Fields.cpp src/FieldSph.cpp src/HartreeFockBogo.cpp src/HartreeFock.cpp \
           src/quadrature.cpp src/ReducedSpBasis.cpp src/Solver.cpp src/SpBasis.cpp src/System.cpp
OBJ_FILES=$(CPP_FILES:.cpp=.o)

DEPS=$(CPP_FILES:.cpp=.d)

all: $(TARGET)

$(TARGET): $(OBJ_FILES)
	$(AR) rcs $@ $(OBJ_FILES)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CFLAGS)

.PHONY: clean mrproper depend

depend: $(DEPS)

clean:
	rm -rf $(OBJ_FILES) $(DEPS)

mrproper: clean
	rm -rf $(TARGET)


%.d: %.cpp
	$(CC) -c -MM $< -o $@ $(CFLAGS)

-include $(DEPS)
