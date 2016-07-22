TARGET1 = ex_quad

CC =gcc-5
CXX=g++-5
DEBUG_CFLAGS  =-g -fno-inline-functions -O0 -D_GLIBCXX_DEBUG
RELEASE_FLAGS =-O3 -march=native -Wall -fPIC -Wextra -DARMA_NO_DEBUG
COMMON_CFLAGS =-Wall -Wextra -std=c++11 -I../include

CFLAGS  =$(DEBUG_CFLAGS) $(COMMON_CFLAGS)
#CFLAGS  =$(RELEASE_CFLAGS) $(COMMON_CFLAGS)
LDFLAGS = -larmadillo

CPP_FILES1=ex_quad.cpp ../src/quadrature.cpp
OBJ_FILES1=$(CPP_FILES1:.cpp=.o)

DEPS1=$(CPP_FILES1:.cpp=.d)

all: $(TARGET1)

$(TARGET1): $(OBJ_FILES1)
	$(CXX) $^ -o $@ $(LDFLAGS) $(CFLAGS)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CFLAGS) -larmadillo

clean:
	rm -rf $(OBJ_FILES1) $(DEPS1)

mrproper: clean
	rm -rf $(TARGET1)


%.d: %.cpp
	$(CC) -c -MM $< -o $@ $(CFLAGS)

-include $(DEPS1)
