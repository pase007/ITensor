CXX := clang++
CXXFLAGS := -O2 -std=c++17

# Folder that CONTAINS the directory "itensor/"
ITENSOR_PARENT := ../..
CXXFLAGS += -I$(ITENSOR_PARENT) -I.
CXXFLAGS += -Wno-deprecated-declarations

# If you built ITensor and it produced lib/libitensor.a under the itensor folder:
LDFLAGS := -L$(ITENSOR_PARENT)/lib
LDLIBS  := -litensor -framework Accelerate

SRC := main.cpp functions.cpp data.cpp TestModel.cpp StructFile.cpp TwoQubitTestModel.cpp BuildingHamiltonians.cpp
OBJ := $(SRC:.cpp=.o)

all: dbqite_itensor

dbqite_itensor: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp functions.h data.h TestModel.h StructFile.h TwoQubitTestModel.h BuildingHamiltonians.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) dbqite_itensor