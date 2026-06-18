CXX := clang++
CXXFLAGS := -O2 -std=c++17

# Folder that CONTAINS the directory "itensor/"
ITENSOR_PARENT := ../..
CXXFLAGS += -I$(ITENSOR_PARENT) -I.
CXXFLAGS += -Wno-deprecated-declarations
CXXFLAGS += -Wno-deprecated-literal-operator

# If you built ITensor and it produced lib/libitensor.a under the itensor folder:
LDFLAGS := -L$(ITENSOR_PARENT)/lib
LDLIBS  := -litensor -framework Accelerate

SRC := main.cpp functions.cpp data.cpp TestModel.cpp StructFile.cpp TwoQubitTestModel.cpp BuildingHamiltonians.cpp OneSigmaSite.cpp BaseModel.cpp BaseModel2.cpp TwoSigmaSite.cpp ChainModel.cpp ChainModelMPS.cpp functionsMPS.cpp SigmaChainAnalysis.cpp GetData.cpp CorrelationScan.cpp
OBJ := $(SRC:.cpp=.o)

all: dbqite_itensor

dbqite_itensor: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp functions.h data.h TestModel.h StructFile.h TwoQubitTestModel.h BuildingHamiltonians.h OneSigmaSite.h BaseModel.h BaseModel2.h TwoSigmaSite.h ChainModel.h ChainModelMPS.h functionsMPS.h SigmaChainAnalysis.h GetData.h CorrelationScan.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) dbqite_itensor
