CXX=g++
CXXFLAGS=-O3 -pthread -lboost_system -lboost_program_options -lcurl
SRC_DIR=newGASS-base
SRCS=$(SRC_DIR)/Individuo.cpp $(SRC_DIR)/GA.cpp $(SRC_DIR)/newGASS.cpp gasscli.cpp
OBJS=$(SRCS:.cpp=.o)
TARGET=gasscli

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)