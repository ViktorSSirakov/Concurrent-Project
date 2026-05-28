CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17

TARGET = main
SOURCES = main.cpp datapoints.cpp
OBJECTS = $(SOURCES:.cpp=.o)

prepare: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

clear: clean