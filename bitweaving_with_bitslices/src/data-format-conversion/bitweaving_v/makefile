  CXX = g++

  CXXFLAGS = -O3

  # the build target executable:
  TARGET = create_bitslices

  all: $(TARGET)

  $(TARGET): $(TARGET).cc
		$(CXX) $(CXXFLAGS) -o $(TARGET) $(TARGET).cc

  clean:
		$(RM) $(TARGET)