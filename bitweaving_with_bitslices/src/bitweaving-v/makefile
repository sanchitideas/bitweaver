  CXX = g++

  CXXFLAGS  = -mavx512f -O3

  # the build target executable:
  TARGET = bitweaver

  all: $(TARGET)

  $(TARGET): bitweaving-v.cc bitweaving_v_test.cc
		$(CXX) $(CXXFLAGS) -o $(TARGET) bitweaving-v.cc bitweaving_v_test.cc

  clean:
		$(RM) $(TARGET)