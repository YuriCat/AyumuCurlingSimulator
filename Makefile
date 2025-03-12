CXX        = c++
CXXFLAGS   = -std=c++20 -MMD -MP -fPIC
OPT        = -O3 -march=native -DNDEBUG
#OPT       := -O0 -g
LDFLAGS    =
LIBS       =
INCLUDES   =
SRC_DIR    = ./
BLD_DIR    = .
OBJ_DIR    = ./obj
SRCS       = $(wildcard $(SRC_DIR)/*.cpp)
TARGETOBJS = $(OBJ_DIR)/sample.o $(OBJ_DIR)/py_fast_simulator.o
OBJS       = $(filter-out $(TARGETOBJS), $(subst $(SRC_DIR),$(OBJ_DIR), $(SRCS:.cpp=.o)))
TARGET     = $(BLD_DIR)/sample
PYTARGET   = $(BLD_DIR)/fast_simulator.so
PYFLAGS    = -fPIC
PYLDFLAGS  = -shared
ifeq ($(shell uname),Darwin)
PYLDFLAGS += -undefined dynamic_lookup
endif
PYINCLUDES = $(INCLUDES) $(shell python3-config --includes) -I./pybind11/include/

DEPENDS    = $(OBJS:.o=.d) $(TARGETOBJS:.o=.d)

default: $(TARGET)
py: $(PYTARGET)
all: $(TARGET) $(PYTARGET)

$(TARGET): $(OBJS) $(OBJ_DIR)/sample.o $(LIBS)
	$(CXX) $(OPT) -o $@ $(OBJ_DIR)/sample.o $(OBJS) $(LIBS) $(LDFLAGS)

$(PYTARGET): $(OBJS) $(OBJ_DIR)/py_fast_simulator.o $(LIBS)
	$(CXX) $(OPT) -o $@ $(OBJ_DIR)/py_fast_simulator.o $(OBJS) $(LIBS) $(LDFLAGS) $(PYLDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@if [ ! -d $(OBJ_DIR) ]; \
		then echo "mkdir -p $(OBJ_DIR)"; mkdir -p $(OBJ_DIR); \
		fi
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDES) -o $@ -c $<

$(OBJ_DIR)/sample.o: $(SRC_DIR)/sample.cpp
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDES) -o $@ -c $<

$(OBJ_DIR)/py_fast_simulator.o: $(SRC_DIR)/py_fast_simulator.cpp
	$(CXX) $(CXXFLAGS) $(OPT) $(PYFLAGS) $(PYINCLUDES) -o $@ -c $<

clean:
	$(RM) -r $(OBJ_DIR) $(TARGET) $(PYTARGET)

-include $(DEPENDS)

.PHONY: all clean
