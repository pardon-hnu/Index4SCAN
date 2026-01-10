# =======================
# Compiler & Flags
# =======================
CXX      := g++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra -march=native
INCLUDES := -IHeadFile

# =======================
# Files
# =======================
BIN  := main
SRC  := main.cpp
OBJ  := main.o

# =======================
# Default Target
# =======================
all: $(BIN)

# =======================
# Link
# =======================
$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

# =======================
# Compile
# =======================
$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# =======================
# Clean
# =======================
clean:
	rm -f $(OBJ) $(BIN)

rebuild: clean all
