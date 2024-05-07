BIN_DIR := $(shell pwd)/bin
LIB_DIR := $(shell pwd)/lib
SRC_DIR := $(shell pwd)/src
INCLUDE_DIR := $(shell pwd)/include

OS := $(shell uname -s)

.PHONY: all clean

all: $(BIN_DIR)/tgsfilter

$(BIN_DIR)/tgsfilter: $(SRC_DIR)/TGSFilter.cpp | $(BIN_DIR)
ifeq ($(OS),Darwin)
	g++ --std=c++11 -g -O3 $< -lhts -lz -lisal -ldeflate -pthread -fPIE -I$(INCLUDE_DIR) -o $@
else
	g++ -L$(LIB_DIR) --std=c++11 -g -O3 $< -lhts -lz -lisal -ldeflate -pthread -fPIE -I$(INCLUDE_DIR) -o $@
endif

$(BIN_DIR):
	mkdir -p $(BIN_DIR)
clean:
	rm -rf $(BIN_DIR)
