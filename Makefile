
TOOL = 
GPP = $(TOOL)g++
BIN = gann
SRC = gann.cpp
CFLAGS =
LDFLAGS = -lm

all:
	@$(GPP) $(CFLAGS) $(SRC) -o $(BIN) $(LDFLAGS)
	@echo ... BUILD DONE ...
	
clean:
	@rm -f $(BIN)
	@echo ... CLEAN DONE ...