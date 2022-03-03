ROOT_DIR=$(shell pwd)
ODIR  = $(ROOT_DIR)/obj
SDIR  = $(ROOT_DIR)/src

CXX   = mpicxx
CFLAG = -std=c++11 -larmadillo -lopenblas -lmpi
 
DEPS  = $(shell ls $(SDIR)/*.h)
SRC   = $(shell ls $(SDIR)/*.cpp)
OBJ   = $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(SRC))

fssh.x : $(OBJ)
	$(CXX) -o $@ $^ $(CFLAG)

$(ODIR)/%.o : $(SDIR)/%.cpp $(DEPS) | $(ODIR)/.
	$(CXX) -c -o $@ $< $(CFLAG)

%/. : 
	mkdir -p $(patsubst %/.,%,$@)
	
.PRECIOUS: %/.
.PHONY: clean clean_dat clean_all
clean:
	rm -rf *.x $(ODIR)
clean_dat:
	rm -rf *.dat *.out
clean_all:
	rm -rf *.x *.dat *.out $(ODIR)
