SRC_DIR = ../NRs/ansi/recipes
OTHER_DIR = ../NRs/ansi/other
SRCS = main.c $(OTHER_DIR)/nrutil.c $(SRC_DIR)/gaussj.c

fit : $(SRCS) ../NRs/ansi/other/nr.h
	cc -o fit $(SRCS) -I $(OTHER_DIR)

clean: 
	rm -f fit