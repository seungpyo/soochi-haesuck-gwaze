CC = cc
CFLAGS = -g -I ../NRs/ansi/other
SRCS = main.c \
	../NRS/ansi/recipes/gaussj.c \
	../NRS/ansi/recipes/ludcmp.c \
	../NRs/ansi/recipes/svdcmp.c \
	../NRs/ansi/recipes/mprove.c \
	../NRs/ansi/recipes/nrutil.c \
	../NRs/ansi/recipes/lubksb.c \
	../NRs/ansi/recipes/pythag.c
OBJS = $(SRCS:.c=.o)
TARGET = runHw03

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(CFLAGS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
