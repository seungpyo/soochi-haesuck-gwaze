CC = cc
CFLAGS = -g -I ../NRs/ansi/other
SRCS = main.c \
	../NRS/ansi/recipes/zbrak.c \
	../NRS/ansi/recipes/rtbis.c \
	../NRS/ansi/recipes/rtflsp.c \
	../NRS/ansi/recipes/rtsec.c \
	../NRS/ansi/recipes/rtnewt.c \
	../NRS/ansi/recipes/rtsafe.c \
	../NRs/ansi/recipes/bessj0.c \
	../NRs/ansi/recipes/bessj1.c \
	../NRs/ansi/other/nrutil.c
OBJS = $(SRCS:.c=.o)
TARGET = runHw02

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(CFLAGS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
