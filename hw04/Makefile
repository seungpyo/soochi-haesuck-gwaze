CC = cc
CFLAGS = -g -I ../NRs/ansi/other
SRCS = main.c \
	../NRS/ansi/recipes/ran1.c \
	../NRS/ansi/recipes/gasdev.c
OBJS = $(SRCS:.c=.o)
TARGET = runHw04

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(CFLAGS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
