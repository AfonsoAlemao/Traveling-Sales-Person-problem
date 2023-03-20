TARGET = tsp
CC = gcc
CFLAGS = -Wall -g -fopenmp

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -fopenmp -Wall -o $@

.PHONY: default all clean

clean:
	-rm -f *.o
	-rm -f $(TARGET)
	-rm -f vgcore*