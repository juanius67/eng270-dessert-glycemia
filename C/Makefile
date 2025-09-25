CC ?= gcc
CFLAGS ?= -O3 -std=c11
PICFLAG = -fPIC

.PHONY: all clean libmodel

all: libmodel model.dll

libmodel: libmodel.so

libmodel.so: model.o
	$(CC) $(CFLAGS) -shared -o $@ $^

model.o: model.c model.h
	$(CC) $(CFLAGS) $(PICFLAG) -c -o $@ $<

model.dll: model.c model.h
	$(CC) $(CFLAGS) -shared -o $@ $<

clean:
	rm -f libmodel.so model.dll model.o
