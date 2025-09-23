CC=gcc
CFLAGS=-O2 -std=c11

all: libmodel.so model.dll

libmodel.so: model.c model.h
	$(CC) $(CFLAGS) -shared -fPIC -o $@ model.c

model.dll: model.c model.h
	$(CC) $(CFLAGS) -shared -o $@ model.c

clean:
	rm -f libmodel.so model.dll
