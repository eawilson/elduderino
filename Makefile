CC          = gcc
CFLAGS      = -Wall -O2
LDFLAGS     = -lz -lm
prefix      = /usr/local
exec_prefix = $(prefix)/bin

src = $(wildcard *.c)
obj = $(src:.c=.o)

elduderino: $(obj)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) elduderino

.PHONY: install
install:
	cp elduderino $(exec_prefix)







