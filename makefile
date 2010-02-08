CC=gcc
CFLAGS=-O3 -Wall -ansi -pedantic-errors
PROG=bk_cliques

$(PROG): $(PROG).c
	$(CC) $(CFLAGS) $? -o $@

install:
	install $(PROG) -t ~/bin

test:
	$(MAKE)
	$(PROG) <toydata

