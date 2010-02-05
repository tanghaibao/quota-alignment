CC=gcc
CFLAGS=-O3 -Wall -ansi
PROG=bk_cliques

$(PROG): $(PROG).c
	$(CC) $(CFLAGS) $? -o $@

test:
	$(MAKE)
	$(PROG) <toy_matrix

