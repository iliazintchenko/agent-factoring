CC = gcc
CXX = g++
CFLAGS = -O3 -march=native -Wall -Wextra
CXXFLAGS = -O3 -march=native -Wall -Wextra
LDFLAGS = -lgmp -lm

ALL_TARGETS = ecm_factor siqs pollard_rho qs mpqs

all: $(ALL_TARGETS)

ecm_factor: ecm_factor.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) -lecm

siqs: siqs.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

pollard_rho: pollard_rho.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

qs: qs.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

mpqs: mpqs.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(ALL_TARGETS)

.PHONY: all clean
