FILES= main.c methods.c methods.h

all: compile run

compile: $(FILES)
	gcc $(FILES) -o res

run: res
	./res