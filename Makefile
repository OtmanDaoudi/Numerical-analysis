FILES= main.c methods.c methods.h

all: clean compile run

compile: $(FILES)
	gcc $(FILES) -o res

run: res
	./res

clean:
	rm -f *.exe