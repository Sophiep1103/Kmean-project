spkmeans: spkmeans.o spkmeans.h
	gcc -o spkmeans spkmeans.o -lm

spkmeans.o: spkmeans.c spkmeans.h
	gcc -c spkmeans.c -ansi -Wall -Wextra -Werror -pedantic-errors -lm

clean: 
	rm*.o


