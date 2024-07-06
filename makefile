symnmf: symnmf.o 
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors symnmf.o -o symnmf -lm

symnmf.o: ccode\symnmf.c ccode\symnmf.h
	gcc -c -ansi -Wall -Wextra -Werror -pedantic-errors symnmf.c -lm