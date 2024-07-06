symnmf: symnmf.o 
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors symnmf.o -o symnmf 

symnmf.o: symnmf.c symnmf.h
	gcc -c -ansi -Wall -Wextra -Werror -pedantic-errors symnmf.c -lm