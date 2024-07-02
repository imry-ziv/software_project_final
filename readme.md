# Project Structure

1. symnmf.py: Python interface of your code. // implemented, but pending alignment with .c signatures
2. symnmf.h: C header file.
3. symnmf.c: C interface of your code.
4. symnmfmodule.c: Python C API wrapper. 
5. analysis.py: Analyze the algorithm.
6. setup.py: The setup file.
7. Makefile: Your make script to build the C interface.
8. *.c/h: Other modules and headers per your design.



# Interface Signatures

- float[,] sym(int n, int d, float[,] points)

- float[,] ddg(int n, int d, float[,] points)

- float[,] norm(int n, int d, float[,] points)

- float[,] symnmf(int n, int k, float[,] w, float[,] h)

- input dimensions:
	- points: always n*d
	- w: n*n
	- h: n*n	
- output dimensions:
	- sym, ddg, norm: n*n
	- symnmf: n*k (n lists of length k, interpreted as n rows and k columns)

The c module will define those. The rest goes in symnmf.c (e.g. the C code needs to be able to parse goal and txt files as well, and pass them into the previous functions just like symnmf.py)

