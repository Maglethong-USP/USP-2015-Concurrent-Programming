Project 1 - Jacobi Method
=========================

Objective
---------

Objective of this project is to execute the [Jacobi Method](https://en.wikipedia.org/wiki/Jacobi_method) for linear equation systems and compare sequential and concurrent executions.


Input
-----

Input files require the following format:


```
size-of-the-system			# Number of equations in the system
row-test					# The row to be used to verify the result
error 						# The precision we are aiming for
max-iterations				# Maximum number of iterations. after this the program will exit
a00 a01 a02 ... a0n			# First line of the coefficients matrix
a10 a11 a12 ... a1n			# 	and so on
...
an0 an1 an2 ... ann
b0							# First line of the constants matrix
b1							# 	and so on
...
bn
```


Example input files can be found in `./project1/doc/inputs/`


Compilation and Execution
---------

To compile run:

	make

or

	gcc -o ./bin/prog -I ./lib/ -lm -pthread ./src/*.c


To execute run:

	./bin/prog input-file-name [ number-of-threads [ number-of-times-to-execute [ output-file-name ] ] ]

Defalut values are:
```
number-of-threads			1
number-of-times-to-execute 	1
output-file-name			none
```