To run the linpack test on your device, you need to install mpi in the Linux terminal, compile the lab code using the command:

```mpic++ ... -o test```

And run with the command:

```mpirun -np <your num of process> ./test.```

Some flags to initialize matrices (for using class Matrix):
- 'r' — random numbers.
- 'e' — the diagonal numbers are equal to one.
- 'n' — all numbers are zero.
- 't' — given conditions for padding (diagonally all numbers are 2.0 and all other numbers are 1.0).
- 'h' — all numbers are equal to the height of the matrix. 