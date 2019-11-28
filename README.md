# CA2019-FP2
Base Code for Final Project Part-2

# Compile using:
nvcc FMIndex.cu -o FMIndex

# Run using:
./FMIndex small.txt

# Arrays that you need to fill and where final results should be stored:
int **SA_Final_student;
int **L_counts_student;
char *L_student;
int F_counts_student;

# Correctness check
checker() will be used to compare the default and the values calculated by you.
100% correctness should be ensured before submission.

# Note:
You can currently use small.txt as input file, but this can change during tests.
Length of each read will be fixed to 64 for tests also, so you can customize your CUDA kernels and shared memory according to it.
