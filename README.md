# CA2019-FP2
Base Code for Final Project Part-2

# Compile using:
nvcc FMIndex.cu -o FMIndex

# Run using:
./FMIndex small.txt

# Arrays that you need to fill and where final results should be stored:
SA_Final_student;

L_counts_student;

L_student;

F_counts_student;

# Correctness check
checker() will be used to compare the default and the values calculated by you.
100% correctness should be ensured before submission.

# Note:
You can currently use small.txt as input file, but this can change during tests.

Length of each read will be fixed to 64 for tests also, so you can customize your CUDA kernels and shared memory according to it.
