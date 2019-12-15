all:
	nvcc bitonic.cu FMIndex.cu -o FMIndex

run:
	@./FMIndex small.txt

clean:
	@rm -fv FMIndex
