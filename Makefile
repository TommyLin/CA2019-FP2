all:
	nvcc bitonic.cu FMIndex.cu -o FMIndex

v2:
	nvcc FMIndex2.cu -o FMIndex

run:
	@CUDA_VISIBLE_DEVICES=2 ./FMIndex small.txt

clean:
	@rm -fv FMIndex
