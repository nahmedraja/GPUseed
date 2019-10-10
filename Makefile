default:
	nvcc -g -O3 -I ./cub/cub -lineinfo -Xcompiler -Wall -Xptxas  -Werror -rdc=true --gpu-architecture=compute_35 --gpu-code=sm_35 --ptxas-options=-v -o seed_gen.exe seed_gen.cu
clean:
	rm -rf *.exe *.o *~
