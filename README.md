
# GPUseed


Initilaize and load the submodules

```
$ git submodule init
$ git submodule update cub
$ git submodule update bwa
```

Compile bwa

```
$ cd bwa
$ make
```


Build the index

```
$ ./build_index.sh -r <compression ratio of suffix postion array>  -p <prefix of the index files> <fasta file>
```
By default compression ratio is 7 and prefix of the index files is same as the fasta file


In the GPUseed directory, run `make` to complie the GPUseed test program. Modify the makefile depending upon the GPU compute capability. The default makefile is for sm_35


Run GPUseed as follows

```
./seed_gen.exe [-k <minimum seed len>] [-s] [-p] <<prefix of the index files> <reads fasta file>

```

the default value of minimum seed length is 20. By default GPUseed computes MEMs. To comput SMEMs use `-s` option. `-p` is for printing the results.




For any issues and suugestions contact Nauman Ahmed (n.ahmed@tudelft.nl)

