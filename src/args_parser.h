#ifndef ARGS_PARSER_H
#define ARGS_PARSER_H

/*
#include <stdint.h>


#include "gasal.h"
*/
#include "gasal.h"
#include <string>
#include <fstream>
#include <iostream>

enum fail_type {
    NOT_ENOUGH_ARGS,
    TOO_MANY_ARGS,
    WRONG_ARG,
    WRONG_FILES,
    WRONG_ALGO
};

class Parameters{

    public: 
        Parameters(int argc, char** argv);
        ~Parameters();
        void print();
        void failure(fail_type f);
        void help();
        void parse();
        void fileopen();





        int print_out;
        int n_threads;
        //int min_seed_size;
        int max_occ;

        bool doRev;
        bool doSmem;
        bool doMem;

        std::string query_batch_fasta_filename;


        std::ifstream query_batch_fasta;


    protected:

    private:
        int argc;
        char** argv;
};


#endif
