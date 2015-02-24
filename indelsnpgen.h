#ifndef INDELGEN_H
#define INDELGEN_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "vcf.h"
#include "hts.h"
#include "faidx.h"
#include "gsl/gsl_rng.h"

#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)<(b)?(a):(b)

//seqreader.cpp
class seqreader{
    char *buffer;
    char *temp; // \0 terminating
    int32_t n, i;
    FILE *fa;
    char getbase();
    int enlarge_buffer();
    bool _eos;
  public:
    char name[300];
    seqreader(FILE *_fa) : fa(_fa) {
        n = 512;
        buffer = (char*) malloc(n*sizeof(char)); 
        temp = (char*) malloc((n+1)*sizeof(char)); 
        // need nextseq()
    }
    bool eos(){ return _eos; };
    bool eof();
    int32_t pos_start, pos_end; // 0 base
    int nextseq(); // 0 or EOF
    int advance(int32_t d);
    int advanceto(int32_t pos);
    char* range(int32_t start, int32_t end); // inclusive
    char at(int32_t pos);
};

//gencore.cpp
class gencore{
    FILE *fin;
    FILE *fout;
    htsFile *BCF;
    seqreader *ref;
    char *sym;
    int nsym;
    int32_t pos;
    int ins(int k);
    int del(int k);
    int M(int k);
    int SNP(int k);
    int write(char c); // c = '\0' to flush
    bcf_hdr_t *h;
    bcf1_t *b;
    int32_t rid;
    gsl_rng *rng;
    int32_t L;
    double mutation_rate();
    int32_t pm;
    double ins_rate(); // given mutation
    double ri; int gi; int ni;
    double del_rate(); // given mutation
    double rd; int gd; int nd;
    double snp_rate(); // given mutation
    double rs; int gs; int ns;

  public:
    gencore(FILE *_fin, FILE *_fout, htsFile *_BCF, double _ins_rate, double _del_rate, double _snp_rate, char *_sym, bcf_hdr_t *_h){
        fin = _fin;
        ref = new seqreader(fin);
        fout = _fout;
        BCF = _BCF;
        sym = strdup(_sym);
        ri = _ins_rate;
        rd = _del_rate;
        rs = _snp_rate;
        b = bcf_init();
        h = _h;
        rid = -1;
        rng = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(rng, time(NULL));

        for(nsym = 0; sym[nsym] != 0; nsym++);
        if(nsym == 0) ri = 0.0;
    }
    ~gencore(){
        delete ref;
        free(sym);
        bcf_hdr_destroy(h);
        bcf_destroy1(b);
        gsl_rng_free(rng);
    }

    int gen();
    int nextseq(); // 0 or EOF
    int writectg();
};



//indelgenhelper.cpp
int helper_writeinit(FILE *stream, int _line_width);
int helper_descputs(char *RNAME);
int helper_write(char *_temp);


#endif
