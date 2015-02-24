#include "indelsnpgen.h"

int main(int argc, char **argv){
    FILE *fin, *fout;
    faidx_t *fai; 
    htsFile *BCF;
    bcf_hdr_t *h;
    bcf_hrec_t *ctg;
    char len[11];
    int nseq;
    double ins_rate, del_rate, snp_rate;
    char *sym;
    int line_width;
    char *d;
    int i, j;

    if(argc != 9){
        fprintf(stderr, "Usage: indelgen <in.fa> <out.fa> <indel.bcf> <ins_rate> <del_rate> <snp_rate> <sym> <line_width>\n");
        exit(1);
    }

    h = bcf_hdr_init("w");
    fai = fai_load(argv[1]);
    if(fai == NULL){
        fprintf(stderr, "Unable to index %s\n", argv[1]);
        fprintf(stderr, "Usage: indelgen <in.fa> <out.fa> <indel.bcf> <ins_rate> <del_rate> <snp_rate> <sym> <line_width>\n");
        exit(1);
    }
    nseq = faidx_fetch_nseq(fai);
    for(i = 0; i < nseq; i++){
        ctg = (bcf_hrec_t*) calloc(1,sizeof(bcf_hrec_t));

        d = faidx_seq_name(fai, i);
        sprintf(len, "%d", faidx_seq_len(fai, i));

        ctg->key = strdup("contig");
        bcf_hrec_add_key(ctg, "ID", strlen("ID"));
        bcf_hrec_set_val(ctg, ctg->nkeys-1, d, strlen(d), 0);
        bcf_hrec_add_key(ctg, "length", strlen("length"));
        bcf_hrec_set_val(ctg, ctg->nkeys-1, len, strlen(len), 0);
        bcf_hdr_add_hrec(h, ctg);
    }
    fai_destroy(fai);
    bcf_hdr_sync(h);

    fin = fopen(argv[1], "r");
    if(fin == NULL){
        fprintf(stderr, "Unable to open %s\n", argv[1]);
        fprintf(stderr, "Usage: indelgen <in.fa> <out.fa> <indel.bcf> <ins_rate> <del_rate> <sym> <line_width>\n");
        exit(1);
    }

    fout = fopen(argv[2], "w");
    if(fout == NULL){
        fprintf(stderr, "Unable to open %s\n", argv[2]);
        fprintf(stderr, "Usage: indelgen <in.fa> <out.fa> <indel.bcf> <ins_rate> <del_rate> <snp_rate> <sym> <line_width>\n");
        exit(1);
    }
    BCF = vcf_open(argv[3], "wu");
    if(BCF == NULL){
        fprintf(stderr, "Unable to open %s\n", argv[3]);
        fprintf(stderr, "Usage: indelgen <in.fa> <out.fa> <indel.bcf> <ins_rate> <del_rate> <snp_rate> <sym> <line_width>\n");
        exit(1);
    }
    ins_rate = atof(argv[4]);
    if((ins_rate < 0.0) || (ins_rate >= 1.0)){
        fprintf(stderr, "Invalid insertion rate %lf\n", ins_rate);
        fprintf(stderr, "Usage: indelgen <in.fa> <out.fa> <indel.bcf> <ins_rate> <del_rate> <snp_rate> <sym> <line_width>\n");
        exit(1);
    }

    del_rate = atof(argv[5]);
    if((del_rate < 0.0) || (del_rate > 1.0)){
        fprintf(stderr, "Invalid deletion rate %lf\n", del_rate);
        fprintf(stderr, "Usage: indelgen <in.fa> <out.fa> <indel.bcf> <ins_rate> <del_rate> <snp_rate> <sym> <line_width>\n");
        exit(1);
    }

    snp_rate = atof(argv[6]);
    if((snp_rate < 0.0) || (snp_rate > 1.0)){
        fprintf(stderr, "Invalid snp rate %lf\n", snp_rate);
        fprintf(stderr, "Usage: indelgen <in.fa> <out.fa> <indel.bcf> <ins_rate> <del_rate> <snp_rate> <sym> <line_width>\n");
        exit(1);
    }

    for(i = 0, j = 0; argv[7][i]; i++){
        d = strchr(argv[7], argv[7][i]);
        if(d >= argv[7] + j){
            argv[7][j] = argv[7][i];
            j++;
        }
    }
    if(i == 0){
        fprintf(stderr, "Invalid sym\n");
        fprintf(stderr, "Usage: indelgen <in.fa> <out.fa> <indel.bcf> <ins_rate> <del_rate> <snp_rate> <sym> <line_width>\n");
        exit(1);
    }

    argv[7][j] = 0;
    sym = (char*) malloc(j*sizeof(char));
    strncpy(sym, argv[7], j+1);

    line_width = atoi(argv[8]);
    if(line_width < 1){
        fprintf(stderr, "Invalid line width %d\n", line_width);
        fprintf(stderr, "Usage: indelgen <in.fa> <out.fa> <indel.bcf> <ins_rate> <del_rate> <snp_rate> <sym> <line_width>\n");
    }


    bcf_hdr_write(BCF, h);

    gencore gc(fin, fout, BCF, ins_rate, del_rate, snp_rate, sym, h);
    helper_writeinit(fout, line_width);

    while(gc.nextseq() != EOF){
        gc.gen();
    }

    fclose(fin);
    fclose(fout);
    bcf_close(BCF);

    return 0;
}

