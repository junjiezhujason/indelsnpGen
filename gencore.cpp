#include "indelsnpgen.h"

int gencore::nextseq(){
    int r;
    rid++;
    r = ref->nextseq(); 
    if(r == EOF) return r;
    helper_descputs(ref->name);

    return 0;
}

int gencore::gen(){
    // Generate the target sequence with simulated SNPs/INDELs
    // 
    pos = 0;
    int rep;
    double r;
    char c;
    bcf_hrec_t *ctg = bcf_hdr_get_hrec(h, BCF_HL_CTG, ref->name);
    if(ctg == NULL){
        fprintf(stderr, "gen cannot find CTG %s\n", ref->name);
        exit(1);
    }
    rep = bcf_hrec_find_key(ctg, "length");
    if((rep < 0) || (rep > ctg->nkeys)){
        fprintf(stderr, "gen cannot find length of %s in CTG\n", ref->name);
        exit(1);
    }
    L = atoi(ctg->vals[rep]);
    ni = 0; gi = L*ri; // number of total insertions to be simulated 
    nd = 0; gd = L*rd; // number of total deletions to be simulated
    ns = 0; gs = L*rs; // number of total snps to be simulated
    pm = 0;

    M(25);
    while(pos < L){
        c = ref->at(pos);
        r = gsl_rng_uniform(rng);
        // NEED TO REDO CALCULATION FOR SNP INDEL
        if(r < mutation_rate()){
          r = gsl_rng_uniform(rng);
          if(r < ins_rate()){ //insertion
            r = gsl_rng_uniform(rng);
            if(r < 0.96){
                rep = 10 * gsl_rng_uniform(rng) + 1;
            }else{
                rep = 4 * gsl_rng_uniform(rng) + 1;
                rep = rep * 20;
            }
            if (c !='N' && c !='n') ins(rep); // do not change N/n
            pm = pos; ni++;
            M(25);
          }else if(r < ins_rate() + del_rate()){ //deletion
            r = gsl_rng_uniform(rng);
            if(r < 0.96){
                rep = 10 * gsl_rng_uniform(rng) + 1;
            }else{
                rep = 4 * gsl_rng_uniform(rng) + 1;
                rep = rep * 20;
            }
            if (c !='N' && c !='n') del(rep); // do not change N/n
            pm = pos; nd++;
            M(25);
          }else{ //mutation
            if (c !='N' && c !='n') SNP(1); // do not change N/n
            pm = pos; ns++;
            M(25);
          }
        }else{ //match
            M(1);
        }
        ref->advanceto(pos-10); // update ref. sequence reading window
    }
    write('\0');
    return 0;
}

int gencore::ins(int k){
    // insert random sequence of length k
    // write this sequence to the target sequence fasta 
    // record the insertion event in the VCF/BCF file
    char altstr[102];
    char refstr[2];
    int i, s;
    const char **alleles = (const char**) malloc(sizeof(char*)*2);;

    b->rid = rid;
    b->pos = pos-1;

    refstr[0] = ref->at(pos-1);
    refstr[1] = '\0';
    
    for(i = 1; i <= k; i++){
        do{
            s = (double)(rand())/(double)(RAND_MAX) * nsym;
        }while(s == nsym);
        altstr[i] = sym[s];
        write(sym[s]);
    }
    altstr[0] = refstr[0];
    altstr[i] = '\0';

    alleles[0] = refstr;
    if(altstr[0]){
        alleles[1] = altstr;
    }else{
        alleles[1] = altstr+1;
    }
    bcf_update_alleles(h, b, alleles, 2);
    bcf_write1(BCF, h, b);

    return 0;
}

int gencore::del(int k){
    // delete sequence of length k; DO NOT write this sequence to the
    // target sequence fasta; record the event in the VCF/BCF file
    char altstr[2];
    char refstr[102];
    int i;
    const char **alleles = (const char**) malloc(sizeof(char*)*2);;

    b->rid = rid;
    b->pos = pos-1;

    altstr[0] = ref->at(pos-1);
    altstr[1] = '\0';
    
    for(i = 1; i <= k; i++){
        refstr[i] = ref->at(pos);
        pos++;
    }
    refstr[0] = altstr[0];
    refstr[i] = '\0';

    alleles[1] = altstr;
    if(refstr[0]){
        alleles[0] = refstr;
    }else{
        alleles[0] = refstr+1;
    }
    bcf_update_alleles(h, b, alleles, 2);
    bcf_write1(BCF, h, b);

    return 0;
}

int gencore::M(int k){
    int i;
    char c;
    for(i = 0; i < k; i++){
        c = ref->at(pos);
        if(c == '\0') break;
        write(c);
        pos++;
    }
    return 0;
}

int gencore::SNP(int k){
    char altstr[2];
    char refstr[2];
    int i;
    char snp, c;
    const char **alleles = (const char**) malloc(sizeof(char*)*2);;
    for(i = 0; i < k; i++){
        c = ref->at(pos);
        if(c == '\0') break;
        // lowercase->uppercase correction
        if(c == 'a') c = 'A';
        if(c == 'c') c = 'C';
        if(c == 'g') c = 'G';
        if(c == 't') c = 'T';
        do{
            snp = sym[(int)(gsl_rng_uniform(rng)*nsym)];
        }while(snp == c);

        b->rid = rid;
        b->pos = pos;

        refstr[0] = c;
        refstr[1] = '\0';
        altstr[0] = snp;
        altstr[1] = '\0';
        alleles[0] = refstr;
        alleles[1] = altstr;

        bcf_update_alleles(h, b, alleles, 2);
        bcf_write1(BCF, h, b);

        write(snp);
        pos++;
    }
    return 0;
}

int gencore::write(char c){
    // write a given base to the output fasta file
    static char temp[81];
    static int x = 0;
    if(c){
        temp[x] = c;
        x++;
    }
    if((x == 80) || (c == '\0')){
        temp[x] = '\0';
        helper_write(temp);
        x = 0;
    }
    return 0;
}

double gencore::mutation_rate(){
    // compute the probablity of mutation at a given position
    // FIX THIS
    if(pos - pm <= 25) return 0;
    if(L - pos - 25*(gi+gd+gs-ni-nd-ns) <= 0) return 0;
    return (double)(gi+gd+gs-ni-nd-ns)/(L-pos-24*(gi+gd+gs-ni-nd-ns));
}

double gencore::ins_rate(){
    return (double)(gi-ni)/(gi+gd+gs-ni-nd-ns);
}

double gencore::del_rate(){
    return (double)(gd-nd)/(gi+gd+gs-ni-nd-ns);
}

double gencore::snp_rate(){
    return (double)(gs-ns)/(gi+gd+gs-ni-nd-ns);
}
