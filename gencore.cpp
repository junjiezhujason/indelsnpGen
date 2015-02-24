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
    ni = 0; gi = L*ri;
    nd = 0; gd = L*rd;
    ns = 0; gs = L*rs;
    pm = 0;

    M(25);
    while(pos < L){
        c = ref->at(pos);
        if(c == 'N'){
            r = 1; // no chance of getting variant
        }else{ 
            r = gsl_rng_uniform(rng);
        }
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
            ins(rep);
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
            del(rep);
            pm = pos; nd++;
            M(25);
          }else{ //mutation
            SNP(1);
            pm = pos; ns++;
            M(25);
          }
        }else{ //match
            M(1);
        }
        ref->advanceto(pos-10);
    }
    write('\0');
    return 0;
}

int gencore::ins(int k){
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
    if(pos - pm <= 25) return 0;
    if(L - pos <= 25) return 0;
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
