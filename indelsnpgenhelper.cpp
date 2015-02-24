#include "indelsnpgen.h"

static int writepos;
//static int current;
static FILE *fout;
static int line_width;


int helper_writeinit(FILE *stream, int _line_width){
    fout = stream;
    writepos = 0;
//  current = 0;
    line_width = _line_width;
    return 0;
}

int helper_descputs(char *RNAME){
    static bool first = true;
    if(!first){
        fputc('\n', fout);
    }
    first = false;
    fputc('>', fout);
    fputs("target", fout);
    fputs(RNAME, fout);
    fputc('\n', fout);
    writepos = 0;
    return 0;
}

int helper_write(char *temp){
    int i;
    for(i = 0; temp[i]; i++){
        if(writepos == line_width) fputc('\n', fout);
        writepos %= line_width;
        fputc(temp[i], fout);
        writepos++;
    }
    return 0;
}

