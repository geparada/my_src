#include <zlib.h>
#include <math.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define MAX_LEN 1024
int main(int argc, char *argv[])
{
    gzFile fp;
    kseq_t *seq;
    long i, cnt[MAX_LEN];
    double err[MAX_LEN], table[128];

    fp = gzopen(argv[1], "r");
    seq = kseq_init(fp);
    for (i = 0; i < 128; ++i) table[i] = pow(10, -i/10.);
    for (i = 0; i < MAX_LEN; ++i) cnt[i] = 0, err[i] = 0.;
    while (kseq_read(seq) >= 0)
        for (i = 0; i < seq->qual.l; ++i)
            ++cnt[i], err[i] += table[seq->qual.s[i] - 33]; // watch out underflow!
    for (i = 0; i < MAX_LEN && cnt[i]; ++i)
        printf("%d\t%ld\t%f\n", i, cnt[i], err[i] / cnt[i]);
    kseq_destroy(seq);
    gzclose(fp);
    return 0;
}

