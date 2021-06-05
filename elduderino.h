#ifndef ELDUDERINO_H
#define ELDUDERINO_H

#include <stdbool.h>
#include <stddef.h>




typedef struct segment_t {
    char *start;
    char *qname;
    char *rname;
    char *cigar;
    char *seq;
    char *qual;
    char *barcode;
    char *next;
    size_t len;
    size_t qname_len;
    size_t rname_len;
    size_t cigar_len;
    size_t seq_len;
    size_t qual_len;
    size_t barcode_len;
    int32_t pos; // Max size 2^31 - 1 according to sam specifications
    uint16_t flag; // Max size 2^16 - 1 according to sam specifications
    } Segment;


typedef struct readpair_t {
    Segment left;
    Segment right;
    } ReadPair;


typedef void (*dedupe_function_t)(ReadPair *family, size_t family_size, int min_family_size);



#endif
