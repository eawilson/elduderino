#ifndef ELDUDERINO_H
#define ELDUDERINO_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>



typedef struct segment_t {
    const char *start;
    const char *qname;
    const char *rname;
    const char *cigar;
    const char *seq;
    const char *qual;
    const char *barcode;
    const char *next;
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
    Segment segment[2];
    } ReadPair;


typedef struct dedupe_t {
    FILE *output_file;
    char *buffer;
    size_t buffer_len;
    } Dedupe;


typedef void (*dedupe_function_t)(Dedupe *dd, ReadPair *family, size_t family_size, size_t min_family_size);



#endif
