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
    char *seq;
    char *qual;
    const char *barcode;
    const char *barcode2;
    const char *next;
    size_t len;
    size_t qname_len;
    size_t rname_len;
    size_t cigar_len;
    size_t seq_len;
    size_t qual_len;
    size_t barcode_len;
    size_t barcode2_len;
    int32_t pos; // Max size 2^31 - 1 according to sam specifications
    uint16_t flag; // Max size 2^16 - 1 according to sam specifications
    } Segment;


typedef struct readpair_t {
    Segment segment[2];
    } ReadPair;


typedef struct optical_t {
    size_t irflt_len; // length of <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:in illumina qname
    int x;
    int y;
    } Optical;


typedef struct dedupe_t {
    size_t min_family_size;
    FILE *output_file;
    char *interread_buffer; // used to create final deduplicated read
    size_t interread_buffer_len;
    char *intraread_buffer; // used to dedupe overlaps
    size_t intraread_buffer_len;
    ReadPair *readpairs; // used by dedupe_all to store readpair family members
    size_t readpair_len;
    Optical *opticals; // used by dedupe_optical to store flowcell position
    size_t optical_len;
    
    size_t *family_sizes; // used to store family size statistics
    size_t family_sizes_len;
    float sequencing_total;
    float sequencing_errors;
    float pcr_total;
    float pcr_errors;
    } Dedupe;


typedef void (*dedupe_function_t)(Dedupe *dd, ReadPair *family, size_t family_size);



#endif
