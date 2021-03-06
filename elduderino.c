#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <getopt.h>
#include <ctype.h>
#include <errno.h>
#include <sys/mman.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>

#include "elduderino.h"
#include "hash.h"
#include "mash.h"


const uint16_t UNMAPPED = 0x4;
const uint16_t MATE_UNMAPPED = 0x8;
const uint16_t REVERSE = 0x10;
const uint16_t MATE_REVERSE = 0x20;
const uint16_t READ1 = 0x40;
const uint16_t READ2 = 0x80;
const uint16_t SECONDARY = 0x100;
const uint16_t FILTERED = 0x200;
const uint16_t SUPPPLEMENTARY = 0x800;
const uint16_t NON_PRIMARY = 0x100 | 0x800; // SECONDARY | SUPPPLEMENTARY;
const uint16_t BOTH_UNMAPPED = 0x4 | 0x8; // UNMAPPED | MATE_UNMAPPED;
const uint16_t READXREVERSEXUNMAPPEDX = 0x40 | 0x80 | 0x10 | 0x20 | 0x4 | 0x8; // READ1 | READ2 | REVERSE | MATE_REVERSE | UNMAPPED | MATE_UNMAPPED;
const uint16_t READX = 0x40 | 0x80; // READ1 | READ2;

const char *CONSUMES_REF = "MDN=X";
const char *CONSUMES_READ = "MIS=X";
const char *bases = "ACGTN";

bool endswith(const char *text, const char *suffix);
const char *parse_segment(const char *sam, const char *sam_end, Segment *segment);
void segment_fprintf(Segment segment, FILE *fp);
int32_t cigar_len(const char *cigar, size_t cigar_len, const char *ops);
int cmp_cigars(const void *p1, const void *p2);
int cmp_barcodes(const void *p1, const void *p2);
int cmp_qnames(const void *p1, const void *p2);
int cmp_irflts(const void *p1, const void *p2);
int cmp_int(const void *p1, const void *p2);
void dedupe_all(Dedupe *dd, MashTable *paired, dedupe_function_t dedupe_function);
void copy_sequence_to_buffer(Dedupe *dd, ReadPair *family, size_t family_size);
void barcode_families(Dedupe *dd, ReadPair *family, size_t family_size);
void connor_families(Dedupe *dd, ReadPair *family, size_t family_size);
void cigar_family(Dedupe *dd, ReadPair *family, size_t family_size);
void tile_families(Dedupe *dd, ReadPair *family, size_t family_size);
size_t coordinate_families(Dedupe *dd, ReadPair *family, size_t family_size);
void dedupe_optical(Dedupe *dd, ReadPair *family, size_t family_size);
void dedupe_pcr(Dedupe *dd, ReadPair *family, size_t family_size);
void trim_family(Dedupe *dd, ReadPair *family, size_t family_size);
int base(char base);
void reverse(char *start, int len);
void reversecomplement(char *start, int len);
char reversebase(char base);
const char *cigar_op(const char *cigar, const char **op, int32_t *num);
void write_stats(const char *stats_filename, Dedupe *dd);
int guess_optical_distance(const char *sam,  const char *sam_end);



int main (int argc, char **argv) {
    /*
     * 
     */
    const char *input_filename = NULL, *output_filename = NULL, *stats_filename = "stats.json";
    
    int sam_fd = -1;
    size_t sam_len = 0;
    const char *sam_start = NULL, *sam_end = NULL, *sam = NULL, *next = NULL;
    
    const char *mate = NULL, *current_rname = "", *sort_check_rname = "";
    char *position = NULL;
    size_t current_rname_len = 0, sort_check_rname_len = 0, len = 0, max_position_len = 0, position_len = 0;
    int32_t max_pos = 0, max_pos2 = 0, sort_check_pos = 0, segment_begin = 0, mate_begin = 0;
    HashTable *unpaired = NULL;
    MashTable *paired = NULL, *paired2 = NULL;
    Segment segment = {0}, mate_segment = {0}, swap_segment = {0};
    ReadPair readpair = {0};
    
    bool disable_optical_duplicates = false;
    dedupe_function_t dedupe_function = cigar_family;
    Dedupe dd = {0};
    
    // variable needed by strtol
    char *endptr = NULL;
    long val = 0;
    // variables needed by getopt_long
    int option_index = 0, c = 0;
    static struct option long_options[] = {{"output", required_argument, 0, 'o'},
                                           {"stats", required_argument, 0, 's'},
                                           {"umi", required_argument, 0, 'u'},
                                           {"min-family-size", required_argument, 0, 'm'},
                                           {"optical-duplicate-distance", required_argument, 0, 'p'},
                                           {"print-family-members", required_argument, 0, 'P'},
                                           {0, 0, 0, 0}};

    // Parse optional arguments
    while (c != -1) {
        c = getopt_long(argc, argv, "o:s:u:m:p:P:", long_options, &option_index);

        switch (c) {
            case 'P':
                dd.print_family_members = optarg;
                break;
                
            case 'o':
                if ((!endswith(optarg, ".fastq")) && (strcmp(optarg, "-") != 0)) {
                    fprintf(stderr, "Error: Output file must be of type fastq\n");
                    exit(EXIT_FAILURE);
                    }
                output_filename = optarg;
                break;
                
            case 's':
                if (!endswith(optarg, ".json")) {
                    fprintf(stderr, "Error: Stats file must be of type json\n");
                    exit(EXIT_FAILURE);
                    }
                stats_filename = optarg;
                break;
            
            case 'm':
                errno = 0;
                val = strtol(optarg, &endptr, 10);
                if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == optarg)) {
                    fprintf(stderr, "Error: Invalid --min-family-size\n");
                    exit(EXIT_FAILURE);
                    }
                dd.min_family_size = (size_t)val;
                break;                
                
            case 'p':
                if (strcmp(optarg, "disable") == 0) {
                    disable_optical_duplicates = true;
                    }
                else {
                    errno = 0;
                    val = strtol(optarg, &endptr, 10);
                    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == optarg) || val > INT_MAX - 1) {
                        fprintf(stderr, "Error: Invalid --optical-duplicate-distance\n");
                        exit(EXIT_FAILURE);
                        }
                    dd.optical_duplicate_distance = (int)val + 1;
                    }
                break;
                
            case 'u':
                if (strcmp(optarg, "thruplex") == 0) {
                    dedupe_function = connor_families;
                    }
                else if (strcmp(optarg, "thruplex_hv") == 0 || strcmp(optarg, "prism") == 0) {
                    dedupe_function = barcode_families;
                    }
                else if (strcmp(optarg, "") != 0) {
                    fprintf(stderr, "Error: Unsupported umi type: %s\n", optarg);
                    exit(EXIT_FAILURE);
                    }
                break;
                
            case '?':
                // unknown option, getopt_long already printed an error message.
                exit(EXIT_FAILURE);
            }
        }
    
    if (argc - optind != 1) {
        fprintf(stderr, "Error: No input file supplied\n");
        exit(EXIT_FAILURE);
        }
    input_filename = *(argv + optind);
    if (!endswith(input_filename, ".sam")) {
        fprintf(stderr, "Error: Input file must be of type sam\n");
        exit(EXIT_FAILURE);
        }
    
    
    if ((sam_fd = open(input_filename, O_RDONLY)) == -1) {
        fprintf(stderr, "Error: Unable to open %s\n", input_filename);
        exit(EXIT_FAILURE);
        }
    
    sam_len = (size_t)lseek(sam_fd, 0, SEEK_END);
    if (sam_len == 0) {
        fprintf(stderr, "Error: Empty sam file\n");
        exit(EXIT_FAILURE);
        }
    lseek(sam_fd, 0, SEEK_SET);
    if ((sam_start = mmap(NULL, sam_len, PROT_READ, MAP_PRIVATE, sam_fd, 0)) == NULL) {
        fprintf(stderr, "Error: Unable to memory map sam file\n");
        exit(EXIT_FAILURE);
        }
    
    if (output_filename == NULL || strcmp(output_filename, "-") == 0) {
        dd.output_file = stdout;
        }
    else {
        if ((dd.output_file = fopen(output_filename, "w")) == NULL) {
            fprintf(stderr, "Error: Unable to open %s for writing\n", output_filename);
            exit(EXIT_FAILURE);
            }
        }
    
    
    unpaired = hash_new(64);
    paired = mash_new(64);
    paired2 = mash_new(64);
    
    // Move past all comments to first read
    sam_end = sam_start + sam_len;
    for (sam = sam_start; sam < sam_end; ++sam) {
        if (*sam != '@') {
            break;
            }
        for (; sam < sam_end && *sam != '\n'; ++sam);
        }
    if (sam == sam_end) {
        fprintf(stderr, "Error: Empty sam file\n");
        exit(EXIT_FAILURE);
        }
    
    if (dd.optical_duplicate_distance == 0 && !disable_optical_duplicates) {
        dd.optical_duplicate_distance = guess_optical_distance(sam, sam_end);
        }
    
    for (;sam < sam_end; sam = next) {
        next = parse_segment(sam, sam_end, &segment);
        // Sanity check to ensure that sam file is sorted by position
        if (segment.rname_len == sort_check_rname_len && memcmp(segment.rname, sort_check_rname, sort_check_rname_len) == 0) {
            if (segment.pos < sort_check_pos) {
                fprintf(stderr, "Error: Sam file must be sorted by position\n");
                exit(EXIT_FAILURE);
                }
            }
        else {
            sort_check_rname = segment.rname;
            sort_check_rname_len = segment.rname_len;
            }
        sort_check_pos = segment.pos;
        
        // Skip secondary, supplementary and completely unmapped reads
        if ((segment.flag & NON_PRIMARY) || ((segment.flag & BOTH_UNMAPPED) == BOTH_UNMAPPED)) {
            continue;
            }
        
        // Do we have a pair of reads yet? If not store this read and move on to the next
        if ((mate = hash_pop(unpaired, segment.qname, segment.qname_len, &len)) == NULL) {
            hash_put(unpaired, segment.qname, segment.qname_len, segment.qname, segment.len);
            continue;
            }
        
        parse_segment(mate, mate + len, &mate_segment);
        
        // This is needed as an unmapped read may be positioned before or after its mate depending
        // on the value of the REVERSE flag in a sorted sam
        // If one read is unmapped then both reads will share the same rname and pos therefore these
        // don't get affected by the swap
        if (mate_segment.flag & UNMAPPED) {
            swap_segment = segment;
            segment = mate_segment;
            mate_segment = swap_segment;
            }
        
        
        // Segment_begin and mate_begin are the positions of the start of a read which will be the
        // segment[1]most position if reverse complemnted
        mate_begin = mate_segment.pos;
        if (mate_segment.flag & REVERSE) {
            mate_begin += cigar_len(mate_segment.cigar, mate_segment.cigar_len, CONSUMES_REF);
            }
        
        if (segment.flag & UNMAPPED) {
            segment_begin = mate_begin;
            }
        else {
            segment_begin = segment.pos;
            if (segment.flag & REVERSE) {
                segment_begin += cigar_len(segment.cigar, segment.cigar_len, CONSUMES_REF);
                }
            }
        
        
        // 2 x max number of decimal digits in an int32 (10 digits) + 
        // max number of decimal digits in a uint16 (5 digits) +
        // 4 x '\t' field separators +
        // terminal '\0' = 30
        position_len = mate_segment.rname_len + segment.rname_len + 30;
        if (position_len > max_position_len) {
            if ((position = realloc(position, position_len)) == NULL) {
                fprintf(stderr, "Error: Unable to alllocate memory for position buffer\n");
                exit(EXIT_FAILURE);
                }
            max_position_len = position_len;
            }
        snprintf(position, position_len, "%.*s\t%010i\t%.*s\t%010i\t%05u",
                                         (int)mate_segment.rname_len, mate_segment.rname,
                                         (int)mate_begin,
                                         (int)segment.rname_len, segment.rname,
                                         (int)segment_begin,
                                         (unsigned)(segment.flag & READXREVERSEXUNMAPPEDX));
        
        readpair.segment[0] = mate_segment;
        readpair.segment[1] = segment;
        
        if (mate_begin > segment_begin) {
            segment_begin = mate_begin;
            }
        
        if (segment.rname_len == current_rname_len && memcmp(segment.rname, current_rname, current_rname_len) == 0 && segment.pos <= max_pos) {
            if (mash_get(paired, position, position_len, &len) != NULL) {
                if (mash_put(paired, position, position_len, &readpair, sizeof(ReadPair)) == -1) {
                    fprintf(stderr, "Error: Unable to to add position to paired hash table\n");
                    exit(EXIT_FAILURE);
                    }
                }
            else {
                if (mash_put(paired2, position, position_len, &readpair, sizeof(ReadPair)) == -1) {
                    fprintf(stderr, "Error: Unable to to add position to paired2 hash table\n");
                    exit(EXIT_FAILURE);
                    }
                if (segment_begin > max_pos2) {
                    max_pos2 = segment_begin;
                    }
                }
            }
        else {
            dedupe_all(&dd, paired, dedupe_function);
            mash_destroy(paired); // this is faster than recycling paired
            paired = paired2;
            max_pos = max_pos2;
            paired2 = mash_new(64);
            max_pos2 = 0;
            if (mash_put(paired, position, position_len, &readpair, sizeof(ReadPair)) == -1) {
                fprintf(stderr, "Error: Unable to to add position to paired hash table\n");
                exit(EXIT_FAILURE);
                }
            if (segment.rname_len != current_rname_len || memcmp(segment.rname, current_rname, current_rname_len) != 0) {
                current_rname = segment.rname;
                current_rname_len = segment.rname_len;
                max_pos = segment_begin;
                }
            else if (segment_begin > max_pos) {
                max_pos = segment_begin;
                }
            }
        
        }
    dedupe_all(&dd, paired, dedupe_function);
    dedupe_all(&dd, paired2, dedupe_function);
    
    if (dd.output_file != stdout) {
        fclose(dd.output_file);
        }
    
    write_stats(stats_filename, &dd);
    
    // Clean up, not really needed but allows confirmation of no memory leaks
    free(position);
    free(dd.readpairs);
    free(dd.buffer);
    free(dd.family_sizes);
    hash_destroy(unpaired);
    mash_destroy(paired);
    mash_destroy(paired2);
    }



void write_stats(const char *stats_filename, Dedupe *dd) {
    FILE *stats_file = NULL;
    char ch = '\0';
    size_t total_reads = 0, total_families = 0, i = 0;
    
    for (i = 1; i < dd->max_family_size + 1; ++i) {
        total_families += dd->family_sizes[i];
        total_reads += i *  dd->family_sizes[i];
        }
    
    if ((stats_file = fopen(stats_filename, "r+")) == NULL) {
        if ((stats_file = fopen(stats_filename, "w")) == NULL) {
            fprintf(stderr, "Error: Unable to open %s\n", stats_filename);
            exit(EXIT_FAILURE);
            }
        fprintf(stats_file, "{\n");
        }
    else {
        fseek(stats_file, 0, SEEK_END);
        if (ftell(stats_file) == 0) {
            fprintf(stats_file, "{}\n");
            }
        
        fseek(stats_file, -1, SEEK_CUR);
        while (true) {
            ch = fgetc(stats_file);
            if (ch != '}' && !isspace(ch)) {
                if (ch != '{') {
                    fprintf(stats_file, ",");
                    }
                break;
                }
            if (fseek(stats_file, -2, SEEK_CUR) == -1) {
                fprintf(stderr, "Error: Malformed stats file\n");
                exit(EXIT_FAILURE);
                }
            }
        }
    
    fprintf(stats_file, "\n    \"family_sizes\": {");
    for (i = 1; i < dd->max_family_size + 1; ++i) {
        if (dd->family_sizes[i]) {
            if (i > 1) {
                fprintf(stats_file, ",");
                }
            fprintf(stats_file, "\n        \"%i\": %f", (int)i, (float)dd->family_sizes[i] / total_families);
            }
        }
    fprintf(stats_file, "\n    },\n");
    fprintf(stats_file, "    \"mean_family_size\": %.2f,\n", (float)total_reads / total_families);
    
    fprintf(stats_file, "    \"duplicate_rate\": %.4f,\n", (float)(dd->optical_duplicates + dd->pcr_duplicates) / dd->total_reads);
    if (dd->optical_duplicate_distance != 0) {
        fprintf(stats_file, "    \"duplicate_rate_optical\": %.4f,\n", (float)dd->optical_duplicates / dd->total_reads);
        fprintf(stats_file, "    \"duplicate_rate_pcr\": %.4f,\n", (float)dd->pcr_duplicates / dd->total_reads);
        }
    
    fprintf(stats_file, "    \"sequencing_error_rate\": %.4f,\n", dd->sequencing_errors / dd->sequencing_total);
    fprintf(stats_file, "    \"pcr_error_rate\": %.4f\n", dd->pcr_errors / dd->pcr_total);
    fprintf(stats_file, "}\n");

    fclose(stats_file);
    }



void dedupe_all(Dedupe *dd, MashTable *paired, dedupe_function_t dedupe_function) {
    char *data = NULL, *key = NULL, *previous_key = NULL;
    size_t data_size = 0, key_size = 0, previous_key_len = 0, readpair_len = 0;
    uint32_t bucket = 0;
    
    while (true) {
        data = mash_popall(paired, (const void **)&key, &key_size, &data_size, &bucket);
        if (data == NULL || key_size != previous_key_len || memcmp(previous_key, key, previous_key_len) != 0) {
            if (readpair_len > 0) {

                if (dd->print_family_members != NULL) {
                    int i = 0, j = 0;
                    size_t len = 0;
                    for (i = 0; i < readpair_len; ++i) {
                        len = dd->readpairs[i].segment[0].qname_len;
                        if (memcmp(dd->readpairs[i].segment[0].qname, dd->print_family_members, len) == 0 && dd->print_family_members[len] == '\0') {
                            for (j = 0; j < 2; ++j) {
                                for (i = 0; i < readpair_len; ++i) {
                                    fprintf(dd->output_file, "%.*s", (int)dd->readpairs[i].segment[j].len, dd->readpairs[i].segment[j].qname);
                                    }
                                }
                            exit(0);
                            }
                        }
                    }
                else {
                    dedupe_function(dd, dd->readpairs, readpair_len);
                    }
                
                readpair_len = 0;
                }
            if (data == NULL) {
                break;
                }
            previous_key_len = key_size;
            previous_key = key;
            }
        
        if (++readpair_len > dd->readpair_len) {
            if ((dd->readpairs = realloc(dd->readpairs, readpair_len * sizeof(ReadPair))) == NULL) {
                fprintf(stderr, "Error: Unable to allocate memory for readpairs buffer\n");
                exit(EXIT_FAILURE);
                }
            dd->readpair_len = readpair_len;
            }
        memcpy(dd->readpairs + readpair_len - 1, data, data_size);
        }
    }



void barcode_families(Dedupe *dd, ReadPair *family, size_t family_size) {
    size_t sub_family_size = 1, i = 1;
    ReadPair *sub_family = family;
    
    if (family_size < 2) {
        cigar_family(dd, family, family_size);
        }
    else {
        qsort(family, family_size, sizeof(ReadPair), cmp_barcodes);
        for (i = 1;; ++i) {
            if (i == family_size || cmp_barcodes(sub_family, family + i) != 0) {
                cigar_family(dd, sub_family, sub_family_size);
                
                if (i == family_size) {
                    break;
                    }
                
                sub_family = family + i;
                sub_family_size = 1;
                }
            else {
                ++sub_family_size;
                }
            }
        }
    }



void connor_families(Dedupe *dd, ReadPair *family, size_t family_size) {
    size_t sub_family_size = 1;
    int i = 0, j = 0;
    bool changed = false;
    ReadPair swap_readpair = {0};

    if (family_size < 2) {
        cigar_family(dd, family, family_size);
        }
    else {
        for (i = 0; i < family_size; ++i) {
            if (family[i].segment[0].barcode == NULL || family[i].segment[0].barcode2 == NULL) {
                fprintf(stderr, "Error: Missing valid barcode tags\n");
                exit(EXIT_FAILURE);
                }
            }
        
        for (; family_size > 0; family_size -= sub_family_size) {
            sub_family_size = 1;
            do {
                changed = false;
                for (i = sub_family_size; i < family_size; ++i) {
                    for (j = 0; j < sub_family_size; ++j) {
                        if (memcmp(family[i].segment[0].barcode, family[j].segment[0].barcode, family[i].segment[0].barcode_len - family[i].segment[0].barcode2_len - 1) == 0 ||
                            memcmp(family[i].segment[0].barcode2, family[j].segment[0].barcode2, family[i].segment[0].barcode2_len) == 0) {
                            if (i > sub_family_size) {
                                swap_readpair = family[sub_family_size];
                                family[sub_family_size] = family[i];
                                family[i] = swap_readpair;
                                }
                            ++sub_family_size;
                            changed = true;
                            break;
                            }
                        }
                    }    
                } while (changed);        
            
            cigar_family(dd, family, sub_family_size);
            family += sub_family_size;
            }
        }
    }   



void cigar_family(Dedupe *dd, ReadPair *family, size_t family_size) {
    size_t sixty_percent_family_size = 0, sub_family_size = 1, i = 1;
    ReadPair *sub_family = family;
    
    if (family_size > dd->max_family_size) {
        if ((dd->family_sizes = realloc(dd->family_sizes, (family_size + 1) * sizeof(size_t))) == NULL) {
            fprintf(stderr, "Error: Unable to allocate memory for family_sizes statistics\n");
            exit(EXIT_FAILURE);
            }
        memset(dd->family_sizes + dd->max_family_size + 1, 0, (family_size - dd->max_family_size) * sizeof(size_t));
        dd->max_family_size = family_size;
        }
    ++dd->family_sizes[family_size];
    
    if (family_size > 1) {
        sixty_percent_family_size = ((family_size * 6) / 10) + !!((family_size * 6) % 10);
        qsort(family, family_size, sizeof(ReadPair), cmp_cigars);
        
        for (i = 1;; ++i) {
            if (i == family_size || cmp_cigars(sub_family, family + i) != 0) {
                if (sub_family_size >= sixty_percent_family_size) {
                    family = sub_family;
                    family_size = sub_family_size;
                    break;
                    }
                
                if (i == family_size) {
                    return;
                    }
                
                sub_family = family + i;
                sub_family_size = 1;
                }
            else {
                ++sub_family_size;
                }
            }
        }
    
    dd->total_reads += family_size;
    copy_sequence_to_buffer(dd, family, family_size);
    if (dd->optical_duplicate_distance > 0) {
        tile_families(dd, family, family_size);
        }
    else {
        trim_family(dd, family, family_size);
        }
    }



void copy_sequence_to_buffer(Dedupe *dd, ReadPair *family, size_t family_size) {
    size_t max_seq_len = 0, required_len = 0, family_size_or_one = 0;
    int i = 0, j = 0;
    char *buffer = NULL;
    
    // Cigars should have already been checked and be identical for all family members, therefore seq_len must also be the same
    // This is NOT TRUE for unmapped segments but as we only process the first family member it does not matter
    // max_seq_len = family->segment[0].seq_len > family->segment[1].seq_len ? family->segment[0].seq_len : family->segment[1].seq_len;
    
    max_seq_len = 0;
    for (j = 0; j < 2; ++j) {
        family_size_or_one = family->segment[j].flag & UNMAPPED ? family_size : 1;
        for (i = 0; i < family_size_or_one; ++i) {
            if (family[i].segment[j].seq_len > max_seq_len) {
                max_seq_len = family[i].segment[j].seq_len;
                }
            }
        }
            
    if ((required_len = max_seq_len * 4 * family_size) > dd->buffer_len) {
        free(dd->buffer);
        if ((dd->buffer = malloc(required_len)) == NULL) {
            fprintf(stderr, "Error: Unable to allocate memory for sequence buffer\n");
            exit(EXIT_FAILURE);
            }
        dd->buffer_len = required_len;
        }
    
    buffer = dd->buffer;
    for (i = 0; i < family_size; ++i) {
        for (j = 0; j < 2; ++j) {
            memcpy(buffer, family[i].segment[j].seq, family[i].segment[j].seq_len);
            family[i].segment[j].seq = buffer;
            buffer += max_seq_len;
            memcpy(buffer, family[i].segment[j].qual, family[i].segment[j].seq_len);
            family[i].segment[j].qual = buffer;
            buffer += max_seq_len;
            }
        }
    }



void tile_families(Dedupe *dd, ReadPair *family, size_t family_size) {
    int colon_count = 0;
    size_t i = 0, j = 0, sub_family_size = 1, removed = 0;
    ReadPair *sub_family = family;

    for (i = 0; i < family_size; ++i) {
        colon_count = 0;
        for (j = 0; j < family[i].segment[0].qname_len; ++j) {
            if (family[i].segment[0].qname[j] == ':' && ++colon_count == 5) {
                family[i].irflt_len = j;
                break;
                }
            }
        if (colon_count != 5) {
            fprintf(stderr, "Error: Invalid illumina read name\n");
            exit(EXIT_FAILURE);
            }
        }
    
    qsort(family, family_size, sizeof(ReadPair), cmp_irflts);
    
    for (i = 1;; ++i) {
        if (i == family_size || cmp_irflts(sub_family, family + i) != 0) {
            
            if (sub_family_size > 1) {
                if ((removed = coordinate_families(dd, sub_family, sub_family_size)) > 0) {
                    memmove(family + i - removed, family + i, (family_size - i) * sizeof(ReadPair));
                    i -= removed;
                    family_size -= removed;
                    }
                 }
            
            if (i == family_size) {
                break;
                }
            
            sub_family = family + i;
            sub_family_size = 1;
            }
        else {
            ++sub_family_size;
            }
        }
    
    trim_family(dd, family, family_size);
    }



size_t coordinate_families(Dedupe *dd, ReadPair *family, size_t family_size) {
    // This function will only be called if there are at least two tile family members
    
    size_t sub_family_size = 1, removed = 0;
    int i = 0, j = 0;
    long val = 0, x = 0, y = 0;
    bool changed = false;
    const char *endptr = NULL, *coordinate = NULL;
    ReadPair swap_readpair = {0};

    for (i = 0; i < family_size; ++i) {
        coordinate = family[i].segment[0].qname + family[i].irflt_len + 1;
        errno = 0;
        val = strtol(coordinate, (char **)&endptr, 10);
        if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == coordinate)) {
            fprintf(stderr, "Error: Invalid illumina read name x coordinate\n");
            exit(EXIT_FAILURE);
            }
        family[i].optical_x = (int)val;

        coordinate = endptr + 1;
        errno = 0;
        val = strtol(coordinate, (char **)&endptr, 10);
        if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == coordinate)) {
            fprintf(stderr, "Error: Invalid illumina read name y coordinate\n");
            exit(EXIT_FAILURE);
            }
        family[i].optical_y = (int)val;
        }

        // Print optical pixel distance to stderr - Used for development only
//         for (i = 0; i < 1; ++i) {
//             for (j = i + 1; j < family_size; ++j) {
//                 x = abs(family[i].optical_x - family[j].optical_x);
//                 y = abs(family[i].optical_y - family[j].optical_y);
//                 
//                 fprintf(stderr, "%i\t%i\n", (int)(family[i].optical_x - family[j].optical_x), (int)(family[i].optical_y - family[j].optical_y));
//                 }
//             }
    

    for (; family_size > 0; family_size -= sub_family_size) {
        sub_family_size = 1;
        do {
            changed = false;
            for (i = sub_family_size; i < family_size; ++i) {
                for (j = 0; j < sub_family_size; ++j) {
                    x = abs(family[i].optical_x - family[j].optical_x);
                    y = abs(family[i].optical_y - family[j].optical_y);
                    if (sqrt((x * x) + (y * y)) < dd->optical_duplicate_distance) {
                        if (i > sub_family_size) {
                            swap_readpair = family[sub_family_size];
                            family[sub_family_size] = family[i];
                            family[i] = swap_readpair;
                            }
                        ++sub_family_size;
                        changed = true;
                        break;
                        }
                    }
                }    
            } while (changed);        
        
        if (sub_family_size > 1) {
            dedupe_optical(dd, family, sub_family_size);
            removed += sub_family_size - 1;
            memmove(family + 1, family + sub_family_size, (sub_family_size - 1) * sizeof(ReadPair));
            }
        family += 1;
        }
    
    return removed;
    }



void dedupe_optical(Dedupe *dd, ReadPair *family, size_t family_size) {
    size_t seq_len = 0, sixty_percent_family_size = 0;
    int read = 0, counts[5] = {0}, quals[5] = {0}, b = 0, q = 0, winner = 0, i = 0, j = 0;
    Segment *first = NULL, *second = NULL;
    
    dd->optical_duplicates += family_size - 1;
    
    if (family_size == 2) {
        for (read = 0; read < 2; ++read) {
            first = &family[0].segment[read];
            second = &family[1].segment[read];
            // seq_len will be the same for both reads as they have identical cigars
            seq_len = first->seq_len;
            for (i = 0; i < seq_len; ++i) {
                // Only correct the first read as this is the one we are going to keep
                if (first->seq[i] != second->seq[i]) {
                    if (second->qual[i] > first->qual[i] + 10) {
                        first->seq[i] = second->seq[i];
                        first->qual[i] = second->qual[i];
                        }
                    else if (first->qual[i] <= second->qual[i] + 10) {
                        first->seq[i] = 'N';
                        first->qual[i] = '!';
                        }
                    }
                }
            }
        }
    else {
        sixty_percent_family_size = ((family_size * 6) / 10) + !!((family_size * 6) % 10);
        
        for (read = 0; read < 2; ++read) {
            // seq_len will be the same for all reads as they have identical cigars
            seq_len = family[0].segment[read].seq_len;
            for (i = 0; i < seq_len; ++i) {
                memset(counts, 0, 5 * sizeof(int));
                memset(quals, 0, 5 * sizeof(int));
                
                for (j = 0; j < family_size; ++j) {
                    b = base(family[j].segment[read].seq[i]);
                    ++counts[b];
                    if ((q = family[j].segment[read].qual[i] - 33) > quals[b]) {
                        quals[b] = q;
                        }
                    }
                quals[4] = 0;
                
                winner = 0;
                for (b = 1; b < 4; ++b) {
                    if (counts[b] > counts[winner]) {
                        winner = b;
                        }
                    }
                
                if (counts[winner] >= sixty_percent_family_size) {
                    family[0].segment[read].seq[i] = bases[winner];
                    family[0].segment[read].qual[i] = quals[winner] + 33;
                    }
                else {
                    family[0].segment[read].seq[i] = 'N';
                    family[0].segment[read].qual[i] = '!';
                    }
                }
            }
        }
    }



int base(char base) {
    switch (base) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return 4;
        }
    }



char reversebase(char base) {
    switch (base) {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        default:
            return 'N';
        }
    }



void trim_family(Dedupe *dd, ReadPair *family, size_t family_size) {
    int32_t lref = 0, rref = 0, lread = 0, rread = 0, num = 0;
    int i = 0, j = 0, l = 0, r = 1, swap = 0, overhang = 0, mismatches = 0;
    const char *cigar = NULL, *cigar_end = NULL, *op = NULL;
    
    // Mappend to the same reference and pointing in different directions therefore may be a concordant pair
    // otherwise skip
    if ((!(family->segment[l].flag & UNMAPPED)) && (!(family->segment[r].flag & UNMAPPED)) &&
        family->segment[l].rname_len == family->segment[r].rname_len && memcmp(family->segment[l].rname, family->segment[r].rname, family->segment[l].rname_len) == 0 &&
        (family->segment[l].flag & REVERSE) != (family->segment[r].flag & REVERSE)) {
    
    
        // Ensure the left segment has the lowest ref pos
        lref = family->segment[l].pos;
        rref = family->segment[r].pos;
        if (rref < lref) {
            l = 1;
            r = 0;
            lref = family->segment[l].pos;
            rref = family->segment[r].pos;
            }
            
        lref -= 1;
        lread = -1;
        cigar = family->segment[l].cigar;
        cigar_end = cigar + family->segment[l].cigar_len;
        for (; cigar < cigar_end;) {
            cigar = cigar_op(cigar, &op, &num);
            if (strchr(CONSUMES_REF, *op) != NULL) {
                if (lref + num > rref) {
                    num = rref - lref;
                    }
                lref += num;
                }
            
            if (strchr(CONSUMES_READ, *op) != NULL) {
                lread += num;
                }
            
            if (lref == rref) {
                break;
                }
            }
        
        if (lref == rref) {
            cigar = family->segment[r].cigar;
            cigar_end = cigar + family->segment[r].cigar_len;
            for (; cigar < cigar_end;) {
                cigar = cigar_op(cigar, &op, &num);
                if (strchr(CONSUMES_REF, *op) != NULL) {
                    break;
                    }
                
                if (strchr(CONSUMES_READ, *op) != NULL) {
                    lread -= num;
                    }
                }
            
            
            // Now ensure the left segment is correctly orientated
            if (family->segment[l].flag & REVERSE) {
                swap = l;
                l = r;
                r = swap;
                lread = -lread;
                }
            
            // lread is the position of the base in the left read that
            // overlaps the first base in the right read. If lread is less
            // than zero then there is readthrough of the right read into
            // the left umi
            if (lread < 0) {
                for (i = 0; i < family_size; ++i) {
                    family[i].segment[r].seq -= lread;
                    family[i].segment[r].seq_len += lread;
                    family[i].segment[r].qual -= lread;
                    }
                lread = 0;
                }
                
            // rread is the position of the base in the right read that
            // overlaps the last base in the left read. If rread is
            // greater than or equal to right read length then there is 
            // readthrough of the left read into the right umi
            if ((rread = family->segment[l].seq_len - lread - 1) >= family->segment[r].seq_len) {
                overhang = rread - family->segment[r].seq_len + 1;
                for (i = 0; i < family_size; ++i) {
                    family[i].segment[l].seq_len -= overhang;
                    }
                rread = family->segment[r].seq_len - 1;
                }
            
            
            for (i = 0; i < family_size; ++i) {
                for (j = 0; j <= rread; ++j) {
                    if (family[i].segment[l].seq[j + lread] != family[i].segment[r].seq[j]) {
                        ++mismatches;
                        if (family[i].segment[l].qual[j + lread] > family[i].segment[r].qual[j] + 10) {
                            family[i].segment[r].seq[j] = family[i].segment[l].seq[j + lread];
                            family[i].segment[r].qual[j] = family[i].segment[l].qual[j + lread];
                            }
                        else if (family[i].segment[r].qual[j] > family[i].segment[l].qual[j + lread] + 10) {
                            family[i].segment[l].seq[j + lread] = family[i].segment[r].seq[j];
                            family[i].segment[l].qual[j + lread] = family[i].segment[r].qual[j];
                            }
                        else {
                            family[i].segment[r].seq[j] = 'N';
                            family[i].segment[r].qual[j] = '!';
                            family[i].segment[l].seq[j + lread] = 'N';
                            family[i].segment[l].qual[j + lread] = '!';
                            }
                        }
                    }
                }

            dd->sequencing_errors += mismatches;
            dd->sequencing_total += rread;
            }
        }

    dedupe_pcr(dd, family, family_size);
   }



void dedupe_pcr(Dedupe *dd, ReadPair *family, size_t family_size) {
    size_t sixty_percent_family_size = 0, len = 0;
    char *corrected_seq = NULL, *corrected_qual = NULL;
    int i = 0, j = 0, counts[5] = {0}, quals[5] = {0}, r = 0, r1r2[2] = {0}, read = 0, b = 0, q = 0, qual = 0, winner = 0, mismatches = 0, total = 0;
    
    dd->pcr_duplicates += family_size - 1;
    
    if (((family->segment[0].flag & READX) == READ1) && ((family->segment[1].flag & READX) == READ2)) {
        r1r2[0] = 0;
        r1r2[1] = 1;
        }
    else if (((family->segment[0].flag & READX) == READ2) && ((family->segment[1].flag & READX) == READ1)) {
        r1r2[0] = 1;
        r1r2[1] = 0;
        }
    else {
        fprintf(stderr, "Error: Invalid first/last segment flags within read pair\n");
        exit(EXIT_FAILURE);
        }
    
    sixty_percent_family_size = ((family_size * 6) / 10) + !!((family_size * 6) % 10);
    
    for (r = 0; r < 2; ++r) {
        read = r1r2[r];
        
        // len will be the same for all family members.
        len = family->segment[read].seq_len;
        // Write corrected data to the first family member
        corrected_seq = family->segment[read].seq;
        corrected_qual = family->segment[read].qual;
        
        if (family_size > 1 && !(family->segment[read].flag & UNMAPPED)) {
            for (i = 0; i < len; ++i) {
                memset(counts, 0, 5 * sizeof(int));
                memset(quals, 0, 5 * sizeof(int));
                
                for (j = 0; j < family_size; ++j) {
                    b = base(family[j].segment[read].seq[i]);
                    ++counts[b];
                    quals[b] += family[j].segment[read].qual[i] - 33;
                    }
                quals[4] = 0;
                
                
                winner = 0;
                total += counts[0];
                for (b = 1; b < 4; ++b) {
                    total += counts[b];
                    if (counts[b] > counts[winner]) {
                        mismatches += counts[winner];
                        winner = b;
                        }
                    else {
                        mismatches += counts[b];
                        }
                    }
                
                if (counts[winner] >= sixty_percent_family_size) {
                    qual = 0;
                    for (q = 0; q < 4; ++q) {
                        if (q == winner) {
                            qual += quals[q];
                            }
                        else {
                            qual -= quals[q];
                            }
                        }
                    if (qual < 0) {
                        qual = 0;
                        }
                    else if (qual > 93) {
                        qual = 93;
                        }
                    
                    corrected_seq[i] = bases[winner];
                    corrected_qual[i] = qual + 33;
                    }
                else {
                    corrected_seq[i] = 'N';
                    corrected_qual[i] = '!';
                    }
                }
            
            dd->pcr_errors += mismatches;
            dd->pcr_total += total;;
            }
        
        if (family->segment[read].flag & REVERSE) {
            reversecomplement(corrected_seq, len);
            reverse(corrected_qual, len);
            }
        
        // write to file
        if (family_size >= dd->min_family_size) {
            fprintf(dd->output_file, "@%.*s XF:i:%i\n%.*s\n+\n%.*s\n", (int)family->segment[read].qname_len, family->segment[read].qname,
                                                                       (int)family_size, 
                                                                       (int)len, corrected_seq,
                                                                       (int)len, corrected_qual);
            }
        }
    }



void reverse(char *start, int len) {
    char *end = start + len - 1;
    char temp = '\0';

    while (end > start) {
        temp = *start;
        *start = *end;
        *end = temp;
        ++start;
        --end;
        }
    }



void reversecomplement(char *start, int len) {
    char *end = start + len - 1;
    char temp = '\0';

    while (end > start) {
        temp = reversebase(*start);
        *start = reversebase(*end);
        *end = temp;
        ++start;
        --end;
        }
    if (start == end) {
        *start = reversebase(*end);
        }
    }



int cmp_qnames(const void *p1, const void *p2) {
    // qnames of both reads in pair will be identical, therefore just compare segment[0]
    Segment *s1 = (Segment *)p1, *s2 = (Segment *)p2;
    size_t min_len = s1->qname_len < s2->qname_len ? s1->qname_len : s2->qname_len;
    int ret = memcmp(s1->qname, s2->qname, min_len);
    
    if (ret == 0) {
        ret = s1->qname_len - s2->qname_len;
        }
    return ret;
    }



int cmp_irflts(const void *p1, const void *p2) {
    size_t min_len = ((ReadPair *)p1)->irflt_len < ((ReadPair *)p2)->irflt_len ? ((ReadPair *)p1)->irflt_len : ((ReadPair *)p2)->irflt_len;
    int ret = memcmp(((Segment *)p1)->qname, ((Segment *)p2)->qname, min_len);
    
    if (ret == 0) {
        ret = ((ReadPair *)p1)->irflt_len - ((ReadPair *)p2)->irflt_len;
        }
    return ret;
    }



int cmp_barcodes(const void *p1, const void *p2) {
    // barcodes of both reads in pair will be identical according to specificatin, therefore just compare segment[0]
    Segment *s1 = (Segment *)p1, *s2 = (Segment *)p2;
    size_t min_len = s1->barcode_len < s2->barcode_len ? s1->barcode_len : s2->barcode_len;
    int ret = memcmp(s1->barcode, s2->barcode, min_len);
    
    if (ret == 0) {
        ret = s1->barcode_len - s2->barcode_len;
        }
    return ret;
    }



int cmp_cigars(const void *p1, const void *p2) {
    ReadPair *r1 = (ReadPair *)p1, *r2 = (ReadPair *)p2;
    size_t len = 0;
    int ret = 0;
    
    len = r1->segment[0].cigar_len < r2->segment[0].cigar_len ? r1->segment[0].cigar_len : r2->segment[0].cigar_len;
    ret = memcmp(r1->segment[0].cigar, r2->segment[0].cigar, len);
    if (ret == 0) {
        ret = r1->segment[0].cigar_len - r2->segment[0].cigar_len;
        if (ret == 0) {
            len = r1->segment[1].cigar_len < r2->segment[1].cigar_len ? r1->segment[1].cigar_len : r2->segment[1].cigar_len;
            ret = memcmp(r1->segment[1].cigar, r2->segment[1].cigar, len);
            if (ret == 0) {
                ret = r1->segment[1].cigar_len - r2->segment[1].cigar_len;
                }
            }
        }
    return ret;
    }



const char *parse_segment(const char *read, const char *sam_end, Segment *segment) {
    int column = 0;
    const char *start = NULL, *endptr = NULL, *beginning = read;
    long val = 0;
    
    // optional fields so must be initialised
    segment->barcode = NULL; 
    segment->barcode2 = NULL;
    segment->barcode_len = 0;
    segment->barcode2_len = 0;
    
    for (start = read; read < sam_end; ++read) {
        if (*read == '\t' || *read == '\n') {
            switch (++column) {
                case 1: // qname
                    segment->qname = start;
                    segment->qname_len = read - start;                    
                    break;
                case 2: // flag
                    errno = 0;
                    val = strtol(start, (char **)&endptr, 10);
                    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == start)) {
                        fprintf(stderr, "Error: Invalid flag in sam file\n");
                        exit(EXIT_FAILURE);
                        }
                    segment->flag = (uint16_t)val;
                    break;
                case 3: // rname
                    segment->rname = start;
                    segment->rname_len = read - start;
                    break;
                case 4: // pos
                    errno = 0;
                    val = strtol(start, (char **)&endptr, 10);
                    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == start)) {
                        fprintf(stderr, "Error: Invalid pos in sam file\n");
                        exit(EXIT_FAILURE);
                        }
                    segment->pos = (int32_t)val;
                    break;
                case 6: // cigar
                    segment->cigar = start;
                    segment->cigar_len = read - start;
                    break;
                case 10: // seq
                    segment->seq = (char *)start;
                    segment->seq_len = read - start;
                    break;
                case 11: // qual
                    segment->qual = (char *)start;
                    // seq and qual must be the same length
                    if (segment->seq_len != read - start) {
                        fprintf(stderr, "Error: Sequence and quality differ in length\n");
                        exit(EXIT_FAILURE);
                        }
                    break;
                default:
                    if (column > 11) { // barcode tag
                        if (memcmp(start, "RX:Z:", 5) == 0) {
                            segment->barcode = start + 5;
                            segment->barcode_len = read - segment->barcode;
                            segment->barcode2 = memchr(segment->barcode, '-', segment->barcode_len) + 1;
                            segment->barcode2_len = read - segment->barcode2;
                            }
                        }
                }
            start = read + 1;
            if (*read == '\n') {
                segment->len = read - beginning + 1; // includes final \n
                break;
                }
            }
        }
    if (column < 11) {
        fprintf(stderr, "Error: Truncated sam file\n");
        exit(EXIT_FAILURE);
        }
    
    
    // seq and cigar must be the same length except for unmapped segment (cigar = *)
    if (memcmp(segment->cigar, "*\t", 2) != 0 && segment->seq_len != cigar_len(segment->cigar, segment->cigar_len, CONSUMES_READ)) {
        fprintf(stderr, "Error: Sequence and cigar differ in length\n");
        exit(EXIT_FAILURE);
        }
    
    return read + 1;
    }



bool endswith(const char *text, const char *suffix) {
    int offset = strlen(text) - strlen(suffix);
    return (offset >= 0 && strcmp(text + offset, suffix) == 0);
    }



int32_t cigar_len(const char *cigar, size_t cigar_len, const char *ops) {
    long val = 0;
    int32_t len = 0;
    const char *endptr = NULL, *cigar_end = cigar + cigar_len;

    if (*cigar == '*') {
        return 0;
        }
    
    for (; cigar < cigar_end; cigar = endptr + 1) {
        errno = 0;
        val = strtol(cigar, (char **)&endptr, 10);
        if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == cigar)) {
            fprintf(stderr, "Error: Invalid cigar string %.*s\n", (int)cigar_len, cigar);
            exit(EXIT_FAILURE);
            }
            
        if (strchr(ops, *endptr) != NULL) {
            len += (int32_t)val;
            }
        }
    return len;
    }



const char *cigar_op(const char *cigar, const char **op, int32_t *num) {
    long val = 0;
    const char *endptr = NULL;

    if (*cigar == '*') {
        fprintf(stderr, "Error: Missing cigar string\n");
        exit(EXIT_FAILURE);
        }
    
    errno = 0;
    val = strtol(cigar, (char **)&endptr, 10);
    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == cigar)) {
        fprintf(stderr, "Error: Invalid cigar string\n");
        exit(EXIT_FAILURE);
        }
    
    *num = (int32_t)val;
    *op = endptr;
    return endptr + 1;
    }



int guess_optical_distance(const char *sam,  const char *sam_end) {
    // sam has already been moved past header before this function is called
    int colon_count = 0, *x_coords = NULL, n = 0, optical_duplicate_distance = 0;
    size_t x_coords_len = 0, i = 0, start = 0;
    bool qname = true;
    const char *preceeding_colon = NULL, *endptr = NULL;
    long val = 0;

    if ((x_coords = calloc(1000, sizeof(int))) == NULL) {
        fprintf(stderr, "Error: Unable to allocate memory for coordinate buffer\n");
        exit(EXIT_FAILURE);
        }
    
    for (; sam < sam_end; ++sam) {
        if (qname) {
            if (*sam == ':') {
                if (++colon_count == 6) {
                    errno = 0;
                    val = strtol(preceeding_colon + 1, (char **)&endptr, 10);
                    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == preceeding_colon + 1)) {
                        // Not an illumina qname
                        break;
                        }
                    x_coords[x_coords_len] = (int)val;
                    if (++x_coords_len == 1000) {
                        break;
                        }
                    }
                preceeding_colon = sam;
                }
            else if (*sam == '\t') {
                qname = false;
                }
            }
        else if (*sam == '\n') {
            if (colon_count < 6) {
                // Not an illumina qname
                break;
                }
            qname = true;
            colon_count = 0;
            }
        }
    
    
    if (x_coords_len == 1000) {
        qsort(x_coords, x_coords_len, sizeof(int), cmp_int);
        
        for (i = 1; i < x_coords_len; ++i) {
            x_coords[i] -= *x_coords;
            }
            
        for (start = 1; start < x_coords_len; ++start) {
            if (x_coords[start]) {
                break;
                }
            }
        
        for (n = 2; n <= x_coords[start]; ++n) {
            for (i = start; i < x_coords_len; ++i) {
                if (x_coords[i] % n) {
                    break;
                    }
                }
            if (i == x_coords_len) {
                break;
                }
            }

        if (n <= x_coords[start]) {
            fprintf(stderr, "elduderino: Illumina query name. Patterned flow cell identified. Optical duplicate distance set to 2500\n");
            optical_duplicate_distance = 2501;
            }
        else {
            fprintf(stderr, "elduderino: Illumina query name. Non-patterned flow cell identified. Optical duplicate distance set to 100\n");
            optical_duplicate_distance = 101;
            }
        }
    else if (sam == sam_end) {
        fprintf(stderr, "elduderino: Illumina query name. Too few reads to identify flow cell type. Optical duplicate detection disabled\n");
        }
    else {
        fprintf(stderr, "elduderino: Non-Illumina query name. Optical duplicate detection disabled\n");
        }

    free(x_coords);
    return optical_duplicate_distance;
    }



int cmp_int(const void *p1, const void *p2) { 
    return *(const int *)p1  - *(const int *)p2; 
    } 



void segment_fprintf(Segment segment, FILE *fp) {
    fprintf(fp, "QNAME:    %.*s\n", (int)segment.qname_len, segment.qname);
    fprintf(fp, "FLAG:    ");
    if (segment.flag & UNMAPPED) fprintf(fp, " UNMAPPED");
    if (segment.flag & MATE_UNMAPPED) fprintf(fp, " MATE_UNMAPPED");
    if (segment.flag & REVERSE) fprintf(fp, " REVERSE");
    if (segment.flag & MATE_REVERSE) fprintf(fp, " MATE_REVERSE");
    if (segment.flag & READ1) fprintf(fp, " READ1");
    if (segment.flag & READ2) fprintf(fp, " READ2");
    if (segment.flag & SECONDARY) fprintf(fp, " SECONDARY");
    if (segment.flag & FILTERED) fprintf(fp, " FILTERED");
    if (segment.flag & SUPPPLEMENTARY) fprintf(fp, " SUPPPLEMENTARY");
    fprintf(fp, "\n");
    fprintf(fp, "RNAME:    %.*s\n", (int)segment.rname_len, segment.rname);
    fprintf(fp, "POS:      %i\n", (int)segment.pos);
    fprintf(fp, "CIGAR:    %.*s\n", (int)segment.cigar_len, segment.cigar);
    fprintf(fp, "SEQ:      %.*s\n", (int)segment.seq_len, segment.seq);
    fprintf(fp, "QUAL:     %.*s\n", (int)segment.seq_len, segment.qual);
    if (segment.barcode != NULL) fprintf(fp, "BARCODE:  %.*s\n", (int)segment.barcode_len, segment.barcode);
    if (segment.barcode2 != NULL) fprintf(fp, "BARCODE2: %.*s\n", (int)segment.barcode2_len, segment.barcode2);
    }

