#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <getopt.h>
#include <errno.h>
#include <sys/mman.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

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


bool endswith(const char *text, const char *suffix);
Segment parse_segment(const char *sam, const char *sam_end);
void segment_fprintf(Segment segment, FILE *fp);
int32_t cigar_len(const char *cigar, size_t cigar_len, const char *ops);
int cmp_cigars(const void *p1, const void *p2);
int cmp_barcodes(const void *p1, const void *p2);
void barcode_families(Dedupe *dd, ReadPair *family, size_t family_size, size_t min_family_size);
void connor_families(Dedupe *dd, ReadPair *family, size_t family_size, size_t min_family_size);
void cigar_family(Dedupe *dd, ReadPair *family, size_t family_size, size_t min_family_size);
void dedupe_family(Dedupe *dd, ReadPair *family, size_t family_size, size_t min_family_size);
inline int base(char base);
void reverse(char *start, int len);
void reversecomplement(char *start, int len);
inline char reversebase(char base);



int main (int argc, char **argv) {
    /*
     * 
     */
    const char *input_filename = NULL, *output_filename = "deduplicated.fastq", *stats_filename = "stats.json";
    size_t min_family_size = 1;
    
    int sam_fd = -1;
    size_t sam_len = 0;
    const char *sam_start = NULL, *sam_end = NULL, *sam = NULL, *next = NULL;
    
    const char *current_rname = "", *sort_check_rname = "";
    char *mate = NULL, *position = NULL;
    size_t current_rname_len = 0, sort_check_rname_len = 0, len = 0, max_position_len = 0, position_len = 0, readpair_len = 0, max_readpair_len = 0;;
    int32_t max_pos = 0, max_pos2 = 0, sort_check_pos = 0, segment_begin = 0, mate_begin = 0;
    HashTable *unpaired = NULL;
    MashTable *paired = NULL, *paired2 = NULL;
    Segment segment = {0}, mate_segment = {0}, swap_segment = {0};
    ReadPair readpair = {0}, *readpairs = NULL;
    
    char *data = NULL, *key = NULL, *previous_key = NULL;
    size_t data_size = 0, key_size = 0, previous_key_len = 0;
    uint32_t bucket = 0;
    
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
                                           {"min-family-size", no_argument, 0, 'm'},
                                           {0, 0, 0, 0}};

    // Parse optional arguments
    while (c != -1) {
        c = getopt_long (argc, argv, "o:s:e:npa:m:l:u:", long_options, &option_index);

        switch (c) {
            case 'o':
                if (!endswith(optarg, ".sam")) {
                    fprintf(stderr, "Error: Output file must be of type sam\n");
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
                val = strtol(optarg, (char **)&endptr, 10);
                if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == optarg)) {
                    fprintf(stderr, "Error: Invalid --min-family-size\n");
                    exit(EXIT_FAILURE);
                    }
                min_family_size = (int)val;
                break;                
                
            case 'u':
                if (strcmp(optarg, "thruplex") == 0) {
                    dedupe_function = connor_families;
                    }
               else if (strcmp(optarg, "thruplex_hv") == 0 || strcmp(optarg, "prism") == 0) {
                    dedupe_function = barcode_families;
                    }
                else {
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
        fprintf(stderr, "Error: No input file supplied.\n");
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
    lseek(sam_fd, 0, SEEK_SET);
    if ((sam_start = mmap(NULL, sam_len, PROT_READ, MAP_PRIVATE, sam_fd, 0)) == NULL) {
        fprintf(stderr, "Error: Unable to memory map sam file\n");
        exit(EXIT_FAILURE);
        }
    
    if ((dd.output_file = fopen(output_filename, "w")) == NULL) {
        fprintf(stderr, "Error: Unable to open %s for writing\n", output_filename);
        exit(EXIT_FAILURE);
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
        for (; *sam != '\n'; ++sam);
        }
    
    for (;sam < sam_end; sam = next) {
        segment = parse_segment(sam, sam_end);
        next = segment.next;
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
        
        mate_segment = parse_segment(mate, mate + len);
        
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
        // terminal '\0'
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
            bucket = 0;
            key = NULL;
            readpair_len = 0;
            previous_key = NULL;
            previous_key_len = 0;
            while (true) {
                data = mash_popall(paired, (const void **)&key, &key_size, &data_size, &bucket);
                if (data == NULL || key_size != previous_key_len || memcmp(previous_key, key, previous_key_len) != 0) {
                    if (readpair_len > 0) {
                        dedupe_function(&dd, readpairs, readpair_len, min_family_size);
                        readpair_len = 0;
                        }
                    if (data == NULL) {
                        break;
                        }
                    previous_key_len = key_size;
                    previous_key = key;
                    }
                
                if (++readpair_len > max_readpair_len) {
                    if ((readpairs = realloc(readpairs, readpair_len * sizeof(ReadPair))) == NULL) {
                        fprintf(stderr, "Error: Unable to allocate memory for readpairs buffer\n");
                        exit(EXIT_FAILURE);
                        }
                    max_readpair_len = readpair_len;
                    }
                memcpy(readpairs + readpair_len - 1, data, data_size);
                }
            
            mash_destroy(paired); // this is faster than recycling the now empty paired
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
    
    
    
    //mash_summary_fprintf(paired, stderr);
    //mash_summary_fprintf(paired2, stderr);
    
    fclose(dd.output_file);
    }



void barcode_families(Dedupe *dd, ReadPair *family, size_t family_size, size_t min_family_size) {
    size_t sub_family_size = 1, i = 1;
    ReadPair *sub_family = family;
    
    if (family_size < 2) {
        cigar_family(dd, family, family_size, min_family_size);
        }
    else {
        qsort(family, family_size, sizeof(ReadPair), cmp_barcodes);
        for (i = 1;; ++i) {
            if (i == family_size || cmp_barcodes((void *)sub_family, (void *)(family + i)) != 0) {
                cigar_family(dd, sub_family, sub_family_size, min_family_size);
                
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



void connor_families(Dedupe *dd, ReadPair *family, size_t family_size, size_t min_family_size) {
    size_t sub_family_size = 1, i = 1;
    ReadPair *sub_family = family;
    
    if (family_size < 2) {
        cigar_family(dd, family, family_size, min_family_size);
        }
    else {
        for (i = 1;; ++i) {
            if (i == family_size || cmp_barcodes((void *)sub_family, (void *)(family + i)) != 0) {
                cigar_family(dd, sub_family, sub_family_size, min_family_size);
                
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



void cigar_family(Dedupe *dd, ReadPair *family, size_t family_size, size_t min_family_size) {
    size_t sixty_percent_family_size = 0, sub_family_size = 1, i = 1;
    ReadPair *sub_family = family;
    
//     if family_size > 1:
//         #collapse_optical(family)
//         family_size = len(family)
//         
//     stats["family_sizes"][family_size] += 1
    
    if (family_size > 1) {
        sixty_percent_family_size = ((family_size * 6) / 10) + !!((family_size * 6) % 10);
        if (sixty_percent_family_size > min_family_size) {
            min_family_size = sixty_percent_family_size;
            }
        qsort(family, family_size, sizeof(ReadPair), cmp_cigars);
        
        for (i = 1;; ++i) {
            if (i == family_size || cmp_cigars((void *)sub_family, (void *)(family + i)) != 0) {
                if (sub_family_size > min_family_size) {
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
    dedupe_family(dd, family, family_size, min_family_size);
    }



inline int base(char base) {
    switch (base) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case 'N':
            return 4;
        case 'a':
            return 0;
        case 'c':
            return 1;
        case 'g':
            return 2;
        case 't':
            return 3;
        default:
            return 4;
        }
    }



inline char reversebase(char base) {
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



void dedupe_family(Dedupe *dd, ReadPair *family, size_t family_size, size_t min_family_size) {
    //fprintf(stdout, "%i\n", (int)family_size);
    size_t sixty_percent_family_size = 0, len = 0;
    char *bases = "ACGTN", *seq_buffer = NULL, *qual_buffer = NULL;
    int i = 0, j = 0, counts[5] = {0}, quals[5] = {0}, r = 0, r1r2[2] = {0}, read = 0, b = 0, q = 0, qual = 0;
    
    
    if (((family->segment[0].flag & READX) == READ1) && ((family->segment[1].flag & READX) == READ2)) {
        r1r2[0] = 0;
        r1r2[1] = 1;
        }
    else if (((family->segment[0].flag & READX) == READ2) && ((family->segment[1].flag & READX) == READ1)) {
        r1r2[0] = 1;
        r1r2[1] = 0;
        }
    else {
        fprintf(stderr, "Error: Invalid first/last segment flags within read pair.\n");
        exit(EXIT_FAILURE);
        }
    
    
    sixty_percent_family_size = ((family_size * 6) / 10) + !!((family_size * 6) % 10);
    
    for (r = 0; r < 2; ++r) {
        read = r1r2[r];
        
        // len will be the same for all family members.
        len = family->segment[read].seq_len;
        if (len > dd->buffer_len) {
            free(dd->buffer);
            if ((dd->buffer = malloc((2 * len))) == NULL) {
                fprintf(stderr, "Error: Unable to allocate memory for deduplication buffer\n");
                exit(EXIT_FAILURE);
                }
            dd->buffer_len = len;
            }
        seq_buffer = dd->buffer;
        qual_buffer = dd->buffer + dd->buffer_len;
        
        
        if (family_size == 1 || (family->segment[read].flag & UNMAPPED)) {
            memcpy(seq_buffer, family->segment[read].seq, len);
            memcpy(qual_buffer, family->segment[read].qual, len);
            }
        else {
            for (i = 0; i < len; ++i) {
                memset(counts, 0, 5 * sizeof(int));
                memset(quals, 0, 5 * sizeof(int));
                
                for (j = 0; j < family_size; ++j) {
                    b = base(family[j].segment[read].seq[i]);
                    ++counts[b];
                    quals[b] += family[j].segment[read].qual[j] - 33;
                    }
                quals[4] = 0;
                
                for (b = 0; b < 4; ++b) {
                    if (counts[b] >= sixty_percent_family_size) {
                        break;
                        }
                    }
                
                qual = 0;
                for (q = 0; q < 4; ++q) {
                    if (q == b) {
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
                
                seq_buffer[i] = bases[b];
                qual_buffer[i] = qual + 33;
                }
            }
        
        if (family->segment[read].flag & REVERSE) {
            reversecomplement(seq_buffer, len);
            reverse(qual_buffer, len);
            }
        
        // write to file
        fwrite(family->segment[read].qname, 1, family->segment[read].qname_len, dd->output_file);
        fwrite("\n", 1, 1, dd->output_file);
        fwrite(seq_buffer, 1, len, dd->output_file);
        fwrite("\n+\n", 1, 3, dd->output_file);
        fwrite(qual_buffer, 1, len, dd->output_file);
        fwrite("\n", 1, 1, dd->output_file);
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



int cmp_barcodes(const void *p1, const void *p2) {
    ReadPair *r1 = (ReadPair *)p1, *r2 = (ReadPair *)p2;
    size_t len = 0;
    int ret = 0;
    
    // barcodes of both reads in pair will be identical according to specificatin, therefore just compare segment[0]
    len = r1->segment[0].barcode_len < r2->segment[0].barcode_len ? r1->segment[0].barcode_len : r2->segment[0].barcode_len;
    ret = memcmp(r1->segment[0].barcode, r2->segment[0].barcode, len);
    if (ret == 0) {
        ret = r1->segment[0].barcode_len - r2->segment[0].barcode_len;
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



Segment parse_segment(const char *read, const char *sam_end) {
    Segment segment;
    int column = 0;
    const char *start = NULL, *endptr = NULL, *beginning = read;
    long val = 0;
    
    segment.start = read;
    segment.barcode = NULL; // optional field so must be initialised
    segment.barcode_len = 0; // optional field so must be initialised
    for (start = read; read < sam_end; ++read) {
        if (*read == '\t' || *read == '\n') {
            switch (++column) {
                case 1: // qname
                    segment.qname = start;
                    segment.qname_len = read - start;
                    break;
                case 2: // flag
                    errno = 0;
                    val = strtol(start, (char **)&endptr, 10);
                    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == start)) {
                        fprintf(stderr, "Error: Invalid flag in sam file.\n");
                        exit(EXIT_FAILURE);
                        }
                    segment.flag = (uint16_t)val;
                    break;
                case 3: // rname
                    segment.rname = start;
                    segment.rname_len = read - start;
                    break;
                case 4: // pos
                    errno = 0;
                    val = strtol(start, (char **)&endptr, 10);
                    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == start)) {
                        fprintf(stderr, "Error: Invalid pos in sam file.\n");
                        exit(EXIT_FAILURE);
                        }
                    segment.pos = (int32_t)val;
                    break;
                case 6: // cigar
                    segment.cigar = start;
                    segment.cigar_len = read - start;
                    break;
                case 10: // seq
                    segment.seq = start;
                    segment.seq_len = read - start;
                    break;
                case 11: // qual
                    segment.qual = start;
                    segment.qual_len = read - start;
                    break;
                default:
                    if (column > 11) { // barcode tag
                        if (memcmp(start, "RX:Z:", 5) == 0) {
                            segment.barcode = start + 5;
                            segment.barcode_len = read - start - 5;
                            }
                        }
                }
            start = read + 1;
            if (*read == '\n') {
                segment.len = read - beginning + 1; // includes final \n
                segment.next = read + 1;
                break;
                }
            }
        }
    if (column < 11) {
        fprintf(stderr, "Error: Truncated sam file.\n");
        exit(EXIT_FAILURE);
        }
    
    // Seq and qual can legally differ in length if qual = "*" but that should not occur here
    if (segment.seq_len != segment.qual_len) {
        fprintf(stderr, "Error: Sequence and quality differ in length.\n");
        exit(EXIT_FAILURE);
        }
    
    // Seq and qual can legally differ in length if qual = "*" but that should not occur here
    if (memcmp(segment.cigar, "*\t", 2) != 0 && segment.seq_len != cigar_len(segment.cigar, segment.cigar_len, CONSUMES_READ)) {
        fprintf(stderr, "Error: Sequence and cigar differ in length.\n");
        exit(EXIT_FAILURE);
        }
    
    return segment;
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



void segment_fprintf(Segment segment, FILE *fp) {
    fprintf(fp, "QNAME: %.*s\n", (int)segment.qname_len, segment.qname);
    fprintf(fp, "FLAG: ");
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
    fprintf(fp, "RNAME: %.*s\n", (int)segment.rname_len, segment.rname);
    fprintf(fp, "POS:   %i\n", (int)segment.pos);
    fprintf(fp, "CIGAR: %.*s\n", (int)segment.cigar_len, segment.cigar);
    }

