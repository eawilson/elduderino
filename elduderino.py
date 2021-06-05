import pdb
import argparse
import sys
import csv
from collections import defaultdict, Counter
from itertools import chain
from multiprocessing import Process, Queue
from collections.abc import Sequence, Mapping
from math import sqrt

from .utils import save_stats, string2cigar


try:
    from contextlib import nullcontext
except ImportError: # python <3.7
    from .utils import nullcontext


QNAME = 0
FLAG = 1
RNAME = 2
POS = 3
MAPQ = 4
CIGAR = 5
RNEXT = 6
PNEXT = 7
TLEN = 8
SEQ = 9
QUAL = 10

UNMAPPED = 0x4
MATE_UNMAPPED = 0x8
RC = 0x10
MATERC = 0x20
READ1 = 0x40
READ2 = 0x80
SECONDARY = 0X100
FILTERED = 0x200
SUPPPLEMENTARY = 0x800
NON_PRIMARY = SECONDARY | SUPPPLEMENTARY
BOTH_UNMAPPED = UNMAPPED | MATE_UNMAPPED
READXRCXUNMAPPEDX = READ1 | READ2 | RC | MATERC | UNMAPPED | MATE_UNMAPPED
READX = READ1 | READ2

CONSUMES_REF = "MDN=X"
CONSUMES_READ = "MIS=X"

L = 0
R = 1

RCOMPLEMENT = str.maketrans("ATGC", "TACG")



def dedupe_process(input_queue, output_queue, stats_queue, dedupe_func, stats, **details):
    while True:
        size_family = input_queue.get()
        if size_family is None:
            break
        output = dedupe_func(size_family, stats=stats, **details)
        if output:
            output_queue.put(output)
    stats_queue.put(stats)
    


def filewriter_process(output_queue, output_file):
    with open(output_file, "wt") as f_out:
        while True:
            text = output_queue.get()
            if text is None:
                break
            f_out.write(text)



def cigar_len(cig, ops):
    return sum(num for num, op in cig if op in ops)



def collapse_optical(family):
    tile_families = defaultdict(list)
    for pair in family:
        try:
            instrument, run, flowcell, lane, tile, x, y = pair[0][QNAME].split(":")
            tile_families[(instrument, run, flowcell, lane, tile)].append((int(x), int(y), pair))
        except (TypeError, ValueError):
            # Not an illumina identifier
            return
        
    for tile_family in tile_families.values():
        if len(tile_family) == 2:
            x0, y0, pair = tile_family[0]
            x1, y1, pair = tile_family[1]
            dist = sqrt(((x0 - x1) ** 2) + ((y0 - y1) ** 2))
            print(int(dist))
    
    
    


def elduderino(input_sam,
               output_file="output.deduped.sam",
               stats_file="stats.json",
               umi="",
               min_family_size=1,
               threads=1):
    
    threads = 1
    print(f"Multithreading not ready, threads = {threads}", file=sys.stderr)
    
    
    non_primary = defaultdict(list)
    with open(input_sam, "rt") as f_in:
        for segment in f_in:
            segment = segment.split("\t")
            if not segment[0].startswith("@") and int(segment[FLAG]) & NON_PRIMARY:
                segment[SEQ] = segment[QUAL] = ""
                non_primary[segment[QNAME]].append(segment)
    
    
    if umi == "thruplex":
        dedupe_func = dedupe_umi_inexact
    elif umi in ("thruplex_hv", "prism"):
        dedupe_func = dedupe_umi_exact
    elif not umi:
        dedupe_func = dedupe
    else:
        sys.exit(f"'{umi}' is not a valid UMI type")
    
    
    details = {"min_family_size": min_family_size,
               "non_primary": non_primary}
    stats = {"family_sizes": Counter()}
    
    
    multithreaded = (threads > 1)
    with (open(output_file, "wt") if not multithreaded else nullcontext()) as f_out:
        
        if multithreaded:
            input_queue = Queue()
            output_queue = Queue()
            stats_queue = Queue()
            workers = []
            for i in range(threads - 1):
                workers.append(Process(target=dedupe_process, args=(input_queue, output_queue, stats_queue, dedupe_func, stats), kwargs=details))
            workers.append(Process(target=filewriter_process, args=(output_queue, output_file)))
            for worker in workers:
                worker.start()
            dedupe_family = input_queue.put
            write = output_queue.put
            
        else:
            def dedupe_family(size_family):
                output = dedupe_func(size_family, stats=stats, **details)
                if output:
                    f_out.write(output)
            write = f_out.write

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        unpaired = {}
        
        current_rname = ""
        paired = defaultdict(list)
        max_pos = 0
        paired2 = defaultdict(list)
        max_pos2 = 0

        with open(input_sam, "rt") as f_in:
            previous_rnames = set()
            sort_check_rname = ""
            sort_check_pos = 0
            for segment in f_in:
                flag = int(segment[FLAG])                
                rname = segment[RNAME]
                pos = int(segment[POS])
                
                

                try:
                    mate = unpaired.pop(segment[QNAME])
                except KeyError:
                    unpaired[segment[QNAME]] = segment
                    continue
                
                mate_flag = int(mate[FLAG])
                if mate_flag & UNMAPPED:
                    segment, mate = mate, segment
                    flag, mate_flag = mate_flag, flag
                
                
                
                
                
                
                
                mate_begin = int(mate[POS])
                if mate_flag & RC:
                    mate_begin += cigar_len(string2cigar(mate[CIGAR]), CONSUMES_REF) - 1
                if flag & UNMAPPED:
                    segment_begin = mate_begin
                else:
                    segment_begin = int(segment[POS])
                    if flag & RC:
                        segment_begin += cigar_len(string2cigar(segment[CIGAR]), CONSUMES_REF) - 1
                location = (mate[RNAME], mate_begin, rname, segment_begin, flag & READXRCXUNMAPPEDX)
                
                
                
                
                
                
                if mate_begin > segment_begin:
                    segment_begin = mate_begin
                if rname == current_rname and pos <= max_pos:
                    if location in paired:
                        paired[location].append([mate, segment])
                    
                    else:
                        paired2[location].append([mate, segment])
                        if segment_begin > max_pos2:
                            max_pos2 = segment_begin
                    
                else:
                    if multithreaded and not all(worker.is_alive() for worker in workers):
                        sys.exit("Worker thread unexpectedly terminated")
                    
                    for size_family in paired.values():
                        dedupe_family(size_family)
                    paired = paired2
                    max_pos = max_pos2
                    paired2 = defaultdict(list)
                    max_pos2 = 0
                    
                    paired[location].append([mate, segment])
                    if rname != current_rname:
                        current_rname = rname
                        max_pos = segment_begin
                    elif segment_begin > max_pos:
                        max_pos = segment_begin
                    
        for size_family in chain(paired.values(), paired2.values()):
            dedupe_family(size_family)


    if multithreaded:
        for worker in workers[:-1]:
            input_queue.put(None)
        for worker in workers[:-1]:
            _stats = stats_queue.get()
            for key, val in _stats.items():
                if isinstance(val, Mapping):
                    for k, v in val.items():
                        stats[key][k] += v
                else:
                    stats[key] += _stats[key]
        output_queue.put(None)
        for worker in workers:
            worker.join()
    
    sizes = stats["family_sizes"]
    total_reads = sum(x*y for x, y in sizes.items())
    total_families = sum(sizes.values())
    stats["mean_family_size"] = total_reads / total_families
    stats["duplicate_rate"] = sum((x-1)*y for x, y in sizes.items()) / total_reads
    stats["triplicate_plus_rate"] = sum(max((x-2),0)*y for x, y in sizes.items()) / total_reads
    
    save_stats(stats_file, stats)



def dedupe_umi_exact(size_family, **kwargs):
    if len(size_family) > 1:
        umi_families = defaultdict(list)
        for pair in size_family:
            for tag in pair[0][11:]:
                if tag.startswith("RX:Z:"):
                    umi_families[tag.rstrip()].append(pair)
                    break
            else:
                sys.exit("Missing RX tags")
        return "".join(dedupe(family, **kwargs) for family in umi_families.values())
    
    else:
        return dedupe(size_family, **kwargs)



def dedupe_umi_inexact(size_family, **kwargs):
    if len(size_family) > 1:
        umi_pairs = []
        for pair in size_family:
            for tag in pair[0][11:]:
                if tag.startswith("RX:Z:"):
                    l_umi, r_umi = tag[5:].rstrip().split("-")
                    umi_pairs.append((l_umi, r_umi, pair))
                    break
            else:
                sys.exit("Missing RX tags")

        families = []
        while umi_pairs:
            l_umi, r_umi, pair = umi_pairs.pop()
            families.append([pair])
            lefts = set([l_umi])
            rights = set([r_umi])
            changed = True
            while changed:
                changed = False
                remaining = []
                for umi_pair in umi_pairs:
                    l_umi, r_umi, pair = umi_pair
                    if l_umi in lefts or r_umi in rights:
                        families[-1].append(pair)
                        lefts.add(l_umi)
                        rights.add(r_umi)
                        changed = True
                    else:
                        remaining.append(umi_pair)
                umi_pairs = remaining        
        return "".join(dedupe(family, **kwargs) for family in families)
    
    else:
        return dedupe(size_family, **kwargs)



def dedupe(family, stats, min_family_size, non_primary):
    family_size = len(family)
    
    if family_size > 1:
        #collapse_optical(family)
        family_size = len(family)
        
    stats["family_sizes"][family_size] += 1
    if family_size < min_family_size:
        return ""
    
    if family_size > 1:
        f6 = family_size * 6
        min_family_size = max(min_family_size, (f6 // 10) + (f6 % 10 > 1))
        cigar_families = defaultdict(list)
        for pair in family:
            cigar_families[(tuple(pair[0][CIGAR]), tuple(pair[1][CIGAR]))].append(pair)
        family = sorted(cigar_families.values(), key=lambda x:len(x))[-1]
        family_size = len(family)
        if family_size < min_family_size:
            return ""
    
    
    if family_size > 1:
        best = Counter()
        f6 = family_size * 6
        sixty_percent = (f6 // 10) + (f6 % 10 > 1)
        seq = [[], []]
        qual = [[], []]
        mapped = (L,) if int(family[0][R][FLAG]) & UNMAPPED else (L, R)
        for lr in mapped:
            for i in range(len(family[0][lr][SEQ])):
                bases = Counter()
                quals = defaultdict(list)
                for member, pair in enumerate(family):
                    try:
                        base = pair[lr][SEQ][i]
                    except IndexError:
                        pdb.set_trace()
                    bases[base] += 1
                    phred = ord(pair[lr][QUAL][i]) - 33
                    quals[base].append(phred)
                    best[member] += phred
                
                base, count = sorted(bases.items(), key=lambda x:x[1])[-1]
                if count < sixty_percent:
                    seq[lr].append("N")
                    qual[lr].append("!")
                else:
                    seq[lr].append(base)
                    phred = sum(quals.pop(base))
                    phred -= sum(chain(*quals.values()))
                    if phred < 0:
                        phred = 0
                    elif phred > 93:
                        phred = 93
                    qual[lr].append(chr(phred + 33))
        
        read = family[sorted(best.items(), key=lambda x:x[1])[-1][0]]
        for lr in mapped:
            read[lr][SEQ] = "".join(seq[lr])
            read[lr][QUAL] = "".join(qual[lr])
        
    else:
        read = family[0]
    
    
    additional = non_primary.get(read[L][QNAME], ())
    if additional:
        read.extend(additional)
        seq = ([read[L][SEQ]], [read[R][SEQ]])
        qual = ([read[L][QUAL]], [read[R][QUAL]])
        flags = (int(read[L][FLAG]), int(read[R][FLAG]))
        for i in range(2, len(read)):
            flag = int(read[i][FLAG])
            # r1r2 = L if corresponds to left primary segment and R if corresponds to right primary segment
            r1r2 = flag & READX == flags[R] & READX
            # rc = 0 if orientated] the same way as the corresponding primaary segment and 1 if reversed
            rc = flag & RC != flags[r1r2] & RC
            
            try:
                read[i][SEQ] = seq[r1r2][rc]
            except IndexError:
                seq[r1r2].append(seq[r1r2][0][::-1].translate(RCOMPLEMENT))
                read[i][SEQ] = seq[r1r2][rc]
                if qual[r1r2][0][-1] != "\n":
                    qual[r1r2].append(qual[r1r2][0][::-1])
                else:
                    qual[r1r2].append("{}\n".format(qual[r1r2][0][len(qual[r1r2][0])-2::-1]))
            read[i][QUAL] = qual[r1r2][rc]
    
    return "".join("\t".join(segments) for segments in read)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_sam', help="Input sam file.")
    parser.add_argument("-o", "--output", help="Output file, may be sam (default) or fastq.", dest="output_file", default=argparse.SUPPRESS)
    parser.add_argument("-m", "--min-family-size", help="Minimum family size.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-s", "--stats", help="Statistics file.", dest="stats_file", default=argparse.SUPPRESS)
    parser.add_argument("-u", "--umi", help="UMI type, allowed = thruplex, thruplex_hv, prism.", default=argparse.SUPPRESS)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", type=int, default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        elduderino(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

