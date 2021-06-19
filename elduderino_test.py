import subprocess
import os
import sys
#from collections import defaultdict



UNMAPPED =  4
MATE_UNMAPPED = 8
REVERSED = 16
MATE_REVERSED = 32
FIRST = 64
LAST = 128
SECONDARY = 256
SUPPLEMENTARY = 2048



class Read(object):
    def __init__(self, seq, cigar=None, qual=None, pos=1, rname="chr1"):
        self.seq = seq.lstrip(" ")
        self.pos = pos + len(seq) - len(self.seq)
        self.rname = rname
        self.qual = qual or ("a" * len(self.seq))
        if len(self.qual) != len(self.seq):
            sys.exit(f"{self.seq} and {self.qual} differ in length")
        self.cigar = cigar or "{}M".format(len(self.seq))
        self.tags = []
    
    def __str__(self):
        return "\t".join([self.qname, str(self.flag), self.rname, str(self.pos), "10", self.cigar, self.rnext, str(self.pnext), "0", self.seq, self.qual] + self.tags + ["\n"])
    
    def __repr__(self):
        return str(self)



class Pair(object):
    i = 0
    
    def __init__(self, read1, read2, barcode=None):
        type(self).i += 1
        read1.qname = read2.qname = "QNAME_{}".format(self.i)
        read1.flag = FIRST | MATE_REVERSED
        read2.flag = LAST | REVERSED
        read1.rnext = read2.rname
        read1.pnext = read2.pos
        read2.rnext = read1.rname
        read2.pnext = read1.pos
        
        if read1.cigar == "*":
            read1.flag = read1.flag | UNMAPPED
            read2.flag = read2.flag | MATE_UNMAPPED
        elif read2.cigar == "*":
            read2.flag = read2.flag | UNMAPPED
            read1.flag = read1.flag | MATE_UNMAPPED
        
        if barcode:
            read1.tags.append(f"RX:Z:{barcode}")
            read2.tags.append(f"RX:Z:{barcode}")
        
        self.read1 = read1
        self.read2 = read2
    
    def __str__(self):
        return "{}{}".format(self.read1, self.read2)
    
    def __repr__(self):
        return str(self)



RC_TRANS = str.maketrans("ATGC", "TACG")
def rc(seq):
    return seq.translate(RC_TRANS)[::-1]



def execute(sam, expected, umi=None, min_family_size=1):
    if os.path.exists("test.sam"):
        sys.exit("test.sam already exists")
    
    reads = []
    for pair in sam:
        reads.extend([pair.read1, pair.read2])
    
    with open("test.sam", "wt") as f:
        for read in sorted(reads, key=lambda x:(x.rname, x.pos, x.flag & REVERSED)):
            f.write(str(read))
    
    cmd = ["./elduderino", "test.sam", "--output", "-", "--min-family-size", str(min_family_size)]
    if umi:
        cmd += ["--umi", umi]
    completed = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, bufsize=1)
    
    if completed.returncode != 0:
        sys.exit(completed.stderr)
    
    os.unlink("test.sam")
    
    n = 7
    result = []
    for i, row in enumerate(completed.stdout.splitlines()):
        n = i % 8
        row = row.strip()
        if n == 0:
            family_size = int(row.split(" XF:i:")[1])
        elif n == 1:
            seq1 = row
        elif n == 3:
            qual1 = row
        elif n == 5:
            seq2 = rc(row)
        elif n == 7:
            qual2 = row[::-1]
            result.append(f"{seq1} {qual1} - {seq2} {qual2} {family_size}")
    
    if n != 7:
        sys.exit("Truncated fastq")

    if sorted(result) != sorted(expected):
        print(result)
        print(expected)
        sys.exit("Failed")



def main():
    print("Overlap perfect match")
    sam = [Pair(Read("AAATTTT"),
                Read("   TTTTCCC"))]
    expected = ["AAATTTT aaaaaaa - TTTTCCC aaaaaaa 1"]
    execute(sam, expected)
    

    print("Just touching")
    sam = [Pair(Read("AAAAAAA"),
                Read("       CCCCCCC"))]
    expected = ["AAAAAAA aaaaaaa - CCCCCCC aaaaaaa 1"]
    execute(sam, expected)
    

    print("Not overlaping")
    sam = [Pair(Read("AAAAAAA"),
                Read("        CCCCCCC"))]
    expected = ["AAAAAAA aaaaaaa - CCCCCCC aaaaaaa 1"]
    execute(sam, expected)
    

    print("Not overlaping, inverted")
    sam = [Pair(Read("        AAAAAAA"),
                Read("CCCCCCC"))]
    expected = ["AAAAAAA aaaaaaa - CCCCCCC aaaaaaa 1"]
    execute(sam, expected)
    

    print("Overlap mismatch right end, q = q")
    sam = [Pair(Read("AAAGTTT"),
                Read("   TTTTCCC"))]
    expected = ["AAANTTT aaa!aaa - NTTTCCC !aaaaaa 1"]
    execute(sam, expected)
    

    print("Overlap mismatch left end, q = q")
    sam = [Pair(Read("AAATTTA"),
                Read("   TTTTCCC"))]
    expected = ["AAATTTN aaaaaa! - TTTNCCC aaa!aaa 1"]
    execute(sam, expected)
    

    print("Overlap mismatch, q = q + 10")
    sam = [Pair(Read("AAATTTA", qual="aaaaaak"),
                Read("   TTTTCCC"))]
    expected = ["AAATTTN aaaaaa! - TTTNCCC aaa!aaa 1"]
    execute(sam, expected)
    

    print("Overlap mismatch, q = q + 11")
    sam = [Pair(Read("AAATTTA", qual="aaaaaal"),
                Read("   TTTTCCC"))]
    expected = ["AAATTTA aaaaaal - TTTACCC aaalaaa 1"]
    execute(sam, expected)
    
    
    print("Overlap mismatch, q + 10 = q")
    sam = [Pair(Read("AAATTTC"),
                Read("   TTTTCCC", qual="aaakaaa"))]
    expected = ["AAATTTN aaaaaa! - TTTNCCC aaa!aaa 1"]
    execute(sam, expected)
    
    
    print("Overlap mismatch, q + 11 = q")
    sam = [Pair(Read("AAATTTC"),
                Read("   TTTTCCC", qual="aaalaaa"))]
    expected = ["AAATTTT aaaaaal - TTTTCCC aaalaaa 1"]
    execute(sam, expected)
    
    
    print("Family size = 2")
    sam = [Pair(Read("AAATTTT"),
                Read("   TTTTCCC")),
           Pair(Read("AAATTTT"),
                Read("   TTTTCCC"))]
    expected = ["AAATTTT ~~~~~~~ - TTTTCCC ~~~~~~~ 2"]
    execute(sam, expected)
    
    
    print("Family size = 3")
    sam = [Pair(Read("AAATTTT"),
                Read("   TTTTCCC")),
           Pair(Read("AAATTTT"),
                Read("   TTTTCCC")),
          Pair(Read("AAATTTT"),
                Read("   TTTTCCC"))]
    expected = ["AAATTTT ~~~~~~~ - TTTTCCC ~~~~~~~ 3"]
    execute(sam, expected)
    
    
    print("Family size = 3, one mismatch")
    sam = [Pair(Read("AACTTTT"),
                Read("   TTTTCCC")),
           Pair(Read("AAATTTT"),
                Read("   TTTTCCC")),
          Pair(Read("AAATTTT"),
                Read("   TTTTCCC"))]
    expected = ["AAATTTT ~~a~~~~ - TTTTCCC ~~~~~~~ 3"]
    execute(sam, expected)
    
    
    print("Family size = 3, one complete mismatch")
    sam = [Pair(Read("AAATTTT"),
                Read("   TTTTCCC")),
           Pair(Read("AAATTTT"),
                Read("   TTTTCCC")),
          Pair(Read("TTTAAAA"),
               Read("   AAAAGGG"))]
    expected = ["AAATTTT aaaaaaa - TTTTCCC aaaaaaa 3"]
    execute(sam, expected)
    
    
    print("Family size = 4, two complete mismatches")
    sam = [Pair(Read("AAATTTT"),
                Read("   TTTTCCC")),
           Pair(Read("AAATTTT"),
                Read("   TTTTCCC")),
          Pair(Read("TTTAAAA"),
               Read("   AAAAGGG")),
          Pair(Read("TTTAAAA"),
               Read("   AAAAGGG"))]
    expected = ["NNNNNNN !!!!!!! - NNNNNNN !!!!!!! 4"]
    execute(sam, expected)
    
        
    print("Readthrough")
    sam = [Pair(Read("   AAAATTT"),
                Read("TTTAAAA"))]
    expected = ["AAAA aaaa - AAAA aaaa 1"]
    execute(sam, expected)
    
        
    print("Readthrough, short r1")
    sam = [Pair(Read("   AAAATTT"),
                Read("TTTAAAATTTT"))]
    expected = ["AAAATTT aaaaaaa - AAAATTTT aaaaaaaa 1"]
    execute(sam, expected)
    
        
    print("Readthrough, short r2")
    sam = [Pair(Read("CAAAATTT"),
                Read(" AAAATT"))]
    expected = ["CAAAATT aaaaaaa - AAAATT aaaaaa 1"]
    execute(sam, expected)
    
    
    print("Same barcodes")
    sam = [Pair(Read("AAATTTT"),
                Read("   TTTTCCC"), barcode="A"),
           Pair(Read("AAATTTT"),
                Read("   TTTTCCC"), barcode="A")]
    expected = ["AAATTTT ~~~~~~~ - TTTTCCC ~~~~~~~ 2"]
    execute(sam, expected)
    
    
    print("Different barcodes prism")
    sam = [Pair(Read("AAATTTT"),
                Read("   TTTTCCC"), barcode="A"),
           Pair(Read("AAATTTT"),
                Read("   TTTTCCC"), barcode="B")]
    expected = ["AAATTTT aaaaaaa - TTTTCCC aaaaaaa 1", "AAATTTT aaaaaaa - TTTTCCC aaaaaaa 1"]
    execute(sam, expected, umi="prism")
    
    
    print("Different barcodes thruplex")
    sam = [Pair(Read("AAATTTT"),
                Read("   TTTTCCC"), barcode="A"),
           Pair(Read("AAATTTT"),
                Read("   TTTTCCC"), barcode="B")]
    expected = ["AAATTTT aaaaaaa - TTTTCCC aaaaaaa 1", "AAATTTT aaaaaaa - TTTTCCC aaaaaaa 1"]
    execute(sam, expected, umi="thruplex")
    
    
    print("Different barcodes 2")
    sam = [Pair(Read("AAATTTT"),
                Read("   TTTTCCC"), barcode="A"),
           Pair(Read("AAATTTT"),
                Read("   TTTTCCC"), barcode="A"),
           Pair(Read("AAATTTT"),
                Read("   TTTTCCC"), barcode="B")]
    expected = ["AAATTTT ~~~~~~~ - TTTTCCC ~~~~~~~ 2", "AAATTTT aaaaaaa - TTTTCCC aaaaaaa 1"]
    execute(sam, expected, umi="prism")
    
    
    print("Same cigars")
    sam = [Pair(Read("AAATTTT", cigar="1I6M"),
                Read("   TTTTCCC", cigar="7M")),
           Pair(Read("AAATTTT", cigar="1I6M"),
                Read("   TTTTCCC", cigar="7M"))]
    expected = ["AAATTTT ~~~~~~~ - TTTTCCC ~~~~~~~ 2"]
    execute(sam, expected)
    
    
    print("Different cigars, 1:1")
    sam = [Pair(Read("AAATTTT", cigar="1I6M"),
                Read("   TTTTCCC", cigar="7M")),
           Pair(Read("AAATTTT", cigar="2I5M"),
                Read("   TTTTCCC", cigar="7M"))]
    expected = []
    execute(sam, expected)
    
    
    print("Different cigars, 2:1")
    sam = [Pair(Read("AAATTTT", cigar="1I6M"),
                Read("   TTTTCCC", cigar="7M")),
           Pair(Read("AAATTTT", cigar="1I6M"),
                Read("   TTTTCCC", cigar="7M")),
           Pair(Read("AAATTTT", cigar="2I5M"),
                Read("   TTTTCCC", cigar="7M"))]
    expected = ["AAATTTT ~~~~~~~ - TTTTCCC ~~~~~~~ 2"]
    execute(sam, expected)
    
    
    


    
    
    
    
    
    















if __name__ == "__main__":
    main()
