#!/usr/bin/env python3
import sys
import re
import gzip

def Reader(fname):
    if fname.endswith("gz"):
        with gzip.open(fname) as f:
            for line in f:
                line2 = f.readline()
                f.readline()
                line4 = f.readline()
                yield Read(line.decode(), line2.decode(), line4.decode())
    else:
        with open(fname) as f:
            for line in f:
                line2 = f.readline()
                f.readline()
                line4 = f.readline()
                yield Read(line, line2, line4)

class Read():
    def __init__(self, line, line2, line4):
        self.rname = line.strip()
        self.seq = line2.strip()
        self.qual = line4.strip()
    def get_umi(self, start, end):
        return self.seq[start: end]
    def barcode(self, start, end):
        return self.seq[start: end]
    def cut(self, start, end):
        self.seq = self.seq[start:end]
        self.qual = self.qual[start:end]
    def set_umi(self, seq):
        arr = self.rname.split()
        arr[0] += ":" + seq
        self.rname = " ".join(arr)
    # def __str__(self):
    #     return self.rname + "\n" + self.seq + "\n" + "+\n" + self.qual
    def qual_gt(self, cutoff):
        qual_sum = sum([ord(i) - 33 for i in self.qual])
        return int(qual_sum / len(self.qual)) >= cutoff
    def write(self, f, compressed):
        if compressed:
            f.write((self.rname + "\n").encode())
            f.write((self.seq + "\n").encode())
            f.write(("+\n").encode())
            f.write((self.qual + "\n").encode())
        else:
            f.write(self.rname + "\n")
            f.write(self.seq + "\n")
            f.write("+\n")
            f.write(self.qual + "\n")

def rev_com(seq):
    base = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
    res = [base[i] for i in seq[::-1]]
    return "".join(res)

def x_hamming(seq1, seq2):
    x = len(seq1) - len(seq2)
    if x <= 0:
        seq2 = rev_com(seq2[x:])
    else:
        seq1 = seq1[: -1*x]
        seq2 = rev_com(seq2)
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def main():
    fname1 = sys.argv[1]
    fname2 = sys.argv[2]
    out1 = sys.argv[3]
    out2 = sys.argv[4]

    read1 = Reader(fname1)
    read2 = Reader(fname2)

    compressed1 = out1.endswith("gz")
    compressed2 = out2.endswith("gz")
    f1 = gzip.open(out1, "w") if compressed1 else open(out1, "w")
    f2 = gzip.open(out2, "w") if compressed2 else open(out2, "w")

    r1_new_end = 0
    for r1 in read1:
        r2 = next(read2)

        b1 = r1.barcode(10, 18)
        b2 = rev_com(b1)

        if b2 in r1.seq:
            r1_new_end = r1.seq.index(b2) - 3
            r1_tail = r1.seq[r1_new_end:]
            r2_head = r2.seq[:21]
            if x_hamming(r1_tail, r2_head) < 3:
                r1.cut(21, r1_new_end)
                r2.cut(21, r1_new_end)
                # # for debug show
                # x = len(r1_tail) - 21
                # if x > 0:
                #     print(1, " " * x + r1.seq)
                #     print(2, rev_com(r2.seq))
                # else:
                #     print(1, r1.seq)
                #     print(2, " " * abs(x) + rev_com(r2.seq))
        if (len(r1.seq) > 0) and r1.qual_gt(30) and r2.qual_gt(30):
            u1 = r1.get_umi(0, 10)
            u2 = r2.get_umi(0, 10)
            r1.set_umi(u1 + "_" + u2)
            r2.set_umi(u1 + "_" + u2)
            r1.write(f1, compressed1)
            r2.write(f2, compressed2)

    f1.close()
    f2.close()

if __name__ == '__main__':
    main()
