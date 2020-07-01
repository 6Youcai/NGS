#!/usr/bin/env python3
from pyensembl import EnsemblRelease
from pyensembl import Genome
import sys

class Region():
    def __init__(self, chromosome, start, end):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.gene = "Unknow"
    def get_symbol(self, data):
        gene_names = data.gene_names_at_locus(contig = self.chromosome[3:],
            position = int(self.start) + 1, # 0 based to 1 based
            end = int(self.end))
        if len(gene_names) != 0:
            self.gene = ",".join(gene_names)
    def __str__(self):
        return "\t".join([self.chromosome, self.start, self.end, "-", self.gene])

def main():
    bed = sys.argv[1]
    gtf = "Homo_sapiens.GRCh38.90.gtf"
    hg38 = Genome(reference_name='GRCh38', annotation_name='my_genome_features', gtf_path_or_url = gtf)
    hg38.index()
    with open(bed) as f:
        for line in f:
            chromosome, start, end, *left = line.strip().split()
            r = Region(chromosome, start, end)
            r.get_symbol(hg38)
            print(line.strip(), r.gene)
            # print(r)

if __name__ == "__main__":
    main()
