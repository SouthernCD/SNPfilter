import argparse
from toolbiox.lib.common.os import cmd_run
from toolbiox.lib.common.genome.genome_feature2 import read_gff_file, read_fasta_by_faidx
import os
import pysam
import time
from interlap import InterLap
from Bio.Seq import Seq


class Job(object):
    def __init__(self):
        pass

    def run_arg_parser(self):
        # argument parse
        parser = argparse.ArgumentParser(
            prog='SNPfilter',
        )

        subparsers = parser.add_subparsers(
            title='subcommands', dest="subcommand_name")

        # argparse for prepare
        parser_a = subparsers.add_parser('prepare',
                                         description='Prepare work environment\n')

        parser_a.add_argument('sample_id', type=str,
                              help='sample id')
        parser_a.add_argument('reference', type=str,
                              help='reference genome in fasta format, should be indexed by bwa')
        parser_a.add_argument('R1', type=str,
                              help='fastq file 1')
        parser_a.add_argument('R2', type=str,
                              help='fastq file 2')
        parser_a.add_argument('-t', '--threads', type=int,
                              help='num of threads', default=8)
        parser_a.add_argument('-q', '--min-MQ', type=int,
                              help='skip alignments with mapQ smaller than INT', default=20)
        parser_a.add_argument('-Q', '--min-BQ', type=int,
                              help='skip bases with baseQ/BAQ smaller than INT', default=20)
        parser_a.set_defaults(func=prepare_main)

        # argparse for qcfilter
        parser_a = subparsers.add_parser('qcfilter',
                                         description='filtering snp from bcf file with min depth and min variant frequency\n')

        parser_a.add_argument('sample_id', type=str,
                              help='sample id')
        parser_a.add_argument('input_bcf', type=str,
                              help='input bcf file')
        parser_a.add_argument('-d', '--min_depth', type=int,
                              help='min depth', default=2)
        parser_a.add_argument('-v', '--min_variant_frequency', type=float,
                              help='min variant frequency', default=0.3)
        parser_a.add_argument('-b', '--background_site', type=str,
                              help='A comma-separated list of bcf files, loci that appear in these bcf files will be filtered out', default=None)
        parser_a.set_defaults(func=qcfilter_main)

        # argparse for codefilter
        parser_a = subparsers.add_parser('codefilter',
                                         description='filtering SNPs based on whether they cause changes in coding amino acids\n')

        parser_a.add_argument('sample_id', type=str,
                              help='sample id')
        parser_a.add_argument('input_bcf', type=str,
                              help='input bcf file')
        parser_a.add_argument('reference', type=str,
                              help='reference genome in fasta format')
        parser_a.add_argument('gff_file', type=str,
                              help='gff file for reference genome')
        parser_a.set_defaults(func=codefilter_main)

        self.arg_parser = parser

        self.args = parser.parse_args()
        
        # parser.set_defaults(func=parser.print_help())

    def run(self):
        self.run_arg_parser()
        # self.args = self.arg_parser.parse_args()
        # if hasattr(self, "subcommand_name"):
        #     self.args.func(self.args)
        # else:
        #     self.arg_parser.print_help()

        args_dict = vars(self.args)

        if args_dict["subcommand_name"] == "prepare":
            prepare_main(self.args)
        elif args_dict["subcommand_name"] == "qcfilter":
            qcfilter_main(self.args)
        elif args_dict["subcommand_name"] == "codefilter":
            codefilter_main(self.args)
        else:
            self.arg_parser.print_help()


def prepare_main(args):
    map_shell_string = """TAG=$1
REF=$2
R1=$3
R2=$4
CPU=$5

bwa mem -o $TAG.sam -t $CPU $REF $R1 $R2
samtools faidx $REF
samtools view -bS -@ $CPU -o $TAG.bam $TAG.sam
samtools fixmate -@ $CPU -r -m $TAG.bam $TAG.fm.bam
samtools sort -@ $CPU -o $TAG.st.fm.bam $TAG.fm.bam
samtools markdup -@ $CPU -r $TAG.st.fm.bam $TAG.ud.st.fm.bam
# samtools mpileup -I -t DP,AD -q 20 -Q 20 -g -o $TAG.bcf -f $REF $TAG.ud.st.fm.bam
samtools mpileup -I -t DP,AD -q %d -Q %d -g -o $TAG.bcf -f $REF $TAG.ud.st.fm.bam
samtools index $TAG.ud.st.fm.bam
bcftools index $TAG.bcf

rm $TAG.sam $TAG.bam $TAG.fm.bam $TAG.st.fm.bam"""

    map_shell_string = map_shell_string % (args.min_MQ, args.min_BQ)
    map_sh_file = args.sample_id + '.map.sh'

    with open(map_sh_file, 'w') as f:
        f.write(map_shell_string)

    cmd_string = f"bash {map_sh_file} {args.sample_id} {args.reference} {args.R1} {args.R2} {args.threads}"
    # print(cmd_string)
    cmd_run(cmd_string)

    cmd_string = "rm %s" % map_sh_file
    cmd_run(cmd_string)


def get_all_site_from_bcf(bcf_file):
    site_list = []
    bcf = pysam.VariantFile(bcf_file)
    for record in bcf.fetch():
        site_list.append((record.chrom, record.pos))
    bcf.close()
    return site_list


def qcfilter_main(args):
    back_site_list = []
    if args.background_site is not None:
        back_bcf_file_list = args.background_site.split(',')
        for back_bcf_file in back_bcf_file_list:
            back_site_list.extend(get_all_site_from_bcf(back_bcf_file))

    back_site_list = set(back_site_list)

    filterd_vcf_file = '%s.d%d.v%.2f.vcf' % (
        args.sample_id, args.min_depth, args.min_variant_frequency)

    input_bcf = pysam.VariantFile(args.input_bcf)

    with open(filterd_vcf_file, 'w') as output_vcf:
        output_vcf.write(str(input_bcf.header.__str__()))

        # input_bcf = pysam.VariantFile(args.input_bcf)

        # chromosome = 'Pjv1Scaffold_0'

        start = time.time()
        num = 0
        # record_list = []
        # for record in input_bcf.fetch(chromosome):
        for record in input_bcf.fetch():
            # report processing
            if time.time() - start > 5:
                print("Processing: ", num)
                start = time.time()
            num += 1

            if (record.chrom, record.pos) in back_site_list:
                continue

            ref_base = record.ref

            alt_count_dict = {}
            for alt in record.alleles:
                if alt == '<*>':
                    continue
                alt_count_dict.setdefault(alt, 0)
                alt_count_dict[alt] += record.samples[0]['AD'][record.alleles.index(
                    alt)]

            depth = sum(alt_count_dict.values())
            variant_ratio = 1 - \
                alt_count_dict[ref_base] / depth if depth > 0 else 0

            if variant_ratio >= args.min_variant_frequency and depth > args.min_depth and depth <= 40:
                output_vcf.write(record.__str__())

    input_bcf.close()

    filterd_bcf_file = '%s.d%d.v%.2f.bcf' % (
        args.sample_id, args.min_depth, args.min_variant_frequency)

    cmd_string = f"bcftools view -O b -o {filterd_bcf_file} {filterd_vcf_file}"
    cmd_run(cmd_string)

    cmd_string = f"bcftools index {filterd_bcf_file}"
    cmd_run(cmd_string)

    cmd_string = f"rm {filterd_vcf_file}"
    cmd_run(cmd_string)


def if_on_gene(chr_id, pos, gene_interlap):
    g_id_list = []
    for s, e, g_id in gene_interlap[chr_id].find((int(pos), int(pos))):
        g_id_list.append(g_id)
    return g_id_list


def build_gene_interlap(gene_dict, ref_dict):
    gene_interlap = {chr_id: InterLap() for chr_id in ref_dict}

    for g_id in gene_dict:
        gene = gene_dict[g_id]
        gene_interlap[gene.chr_id].add((gene.start, gene.end, g_id))

    return gene_interlap


def get_coding_pos_base_map_dict(mRNA, ref_dict):
    """
    This function is used to get the coding position and base map dict for a gene.
    """
    pos_base_map_dict = {}

    for cds in mRNA.sub_features:
        if cds.type == 'CDS':
            # print(cds.id,cds.start, cds.end)
            for i in range(cds.start, cds.end + 1):
                if mRNA.strand == "+":
                    pos_base_map_dict[i] = ref_dict[mRNA.chr_id].faidx[i - 1].seq
                elif mRNA.strand == "-":
                    pos_base_map_dict[i] = ref_dict[mRNA.chr_id].faidx[i -
                                                                       1].complement.seq

    return pos_base_map_dict


def get_complement(seq):
    if seq == "A":
        return "T"
    elif seq == "T":
        return "A"
    elif seq == "C":
        return "G"
    elif seq == "G":
        return "C"


def if_on_cds_and_change_aa(var_pos, var_base, gene, ref_dict):
    """
    This function is used to judge if a variant is on the CDS of a gene and change the aa.
    """
    output_list = []

    for mRNA in gene.sub_features:
        if mRNA.type == 'mRNA':
            on_cds = False
            change_aa = False

            pos_base_map_dict = get_coding_pos_base_map_dict(mRNA, ref_dict)

            if var_pos in pos_base_map_dict:
                on_cds = True
                pos_sort_list = sorted(
                    pos_base_map_dict.keys(), reverse=mRNA.strand == "-")

                mutative_seq = ""
                wild_seq = ""

                for i in pos_sort_list:
                    if i == var_pos:
                        mutative_seq += var_base if mRNA.strand == "+" else get_complement(
                            var_base)
                    else:
                        mutative_seq += pos_base_map_dict[i]
                    wild_seq += pos_base_map_dict[i]

                mutative_aa = str(Seq(mutative_seq).translate())
                wild_aa = str(Seq(wild_seq).translate())

                if mutative_aa != wild_aa:
                    change_aa = True

            if change_aa:
                output_list.append(
                    [gene.id, mRNA.id, on_cds, change_aa, mutative_aa, wild_aa])
            else:
                output_list.append(
                    [gene.id, mRNA.id, on_cds, change_aa, None, None])

    return output_list


def codefilter_main(args):
    gff_file = args.gff_file
    reference = args.reference
    input_bcf_file = args.input_bcf
    output_tsv = args.sample_id + ".codefilter.tsv"

    ref_dict = read_fasta_by_faidx(reference)
    gff_dict = read_gff_file(gff_file)

    gene_dict = gff_dict['gene']
    gene_interlap = build_gene_interlap(gene_dict, ref_dict)

    with open(output_tsv, 'w') as f:

        f.write("\t".join(["chr_id", "pos", "ref", "alt", "depth", "freq", "gene_id",
                           "mRNA_id", "on_CDS", "AA_changed", "changed_site","wild_aa", "mutative_aa"]) + "\n")
        input_bcf = pysam.VariantFile(input_bcf_file)

        for record in input_bcf.fetch():
            alts_list = [i for i in record.alts if i != "<*>"]
            if len(alts_list) > 1:
                continue
            alts = alts_list[0]

            if record.ref == 'N':
                continue

            ref_base = record.ref

            alt_count_dict = {}
            for alt in record.alleles:
                if alt == '<*>':
                    continue
                alt_count_dict.setdefault(alt, 0)
                alt_count_dict[alt] += record.samples[0]['AD'][record.alleles.index(
                    alt)]

            depth = sum(alt_count_dict.values())
            variant_ratio = 1 - \
                alt_count_dict[ref_base] / depth if depth > 0 else 0

            gene_list = if_on_gene(record.chrom, record.pos, gene_interlap)
            if len(gene_list) > 0:
                for gene_id in gene_list:
                    gene = gene_dict[gene_id]
                    output_list = if_on_cds_and_change_aa(
                        record.pos, alts, gene, ref_dict)
                    for g_id, m_id, on_cds, change_aa, mutative_aa, wild_aa in output_list:
                        change_site_info = "None"
                        if change_aa:
                            for i in range(len(wild_aa)):
                                if wild_aa[i] != mutative_aa[i]:
                                    change_site_info = "%s%d%s" % (wild_aa[i], (i+1), mutative_aa[i])

                        f.write("\t".join([record.chrom, str(record.pos), record.ref, alts, str(depth), str(
                            variant_ratio), g_id, m_id, str(on_cds), str(change_aa), change_site_info, str(wild_aa), str(mutative_aa)]) + "\n")

        input_bcf.close()


def main():
    job = Job()
    job.run()


if __name__ == '__main__':
    main()
