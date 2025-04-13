import sys
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
# Trim R1 starts from AAAG. Only left R1 starts from AAAG! 
def process_reads(r1_path, r2_path, r1_out_path, r2_out_path):
    with gzip.open(r1_path, 'rt') as r1_file, gzip.open(r2_path, 'rt') as r2_file, \
         gzip.open(r1_out_path, 'wt') as r1_out, gzip.open(r2_out_path, 'wt') as r2_out:

        r1_seqs = SeqIO.parse(r1_file, 'fastq')
        r2_seqs = SeqIO.parse(r2_file, 'fastq')

        count_processed, count_trimmed = 0, 0
        for r1_record, r2_record in zip(r1_seqs, r2_seqs):
            count_processed += 1
            if r1_record.seq.startswith("AAA"):
                # Create a new SeqRecord with the trimmed sequence and adjusted quality scores
                # trimmed_seq = r1_record.seq[3:]
                # trimmed_qualities = r1_record.letter_annotations["phred_quality"][3:]
                trimmed_record = SeqRecord(r1_record.seq[3:],
                                        id=r1_record.id,
                                        name=r1_record.name,
                                        description=r1_record.description,
                                        letter_annotations={"phred_quality": r1_record.letter_annotations["phred_quality"][3:]})
                SeqIO.write(trimmed_record, r1_out, "fastq")
                SeqIO.write(r2_record, r2_out, "fastq")
                count_trimmed += 1

        print(f"Processed {count_processed} read pairs. Trimmed {count_trimmed} R1 reads.\n AAA reads ratio: {count_trimmed/count_processed}")

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python filter_trim_reads.py <R1_input> <R2_input> <R1_output> <R2_output>")
        sys.exit(1)

    process_reads(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
