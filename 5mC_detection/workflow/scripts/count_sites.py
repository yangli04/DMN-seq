import pysam
import sys

bam_file = sys.argv[1]
output_bed = sys.argv[2]

with pysam.AlignmentFile(bam_file, "rb") as bam, open(output_bed, "w") as bed:
    header_bed_line = "ref\tstart\tend\tid\tquality\tstrand\tCIGAR\n"
    bed.write(header_bed_line)
    for read in bam.fetch():
        if read.is_unmapped or not read.is_read1:  # Skip unmapped reads and non-R1 reads
            continue
        if 'S' in read.cigarstring:
            continue  # Skip this read if there's any soft clipping

        # Determine the strand
        strand = '-' if read.is_reverse else '+'

        start = read.reference_start + 1  # Convert to 1-based position

        if read.is_reverse:
            bed_line = f"{read.reference_name}\t{start - 1 + read.reference_length}\t{start - 1}\t{read.query_name}\t{read.query_sequence}\t{read.mapping_quality}\t{strand}\t{read.cigarstring}\n"
            bed.write(bed_line)
        else:
            bed_line = f"{read.reference_name}\t{start - 1}\t{start - 1 + read.reference_length}\t{read.query_name}\t{read.query_sequence}\t{read.mapping_quality}\t{strand}\t{read.cigarstring}\n"
            bed.write(bed_line)
