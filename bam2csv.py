import os
import glob
import logging
import argparse
import pysam
import csv

# Configure logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    filename='ngs_analysis.log')

def read_bed_file(bed_file):
    """
    Read BED file and return a list of features.
    """
    if not bed_file or not os.path.exists(bed_file):
        return []
    
    features = []
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    feature = {
                        'chrom': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'name': parts[3] if len(parts) > 3 else 'Unknown'
                    }
                    features.append(feature)
   
    except Exception as e:
        logging.error(f"Error reading BED file: {e}")
    
    return features

def find_overlapping_features(chrom, pos, bed_features):
    """
    Find features that overlap with a given position.
    """
    overlapping = []
    for feature in bed_features:
        if (feature['chrom'] == chrom and 
            feature['start'] <= pos < feature['end']):
            overlapping.append(feature['name'])
    return overlapping

def process_bam_file(bam_path, ref_fasta, bed_file=None, use_pandas=True):
    """
    Process a single BAM file and generate comprehensive coverage analysis.
    """
    try:
        # Load BED features if provided
        bed_features = read_bed_file(bed_file)
        
        # Load reference FASTA
        ref_genome = pysam.FastaFile(ref_fasta)
        
        # Open BAM file
        bam_file = pysam.AlignmentFile(bam_path, "rb")
        
        # Collect results
        results = []
        
        # Iterate through reference sequences
        for chrom in ref_genome.references:
            chrom_seq = ref_genome.fetch(chrom)
            
            for pileupcolumn in bam_file.pileup(chrom):
                pos = pileupcolumn.pos
                
                base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'Del': 0, 'Ins': 0}
                five_prime_count = 0
                three_prime_count = 0
                
                
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        base_counts[base] += 1
                    
                    if pileupread.is_del:
                        base_counts['Del'] += 1
                    if pileupread.indel != 0:
                        base_counts['Ins'] += 1
                    
                    if pileupread.alignment.is_reverse:
                        if pileupread.alignment.reference_end == pos + 1:
                            three_prime_count += 1
                    else:
                        if pileupread.alignment.reference_start == pos:
                            five_prime_count += 1
                
                total_coverage = sum(base_counts.values())
                ref_base = chrom_seq[pos]
                
                gene_feature = 'Unknown'
                if bed_features:
                    overlapping = find_overlapping_features(chrom, pos, bed_features)
                    if overlapping:
                        gene_feature = ','.join(overlapping)
                
                # Calculate percentages for A, C, G, T, Del, and Ins
                def calculate_percentage(count):
                    return round((count / total_coverage * 100), 2) if total_coverage > 0 else 0
                
                result = {
                    'Coordinate': pos + 1,
                    'Gene': gene_feature,
                    'RefSeq': ref_base,
                    'Coverage': total_coverage,
                    'A': base_counts['A'],
                    'C': base_counts['C'],
                    'G': base_counts['G'],
                    'T': base_counts['T'],
                    'Del': base_counts['Del'],
                    'Ins': base_counts['Ins'],
                    '5_Prime_End': five_prime_count,
                    '3_Prime_End': three_prime_count,
                    'A_perc': calculate_percentage(base_counts['A']),
                    'C_perc': calculate_percentage(base_counts['C']),
                    'G_perc': calculate_percentage(base_counts['G']),
                    'T_perc': calculate_percentage(base_counts['T']),
                    'Del_perc': calculate_percentage(base_counts['Del']),
                    'Ins_perc': calculate_percentage(base_counts['Ins'])
                }
                results.append(result)
        
        return results
    
    except Exception as e:
        logging.error(f"Error processing BAM file {bam_path}: {str(e)}")
        raise

def write_csv_output(results, output_path):
    """
    Write results to a CSV file using pure Python.
    """
    try:
        headers = [
            'Coordinate', 'Gene', 'RefSeq', 'Coverage',
            'A', 'C', 'G', 'T', 'Del', 'Ins', '5_Prime_End', '3_Prime_End',
            'A_perc', 'C_perc', 'G_perc', 'T_perc', 'Del_perc', 'Ins_perc'
        ]
        with open(output_path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=headers)
            writer.writeheader()
            writer.writerows(results)
    except Exception as e:
        logging.error(f"Error writing output CSV file: {e}")

def main():
    parser = argparse.ArgumentParser(description='NGS Analysis Script')
    parser.add_argument('--input_folder', type=str, required=True, help='Input folder containing BAM files')
    parser.add_argument('--reference_fasta', type=str, required=True, help='Reference FASTA file')
    parser.add_argument('--output_folder', type=str, help='Output folder for analysis results (optional)')
    parser.add_argument('--bed_file', type=str, help='Optional BED file for gene annotation')
    args = parser.parse_args()
    
    input_folder = args.input_folder
    ref_fasta = args.reference_fasta
    output_folder = args.output_folder or input_folder
    bed_file = args.bed_file
    
    os.makedirs(output_folder, exist_ok=True)
    bam_files = glob.glob(os.path.join(input_folder, '*.bam'))
    
    if not bam_files:
        logging.warning(f"No BAM files found in {input_folder}")
        return
    
    for bam_path in bam_files:
        try:
            results = process_bam_file(bam_path, ref_fasta, bed_file)
            output_csv = os.path.join(output_folder, os.path.basename(bam_path).replace('.bam', '_analysis.csv'))
            write_csv_output(results, output_csv)
            logging.info(f"Analysis results written to {output_csv}")
        except Exception as e:
            logging.error(f"Failed to process BAM file {bam_path}: {e}")

if __name__ == "__main__":
    main()
