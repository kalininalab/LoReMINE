import os
import glob
import pandas as pd

# Set your expected genome size here (in bp)
EXPECTED_GENOME_SIZE = 5000000  # e.g., ~5 Mb for many bacteria
CHROMOSOME_SIZE_THRESHOLD = 0.95 * EXPECTED_GENOME_SIZE  # 90% of genome size

def generate_contig_circularity_info(quast_path, output_path):
    write_file = open(quast_path + '/circularity.tsv', 'x')
    write_file.write("contig_name\tcircular\tLength\n")


    with open(output_path + "/assembly_info.txt") as assembly_info_file:
        counter = 0
        for line in assembly_info_file.readlines():
            counter += 1
            if counter == 1:
                continue
            split_array = line.split()
            if split_array[3].strip() == "Y":
                write_file.write(split_array[0].strip() + '\tyes' + '\t' + str(split_array[1].strip()) + '\n')
            elif split_array[3].strip() == "N":
                write_file.write(split_array[0].strip() + '\tno' + '\t' + str(split_array[1].strip()) + '\n')

def parse_quast_report(report_path):
    """Extract N50, #contigs, and total length from QUAST report.tsv."""
    stats = {}
    with open(report_path) as f:
        for line in f:
            if line.startswith("N50"):
                stats["N50"] = int(line.strip().split('\t')[1])
            elif line.startswith("# contigs"):
                stats["num_contigs"] = int(line.strip().split('\t')[1])
            elif line.startswith("Total length"):
                stats["total_length"] = int(line.strip().split('\t')[1])
    return stats

def parse_circular_contigs(circularity_report_path, genome_size):
    """Return total circular contigs and if a likely circular chromosome is present."""
    CHROMOSOME_SIZE_THRESHOLD = 0.95 * float(genome_size)
    df = pd.read_csv(circularity_report_path, sep='\t')
    if "circular" not in df.columns:
        return 0, False

    circular_df = df[df["circular"].str.lower() == "yes"]
    num_circular = len(circular_df)


    # Check if any circular contig is large enough to be the chromosome
    has_circular_chromosome = (circular_df["Length"] >= CHROMOSOME_SIZE_THRESHOLD).any()

    return num_circular, has_circular_chromosome

def evaluate_assemblies(quast_dirs, genome_size):
    rankings = []

    for quast_dir in quast_dirs:
        report_path = os.path.join(quast_dir, "report.tsv")
        circularity_report_path = os.path.join(quast_dir, "circularity.tsv")

        if not os.path.exists(report_path) or not os.path.exists(circularity_report_path):
            continue

        stats = parse_quast_report(report_path)
        stats["num_circular"], stats["has_circular_chromosome"] = parse_circular_contigs(circularity_report_path, genome_size)
        stats["assembly_name"] = os.path.basename(quast_dir)
        rankings.append(stats)

    df = pd.DataFrame(rankings)
    if df.empty:
        print("No valid QUAST outputs found.")
        return None

    # Scoring formula (customized):
    # High bonus for circular chromosome
    # Moderate bonus for plasmid-like circular contigs
    # Heavy penalty for fragmented assemblies
    # Small bonus for N50
    df["score"] = (
        df["has_circular_chromosome"].astype(int) * 10 +  # very high priority
        df["num_circular"] * 2 -                          # minor bonus for other circular contigs
        df["num_contigs"] * 2 +                           # strong penalty for fragmented assemblies
        df["N50"] / 1e6                                   # small bonus for continuity
    )

    df_sorted = df.sort_values( by=["score", "N50", "total_length"], ascending=[False, False, False])

    #df_sorted = df.sort_values(by=["score", "total_length"], ascending=[False, False])

    print("\n📊 Assembly Ranking Summary:\n")
    # print(df_sorted[["assembly_name", "has_circular_chromosome", "num_circular", "num_contigs", "N50", "total_length", "score"]])
    print(df_sorted[["assembly_name", "has_circular_chromosome", "num_circular", "num_contigs", "N50", "total_length"]])

    best = df_sorted.iloc[0]
    print(f"\n✅ Best Assembly: {best['assembly_name']}")
    return best["assembly_name"]


min_overlap = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
basepath = "/wibicomfs/LTBS/HIPS_genome_assembly_miscellaneous/anna_nanopore/barcode10/flye/"
output_path = "/wibicomfs/LTBS/HIPS_genome_assembly_miscellaneous/anna_nanopore/quast_outputs/barcode10/"

for overlap in min_overlap:
    os.makedirs(output_path + str(overlap))
    os.system("quast -o " + output_path + str(overlap) + " " + basepath + str(overlap) + "/assembly.fasta")
    generate_contig_circularity_info(output_path + str(overlap), basepath)

evaluate_assemblies(glob.glob("/wibicomfs/LTBS/HIPS_genome_assembly_miscellaneous/anna_nanopore/quast_outputs/barcode10/*"), 5000000)