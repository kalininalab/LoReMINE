import os
import glob
import pandas as pd

# Set your expected genome size here (in bp)
EXPECTED_GENOME_SIZE = 5000000  # e.g., ~5 Mb for many bacteria
CHROMOSOME_SIZE_THRESHOLD = 0.95 * EXPECTED_GENOME_SIZE  # 95% of genome size

def parse_quast_report(report_path):

    """Extract N50, #contigs, and total length from QUAST report.tsv"""

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

    """Return total circular contigs and if a likely circular chromosome is present"""

    CHROMOSOME_SIZE_THRESHOLD = 0.95 * float(genome_size)
    df = pd.read_csv(circularity_report_path, sep='\t')
    if "circular" not in df.columns:
        return 0, False

    circular_df = df[df["circular"].str.lower() == "yes"]
    num_circular = len(circular_df)


    # Check if any circular contig is large enough to be the chromosome
    has_circular_chromosome = (circular_df["Length"] >= CHROMOSOME_SIZE_THRESHOLD).any()

    return num_circular, has_circular_chromosome

def evaluate_assemblies(quast_dirs, genome_size, best_assembly_store_path):

    """Evaluate multiple assemblies and pick the best assembly"""

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

    # The scoring function deciding the score for each assembly based on circularity of chromosome, no. of circular contigs, total no. of contigs and size of chromosome (bigger the better)

    df["score"] = (
        df["has_circular_chromosome"].astype(int) * 10 +  # very high priority
        df["num_circular"] * 2 -                          # minor bonus for other circular contigs
        df["num_contigs"] * 2 +                           # strong penalty for fragmented assemblies
        df["N50"] / 1e6                                   # small bonus for continuity
    )

    df_sorted = df.sort_values( by=["score", "N50", "total_length"], ascending=[False, False, False])

    # Storing the best assembly in a file
    best_assembly_file = open(best_assembly_store_path + "chosen_best_assembly.txt", 'x')

    #df_sorted = df.sort_values(by=["score", "total_length"], ascending=[False, False])

    print("\n📊 Assembly Ranking Summary:\n")
    best_assembly_file.write("\n📊 Assembly Ranking Summary:\n\n")
    # print(df_sorted[["assembly_name", "has_circular_chromosome", "num_circular", "num_contigs", "N50", "total_length", "score"]])
    print(df_sorted[["assembly_name", "has_circular_chromosome", "num_circular", "num_contigs", "N50", "total_length"]])
    best_assembly_file.write(str(df_sorted[["assembly_name", "has_circular_chromosome", "num_circular", "num_contigs", "N50", "total_length"]]) + '\n\n')

    best = df_sorted.iloc[0]
    print(f"\n✅ Best Assembly: {best['assembly_name']}")
    best_assembly_file.write(f"\n✅ Best Assembly according to us: {best['assembly_name']}")
    return best["assembly_name"]
