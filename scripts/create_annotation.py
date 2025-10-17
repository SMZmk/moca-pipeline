import os
import argparse
import pyfaidx


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create annotation file from FASTA.")
    parser.add_argument("--fasta_path", "-f", type=str, help="Path to the input FASTA file.",)
    return parser.parse_args()


def create_annotation_file(fasta_path: str) -> None:
    output_path = os.path.splitext(fasta_path)[0] + ".gff3"
    print(f"Creating annotation file at {output_path} from {fasta_path}...")
    fasta = pyfaidx.Fasta(fasta_path)
    with open(output_path, "w") as out_file:
        out_file.write("##gff-version 3\n")
        for seq_id, sequence in fasta.items():
            region = seq_id
            source = "DREAM"
            type_ = "region"
            start = 1001
            end = len(sequence) - 1000
            score = "."
            strand = "+"
            phase = "."
            attributes = f"ID=chr:{seq_id}"
            line = "\t".join([region, source, type_, str(start), str(end), score, strand, phase, attributes])
            out_file.write(f"{line}\n")
            type_ = "gene"
            attributes = f"ID=gene:{seq_id}"
            line = "\t".join([region, source, type_, str(start), str(end), score, strand, phase, attributes])
            out_file.write(f"{line}\n")
    print(f"Annotation file created at {output_path}")


def check_duplicate_keys(fasta_path: str) -> str:
    with open(fasta_path, "r") as fasta_file:
        lines = fasta_file.readlines()
    keys = {}
    deduped = {}
    duplicates = False
    for i, line in enumerate(lines):
        if line.startswith(">"):
            key = line[1:].strip()
            keys[key] = keys.get(key, []) + [i]
    for key, occurrences in keys.items():
        different_sequences = set([0])
        if len(occurrences) > 1:
            print(f"Found duplicate key: {key} at lines {occurrences}")
            duplicates = True
            sequences = []
            for occurrence in occurrences:
                line = lines[occurrence + 1]
                sequences.append(line)
            for i, seq in enumerate(sequences):
                for j, (letter_a, letter_b) in enumerate(zip(sequences[0], seq)):
                    if letter_a != letter_b:
                        different_sequences.add(i)
                        print(f"sequences are not identical at position {j}.")
            for i, difference_index in enumerate(different_sequences):
                deduped[f"{key}_{i}"] = [occurrences[difference_index]]
        else:
            deduped[key] = occurrences
        
    if not duplicates:
        return fasta_path

    deduplicate_path = os.path.splitext(fasta_path)[0] + "_deduplicated.fa"
    print(f"Found duplicate keys in {fasta_path}. Deduplicating...")
    with open(deduplicate_path, "w") as dedup_file:
        for key, occurrences in deduped.items():
            dedup_file.write(f">{key}\n")
            dedup_file.write(lines[occurrences[0] + 1])
    print(f"Deduplicated FASTA file saved to {deduplicate_path}")
    return deduplicate_path



def main() -> None:
    args = parse_args()
    deduplicated_path = check_duplicate_keys(args.fasta_path)
    create_annotation_file(deduplicated_path)

if __name__ == "__main__":
    main()
    # moca_blue/ref_seq/orig_genomes/selected_genes_simon_min_msr_250711_110304_129027_pareto_fronts.fa