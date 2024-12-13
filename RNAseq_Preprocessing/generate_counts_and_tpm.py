import os
import pandas as pd

# Directory containing the counts files
counts_dir = "/home/cbb575_cl2658/palmer_scratch/final_project/processed_data/counts/"  # Update this to your directory
output_file_counts = "counts_with_chr_and_length.csv"
output_file_tpm = "tpm_matrix.csv"

# Initialize dictionary for storing sample data and a list for annotations
data = {}
annotations_list = []  # To store annotations from all files

# Traverse all *_counts.txt files
for file in os.listdir(counts_dir):
    if file.endswith("_counts.txt") and not file.endswith(".summary"):
        file_path = os.path.join(counts_dir, file)
        
        # Extract sample name (removing "_counts.txt")
        sample_name = file.replace("_counts.txt", "")
        
        # Read the counts file (skip header comments starting with "#")
        df = pd.read_csv(file_path, sep="\t", comment="#")
        
        # Set Geneid as the index and extract the counts column
        data[sample_name] = df.set_index("Geneid")["bam_data/" + sample_name + "_sorted.bam"]
        
        # Extract annotations (Chr and Length) and append to the list
        annotations_list.append(df.set_index("Geneid")[["Chr", "Length"]])

# Verify annotations are consistent across files
annotations = annotations_list[0]  # Start with annotations from the first file
for i, annotation in enumerate(annotations_list[1:], start=1):
    if not annotations.equals(annotation):
        raise ValueError(f"Inconsistent annotations in file {i+1}. Ensure all files have the same Chr and Length.")

# Combine all samples into a single DataFrame
counts_matrix = pd.DataFrame(data)

# Merge annotations (Chr and Length) into the counts matrix
counts_with_annotations = annotations.join(counts_matrix)

# Save the counts matrix with annotations to a CSV file
counts_with_annotations.to_csv(output_file_counts)
print(f"Counts matrix with Chr and Length saved to: {output_file_counts}")

# Calculate TPM
gene_lengths_kb = counts_with_annotations["Length"] / 1000  # Convert Length to kilobases
counts = counts_with_annotations.drop(columns=["Chr", "Length"])  # Drop Chr and Length for TPM calculation
rpk = counts.div(gene_lengths_kb, axis=0)  # Calculate RPK
rpk_sum = rpk.sum(axis=0)  # Sum RPKs per sample
tpm = rpk.div(rpk_sum, axis=1) * 1e6  # Calculate TPM

# Save TPM matrix to a CSV file
tpm.to_csv(output_file_tpm)
print(f"TPM matrix saved to: {output_file_tpm}")

