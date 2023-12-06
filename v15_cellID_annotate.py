import argparse
import gzip
import pandas as pd

parser = argparse.ArgumentParser(description='Annotate cell types')
parser.add_argument('--input-UMI-fastq', required=True, help='input h5ad file path')
parser.add_argument('--input-perf-reads', help='gene class file in CSV format', required=True)
parser.add_argument('--output', required=True, help='output file path')
args = parser.parse_args()

# Function to extract read ID from the perf_reads_file or umi_fastq (R1) file
def extract_read_id_perf(line):
    return line.split("_")[0][1:]

def extract_read_id_umi(line):
    return line.split(" ")[0][1:]

# Paths to input files
perf_reads_file = args.input_perf_reads
umi_barcode_file = args.input_UMI_fastq

# Load the CSV file with UMI and cell type information
umi_cell_info = pd.read_csv("umi_cell_info.csv")

# Create a dictionary to store UMIs and cell types from the CSV file
umi_cell_dict = dict(zip(umi_cell_info["UMI"], umi_cell_info["Cell_Type"]))
# Create a list to store the results
results = []

# Read the perf_reads_file and process each line
with open(perf_reads_file, "r") as file:
    for line in file:
        read_id = extract_read_id_perf(line)
        #print("Perf_Read # ", read_id)
        with gzip.open(umi_barcode_file, "rt") as r1_file:
            while True:
                read_id_line = r1_file.readline().strip()
                if not read_id_line:
                    break
                if read_id_line.startswith('@'):
                    read_id_l = extract_read_id_umi(read_id_line)
                    if read_id_l != read_id:
                        continue
                    print("Matching read ID found: ", read_id_l)
                    umi = r1_file.readline().strip()[:16]
                    #barcode = r1_file.readline().strip()[-10:]
                    #print("Matching UMI found: ", umi)

                    if not umi:
                        #print(f"Warning: UMI not found for read ID {read_id}.")
                        continue
                    # Look up the cell type for the extracted UMI from the CSV file
                    cell_type = umi_cell_dict.get(umi, "Unknown")
                    results.append({"Read_ID": read_id, "Cell_Type": cell_type, "UMI": umi})
                    
                                #print(results)
                    # Flush the output to ensure prints are displayed promptly
                    import sys
                    sys.stdout.flush()

                    break  # Exit the inner loop after finding the UMI for the read ID


# Convert the list of results to a DataFrame
result_df = pd.DataFrame(results)

# Add to CSV
output_file = args.output
result_df.to_csv(output_file, sep="\t", index=False)
