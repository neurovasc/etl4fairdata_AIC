import csv
import random

# File paths
input_csv = "fake-GAIA-extraction.csv"  # Path to the input CSV file
output_csv = "fake-GAIA-extraction-clean.csv"  # Path to save the output CSV file

# Read the CSV file with UTF-8 encoding
with open(input_csv, 'r', newline='', encoding='utf-8') as infile:
    reader = list(csv.reader(infile, delimiter=',', quotechar='"'))
    
    # Separate header and data
    header = reader[0]  # First row is the header
    data = reader[1:]   # Remaining rows are the data
    
    # Transpose the data (columns to rows for easy shuffling)
    columns = list(zip(*data))
    
    # Shuffle each column individually
    shuffled_columns = [list(col) for col in columns]
    for col in shuffled_columns:
        random.shuffle(col)
    
    # Transpose the columns back to rows
    shuffled_data = list(zip(*shuffled_columns))

# Write the shuffled data to a new CSV file with UTF-8 encoding
with open(output_csv, 'w', newline='', encoding='utf-8') as outfile:
    writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    # Write header
    writer.writerow(header)
    
    # Write shuffled data
    writer.writerows(shuffled_data)

print(f"Data shuffled and saved to {output_csv}")
