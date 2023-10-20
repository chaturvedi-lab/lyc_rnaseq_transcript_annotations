import csv
import argparse

# Load IPR entries from a file into a dictionary
def load_entries(file_path):
    entries = {}
    with open(file_path, "r") as entry_file:
        for line in entry_file:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                entry_key, entry_value = parts[0], parts[1]
                entries[entry_key] = entry_value
    return entries

# Load GO annotations from a file into a dictionary
def load_go_annotations(go_file):
    go_annotations = {}
    current_go_id = None
    current_go_name = None
    current_go_namespace = None
    is_reading_entry = False

    with open(go_file, "r") as go_csv:
        for line in go_csv:
            line = line.strip()

            if line == "[Term]":
                is_reading_entry = True
                current_go_id = None
                current_go_name = None
                current_go_namespace = None
            elif is_reading_entry:
                if line.startswith("id: "):
                    current_go_id = line[4:]
                elif line.startswith("name: "):
                    current_go_name = line[6:]
                elif line.startswith("namespace: "):
                    current_go_namespace = line[11:]
                elif line.startswith("is_a: "):
                    go_annotations[current_go_id] = {"name": current_go_name, "namespace": current_go_namespace}
                    current_go_id = None
                    current_go_name = None
                    current_go_namespace = None

    return go_annotations

# Function to split and process the data
def split_and_process_data(input_file, output_file):
    ipr_file = "entry.list"  # Hardcoded IPR entry file
    go_file = "go-basic.obo"  # Hardcoded GO annotation file
    ipr_entries = load_entries(ipr_file)
    go_annotations = load_go_annotations(go_file)
    #print(go_annotations)

    with open(input_file, "r") as input_csv, open(output_file, "w", newline="") as output_csv:
        csv_reader = csv.reader(input_csv, delimiter='\t')
        csv_writer = csv.writer(output_csv, delimiter='\t')

        for row in csv_reader:
            #print(row)
            # Split "Column 1" into multiple columns
            column1_parts = row[0].split(':')
            orf_parts = column1_parts[1].split('_')
            scaffold = column1_parts[0]
            #print(orf_parts)
            start_position, stop_position, orf_value = orf_parts[0].split('-')[0], orf_parts[0].split('-')[1], orf_parts[1]
            
            # Extract IPR and GO numbers from the respective columns
            ipr_numbers = row[11].split('|') if row[11] else []  # Split IPR numbers if present
            print(ipr_numbers)
            if len(row) >= 14:  # Ensure the row contains at least 14 columns
                go_numbers = row[13].split('|') if row[13] else []  # Split GO numbers if present
            else:
                go_numbers = []
            #print(go_numbers)

            # Create a list with the selected columns
            selected_columns = [scaffold, start_position, stop_position, orf_value, row[1], "|".join(ipr_numbers), "|".join(go_numbers)]

            # Check if IPR number exists in the dictionary and append the corresponding data
            ipr_data = row[12]
            go_data = []  # Define empty lists for GO data
            go_name_data = []  # Define empty lists for GO term names
            go_namespace_data = []  # Define empty lists for GO namespaces

            for ipr_number in ipr_numbers:
                 if ipr_number in ipr_entries:
                     selected_columns.append(ipr_entries[ipr_number])
                 else:
                     selected_columns.append("")  # Empty value if no match

            for go_number in go_numbers:
                #print(go_number)
                if go_number in go_annotations:
                    go_name_data.append(go_annotations[go_number]["name"])
                    go_namespace_data.append(go_annotations[go_number]["namespace"])

            # Create a list with the selected columns
            selected_columns = [scaffold, start_position, stop_position, orf_value, ipr_number, ipr_data,go_number,"|".join(go_data), "|".join(go_name_data), "|".join(go_namespace_data)]

            # Write the extracted data to the output file
            csv_writer.writerow(selected_columns)

def main():
    # Create a command-line argument parser
    parser = argparse.ArgumentParser(description="Process input data and create a tab-separated output file")
    parser.add_argument("input_file", help="Path to the input file")
    parser.add_argument("output_file", help="Path to the output file")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the function to process the data
    split_and_process_data(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
