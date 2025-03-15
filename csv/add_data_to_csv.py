import csv

# Sample list to add
new_data = ["Eve", 28]

# Specify the CSV file name (make sure this file already exists)
csv_file = "output_lists.csv"

# Open the CSV file for appending
with open(csv_file, mode='a', newline='') as file:
    writer = csv.writer(file)

    # Write the new data as a new row
    writer.writerow(new_data)

print(f"New data added to {csv_file} successfully!")