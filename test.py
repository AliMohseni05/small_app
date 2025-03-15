import csv

# Sample dictionary
data = {
    "Name": ["Alice", "Bob", "Charlie"],
    "Age": [25, 30, 35],
    "City": ["New York", "Los Angeles", "Chicago"]
}

# Specify the CSV file name
csv_file = "output.csv"

# Open the CSV file for writing
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)

    # Write the header (keys of the dictionary)
    writer.writerow(data.keys())

    # Write the values (zip combines the lists together)
    writer.writerows(zip(*data.values()))

print(f"Data exported to {csv_file} successfully!")