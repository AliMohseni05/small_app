import csv

# Sample lists
add_name=input("add name: ")
add_age=input("add age ")
names = ["Alice", "Bob", "Charlie", "David"]
ages = [25, 30, 35, 40]

# Specify the CSV file name
csv_file = "output_lists.csv"

# Open the CSV file for writing
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)

    # Write the header row
    writer.writerow(["Name", "Age"])

    # Write the rows using zip to combine the two lists
    writer.writerows(zip(names, ages))

print(f"Lists exported to {csv_file} successfully!")