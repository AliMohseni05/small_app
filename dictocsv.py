import csv

def create_dictionary():
    my_dict = {}
    
    while True:
        key = input("Enter the key (or type 'exit' to finish): ")
        if key.lower() == 'exit':
            break
        
        value = input(f"Enter the value for '{key}': ")
        
        my_dict[key] = value
        print(f"Added: {key}: {value}")
    
    return my_dict
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


if __name__ == "__main__":
    dictocsv = create_dictionary()
    print("\nYour final dictionary:")
    print(my_dict)


