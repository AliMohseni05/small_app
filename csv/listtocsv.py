import csv

def list_name():
    names = []
    count_n=0
    while True:
        add_name = input(f"you have {count_n} name in this list Add name or type 'fin' to finish: ")
        if add_name.lower() == 'fin':
            break
        names.append(add_name)
        count_n +=1
    return names

def list_age():
    ages = []
    count_g=0
    while True:
        add_age = input(f"you have {count_g} age in this list Add age or type 'fin' to finish: ")
        if add_age.lower() == 'fin':
            break
        ages.append(add_age)
        count_g +=1

    return ages

# Collect the names and ages
names = list_name()
ages = list_age()

# Check if both lists are of the same length
if len(names) != len(ages):
    print("The number of names and ages must match. Please try again.")
else:
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