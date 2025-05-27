import pandas as pd
# coumt indevitual item in colome for this example traget 
file_path = 'C:/Users/ARASH-STORE/Documents/ALL.xlsx'  # Your Excel file path
output_txt = 'unique_value_counts2.txt'  # Output text file path

try:
    # Load Excel file
    df = pd.read_excel(file_path)
    
    # Count unique values in the "Target" column
    counts = df["Target"].value_counts()
    
    # Export counts to a text file
    with open(output_txt, 'w', encoding='utf-8') as f:
        for value, count in counts.items():
            f.write(f"{value}\t{count}\n")
    
    print(f"Unique value counts saved to '{output_txt}'")

except FileNotFoundError:
    print(f"Error: File '{file_path}' not found.")
except Exception as e:
    print(f"An error occurred: {str(e)}")

# Optional: print first few rows to check data
df = pd.read_excel(file_path)
print(df.head())
