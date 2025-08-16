import os
import gzip
import ftplib
import shutil
from concurrent.futures import ThreadPoolExecutor

def download_plant_proteins():
    # Create directory structure
    os.makedirs("plant", exist_ok=True)
    
    # Connect to NCBI FTP server
    print("Connecting to NCBI FTP server...")
    with ftplib.FTP("ftp.ncbi.nlm.nih.gov") as ftp:
        ftp.login()  # Anonymous login
        ftp.cwd("refseq/release/plant")
        
        # Get list of plant protein files
        files = []
        ftp.retrlines('NLST', files.append)
        plant_files = [f for f in files if f.startswith("plant.") 
                      and f.endswith(".protein.faa.gz")]
        
        if not plant_files:
            raise Exception("No plant protein files found on server")
        
        print(f"Found {len(plant_files)} plant protein files")

        # Download files in parallel
        def download_file(filename):
            print(f"Downloading {filename}...")
            local_path = os.path.join("plant", filename)
            with open(local_path, 'wb') as f:
                ftp.retrbinary(f"RETR {filename}", f.write)
            return filename

        print("Starting downloads...")
        with ThreadPoolExecutor(max_workers=5) as executor:
            executor.map(download_file, plant_files)

    # Process files
    #decompress_files()
    #combine_faa_files()
    
    print("Operation completed successfully!")
    final_file = "all_refseq_plant_proteins.faa"
    print(f"Combined file created: {final_file}")
    print(f"Size: {os.path.getsize(final_file)/1024/1024:.2f} MB")

def decompress_files():
    print("Decompressing files...")
    for filename in os.listdir("plant"):
        if filename.endswith(".gz"):
            gz_path = os.path.join("plant", filename)
            faa_path = os.path.join("plant", filename[:-3])  # Remove .gz
            
            with gzip.open(gz_path, 'rb') as f_in:
                with open(faa_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            os.remove(gz_path)

def combine_faa_files():
    print("Combining FASTA files...")
    with open("all_refseq_plant_proteins.faa", 'wb') as outfile:
        for filename in os.listdir("plant"):
            if filename.endswith(".faa"):
                with open(os.path.join("plant", filename), 'rb') as infile:
                    shutil.copyfileobj(infile, outfile)

if __name__ == "__main__":
    download_plant_proteins()