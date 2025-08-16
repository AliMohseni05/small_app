import os
import ftplib
from concurrent.futures import ThreadPoolExecutor
from threading import Lock

def download_plant_proteins():
    os.makedirs("plant", exist_ok=True)

    print("Connecting to NCBI FTP server...")
    with ftplib.FTP("ftp.ncbi.nlm.nih.gov") as ftp:
        ftp.login()
        ftp.cwd("refseq/release/plant")

        # Get list of protein files
        files = []
        ftp.retrlines('NLST', files.append)
        plant_files = [f for f in files if f.startswith("plant.") and f.endswith(".protein.faa.gz")]

        if not plant_files:
            raise Exception("No plant protein files found on server")

        print(f"Found {len(plant_files)} plant protein files")

        ftp.voidcmd('NOOP')  # keep connection alive

        # Lock for thread-safe print statements
        print_lock = Lock()

        def download_file(filename):
            with ftplib.FTP("ftp.ncbi.nlm.nih.gov") as ftp_thread:
                ftp_thread.login()
                ftp_thread.cwd("refseq/release/plant")

                local_path = os.path.join("plant", filename)
                total_size = ftp_thread.size(filename)
                downloaded = 0

                def callback(data):
                    nonlocal downloaded
                    downloaded += len(data)
                    f.write(data)
                    percent = downloaded * 100 // total_size
                    with print_lock:
                        print(f"\rDownloading {filename}: {percent}% complete", end='')

                with open(local_path, 'wb') as f:
                    ftp_thread.retrbinary(f"RETR {filename}", callback)

                with print_lock:
                    print(f"\nFinished downloading {filename}")

        print("Starting downloads...")
        with ThreadPoolExecutor(max_workers=5) as executor:
            executor.map(download_file, plant_files)

if __name__ == "__main__":
    download_plant_proteins()
