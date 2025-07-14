from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
import time

def run_cpc2_analysis(fasta_file, output_file):
    """
    Submits FASTA sequences to CPC2 website and saves results
    """
    # Read FASTA sequences
    with open(fasta_file, 'r') as f:
        fasta_data = f.read()
    
    # Configure Chrome options
    chrome_options = Options()
    chrome_options.add_argument("--disable-gpu")
    chrome_options.add_argument("--no-sandbox")
    chrome_options.add_argument("--disable-dev-shm-usage")
    chrome_options.add_argument("--disable-infobars")
    chrome_options.add_argument("--disable-extensions")
    chrome_options.add_argument("--disable-logging")
    chrome_options.add_argument("--log-level=3")
    
    # Initialize WebDriver
    driver = webdriver.Chrome(options=chrome_options)
    
    try:
        # Open CPC2 website
        print("Opening CPC2 website...")
        driver.get('https://cpc2.gao-lab.org/')
        
        # Switch to text input mode
        print("Switching to text input mode...")
        text_btn = WebDriverWait(driver, 20).until(
            EC.element_to_be_clickable((By.XPATH, "//button[contains(text(), 'Text')]"))
        )
        text_btn.click()
        
        # Enter FASTA sequences
        print("Entering FASTA sequences...")
        text_area = WebDriverWait(driver, 20).until(
            EC.visibility_of_element_located((By.ID, "sequence_text"))
        )
        text_area.clear()
        text_area.send_keys(fasta_data)
        
        # Submit the form
        print("Submitting sequences for analysis...")
        submit_btn = driver.find_element(By.XPATH, "//input[@type='submit' and @value='Run']")
        submit_btn.click()
        
        # Wait for results to load
        print("Waiting for results (this may take several minutes)...")
        WebDriverWait(driver, 300).until(
            EC.presence_of_element_located((By.ID, "result-table"))
        )
        
        # Get result table
        print("Results loaded, saving to file...")
        result_table = driver.find_element(By.ID, "result-table")
        
        # Save results to file
        with open(output_file, 'w') as f:
            f.write(result_table.text)
            
        print(f"Success! Results saved to {output_file}")
        
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        # Save page source for debugging
        with open('error_page.html', 'w', encoding='utf-8') as f:
            f.write(driver.page_source)
        print("Saved error page to 'error_page.html'")
        
    finally:
        driver.quit()

if __name__ == "__main__":
    input_file = "input.fasta"  # Replace with your input file
    output_file = "cpc2_results.txt"  # Replace with desired output file
    
    run_cpc2_analysis(input_file, output_file)