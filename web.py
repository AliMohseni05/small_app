from selenium import webdriver
import keyboard  # Requires 'keyboard' module: pip install keyboard

def close_app():
    """Closes the browser and exits the script."""
    driver.quit()
    exit()

# Initialize the WebDriver (modify if using a different browser)
driver = webdriver.Chrome()

# Open the desired website
driver.get("https://realpython.com/")  # Replace with your target URL

# Register the Ctrl+Q hotkey to close the app
keyboard.add_hotkey('ctrl+q', close_app)

print("Browser opened. Press Ctrl+Q and Ctrl+C to close the browser and exit.")

try:
    # Keep the script running to listen for the hotkey
    keyboard.wait()
except KeyboardInterrupt:
    # Handle Ctrl+C to ensure browser closes
    print("\nClosing browser...")
    driver.quit()