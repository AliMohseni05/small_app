import pyautogui
import keyboard
import time
import tkinter as tk
from tkinter import messagebox

# Function to start the clicking action
def start_clicking():
    duration = 10
    start_time = time.time()  # Record the start time
    print("Press 'Q' to stop the script early.")
    print(f"The script will run for {duration} seconds.")
    
    try:
        while True:
            current_time = time.time()  # Get the current time
            elapsed_time = current_time - start_time  # Calculate elapsed time
            
            # Click at the current mouse position
            pyautogui.click()
            
            # Check if 'Q' is pressed
            if keyboard.is_pressed('q'):
                print("Stopping the script.")
                break
            
            # Stop the script after the specified duration
            if elapsed_time >= duration:
                print("{duration} seconds have passed. Stopping the script.")
                break
            
            time.sleep(0.02)  # Adjust to control the speed of clicking

    except KeyboardInterrupt:
        print("Script interrupted manually.")

# Function to create a popup window
def show_popup():
    # Create a new tkinter window
    popup = tk.Tk()
    popup.title("Clicker App")
    popup.geometry("300x150")

    # Informational message
    messagebox.showinfo("Welcome", "Click 'OK' to start the clicking process.")

    # Start the clicking process after the popup is closed
    start_clicking()
    popup.destroy()  # Close the popup window

# Run the popup window
show_popup()