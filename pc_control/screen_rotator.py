import sys
import platform
import keyboard
import time
import pyautogui

def rotate_windows():
    try:
        import ctypes
        from ctypes import wintypes

        while True: 
            #current_time = time.time()  # Get the current time
            #elapsed_time = current_time - start_time  # Calculate elapsed time
            
            # Click at the current mouse position
            pyautogui.click()

            # Check if 'Q' is pressed
            if keyboard.is_pressed('q'):
                print("Stopping the script.")
                break
            
            # Stop the script after the specified duration
           # if elapsed_time >= duration:
                #print("{duration} seconds have passed. Stopping the script.")
                #break
            
            #time.sleep(0.002)

    except ImportError:
        print("This script requires ctypes to run on Windows.")
        return 1

    # Constants
    DM_DISPLAYORIENTATION = 0x00000080
    DMDO_90 = 1
    ENUM_CURRENT_SETTINGS = -1
    CDS_UPDATEREGISTRY = 0x01

    class DEVMODEW(ctypes.Structure):
        _fields_ = [
            # Corrected WCHAR references:
            ("dmDeviceName", ctypes.c_wchar * 32),
            ("dmSpecVersion", wintypes.WORD),
            ("dmDriverVersion", wintypes.WORD),
            ("dmSize", wintypes.WORD),
            ("dmDriverExtra", wintypes.WORD),
            ("dmFields", wintypes.DWORD),
            ("dmOrientation", ctypes.c_short),
            ("dmPaperSize", ctypes.c_short),
            ("dmPaperLength", ctypes.c_short),
            ("dmPaperWidth", ctypes.c_short),
            ("dmScale", ctypes.c_short),
            ("dmCopies", ctypes.c_short),
            ("dmDefaultSource", ctypes.c_short),
            ("dmPrintQuality", ctypes.c_short),
            ("dmColor", ctypes.c_short),
            ("dmDuplex", ctypes.c_short),
            ("dmYResolution", ctypes.c_short),
            ("dmTTOption", ctypes.c_short),
            ("dmCollate", ctypes.c_short),
            ("dmFormName", ctypes.c_wchar * 32),  # Fixed here too
            ("dmLogPixels", wintypes.WORD),
            ("dmBitsPerPel", wintypes.DWORD),
            ("dmPelsWidth", wintypes.DWORD),
            ("dmPelsHeight", wintypes.DWORD),
            ("dmDisplayFlags", wintypes.DWORD),
            ("dmDisplayFrequency", wintypes.DWORD),
            ("dmDisplayOrientation", wintypes.DWORD),
        ]
    devmode = DEVMODEW()
    devmode.dmSize = ctypes.sizeof(DEVMODEW)
    
    if not ctypes.windll.user32.EnumDisplaySettingsW(None, ENUM_CURRENT_SETTINGS, ctypes.byref(devmode)):
        print("Failed to get display settings.")
        return 1

    original_orientation = devmode.dmDisplayOrientation

    if original_orientation == DMDO_90:
        print("Screen is already rotated 90 degrees.")
        return 0

    devmode.dmDisplayOrientation = DMDO_90
    devmode.dmFields |= DM_DISPLAYORIENTATION

    result = ctypes.windll.user32.ChangeDisplaySettingsExW(
        None,
        ctypes.byref(devmode),
        None,
        CDS_UPDATEREGISTRY,
        None
    )

    if result == 0:
        print("Screen rotated 90 degrees successfully.")
        return 0
    else:
        print(f"Failed to rotate screen. Error code: {result}")
        return 1
    
def rotate_linux():

    try:
        import subprocess
    except ImportError:
        print("This script requires subprocess to run on Linux.")
        return 1

    try:
        while True: 
            current_time = time.time()  # Get the current time
            elapsed_time = current_time - start_time  # Calculate elapsed time
            
            # Click at the current mouse position
            pyautogui.click()

            # Get connected display name
            output = subprocess.check_output(['xrandr']).decode()
            displays = [line.split()[0] for line in output.splitlines() if ' connected' in line]
        
            if not displays:
                print("No connected displays found.")
                return 1
            
            display = displays[0]
        
            # Rotate display
            subprocess.run(
                ['xrandr', '--output', display, '--rotate', 'left', '--auto'],
                check=True
            )
            print(f"Rotated {display} 90 degrees left.")
            # Check if 'Q' is pressed
            if keyboard.is_pressed('q'):
                print("Stopping the script.")
                break
            
            # Stop the script after the specified duration
            if elapsed_time >= duration:
                print("{duration} seconds have passed. Stopping the script.")
                break
            
            time.sleep(0.002)
        
            return 0
    except subprocess.CalledProcessError as e:
        print(f"Error rotating screen: {e}")
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}")
        return 1

if __name__ == "__main__":
    system = platform.system()
    
    if system == "Windows":
        sys.exit(rotate_windows())
    elif system == "Linux":
        sys.exit(rotate_linux())
    else:
        print(f"Unsupported operating system: {system}")
        sys.exit(1)