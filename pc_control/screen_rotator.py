import sys
import ctypes
import ctypes.wintypes
import keyboard

# Constants
DM_DISPLAYORIENTATION = 0x00000080
DMDO_DEFAULT = 0
DMDO_90 = 1
DMDO_180 = 2
DMDO_270 = 3

class DEVMODEW(ctypes.Structure):
    _fields_ = [
        ("dmDeviceName", ctypes.c_wchar * 32),
        ("dmSpecVersion", ctypes.wintypes.WORD),
        ("dmDriverVersion", ctypes.wintypes.WORD),
        ("dmSize", ctypes.wintypes.WORD),
        ("dmDriverExtra", ctypes.wintypes.WORD),
        ("dmFields", ctypes.wintypes.DWORD),
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
        ("dmFormName", ctypes.c_wchar * 32),
        ("dmLogPixels", ctypes.wintypes.WORD),
        ("dmBitsPerPel", ctypes.wintypes.DWORD),
        ("dmPelsWidth", ctypes.wintypes.DWORD),
        ("dmPelsHeight", ctypes.wintypes.DWORD),
        ("dmDisplayFlags", ctypes.wintypes.DWORD),
        ("dmDisplayFrequency", ctypes.wintypes.DWORD),
        ("dmDisplayOrientation", ctypes.wintypes.DWORD),
    ]

def set_display_orientation(orientation):
    devmode = DEVMODEW()
    devmode.dmSize = ctypes.sizeof(DEVMODEW)
    
    if not ctypes.windll.user32.EnumDisplaySettingsW(None, -1, ctypes.byref(devmode)):
        return False

    original = devmode.dmDisplayOrientation
    if devmode.dmDisplayOrientation == orientation:
        return original
    
    devmode.dmDisplayOrientation = orientation
    devmode.dmFields |= DM_DISPLAYORIENTATION
    
    result = ctypes.windll.user32.ChangeDisplaySettingsExW(
        None, ctypes.byref(devmode), None, 0x01, None)
    return original if result == 0 else None

def is_admin():
    try:
        return ctypes.windll.shell32.IsUserAnAdmin()
    except:
        return False

def main():
    if not is_admin():
        ctypes.windll.shell32.ShellExecuteW(
            None, "runas", sys.executable, __file__, None, 1)
        sys.exit(0)

    try:
        original = set_display_orientation(DMDO_90)
        if original is None:
            print("Failed to rotate display")
            sys.exit(1)
            
        print(f"Display rotated 90° (Original: {original}). Press Ctrl+Q to restore")
        keyboard.add_hotkey('ctrl+q', lambda: [
            set_display_orientation(original),
            sys.exit(0)
        ])
        
        keyboard.wait()
        
    except KeyboardInterrupt:
        set_display_orientation(original)
        sys.exit(0)
    finally:
        keyboard.unhook_all()

if __name__ == "__main__":
    main()