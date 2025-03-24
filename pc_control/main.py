import sys
import ctypes
import ctypes.wintypes
import keyboard

def set_orientation(orientation):
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

    devmode = DEVMODEW()
    devmode.dmSize = ctypes.sizeof(DEVMODEW)
    
    if not ctypes.windll.user32.EnumDisplaySettingsW(None, -1, ctypes.byref(devmode)):
        return False

    original = devmode.dmDisplayOrientation
    if devmode.dmDisplayOrientation == orientation:
        return original  # Already in requested orientation
    
    devmode.dmDisplayOrientation = orientation
    devmode.dmFields |= 0x00000080  # DM_DISPLAYORIENTATION
    
    result = ctypes.windll.user32.ChangeDisplaySettingsExW(
        None, ctypes.byref(devmode), None, 0x01, None
    )
    return original if result == 0 else None

def is_admin():
    try:
        return ctypes.windll.shell32.IsUserAnAdmin()
    except:
        return False

def main():
    if not is_admin():
        print("Please run as Administrator")
        sys.exit(1)

    try:
        original_orientation = set_orientation(1)  # 1 = DMDO_90
        if original_orientation is None:
            print("Rotation failed")
            sys.exit(1)
            
        print("Screen rotated 90°. Press Ctrl+Q to restore and exit")
        keyboard.add_hotkey('ctrl+q', lambda: set_orientation(original_orientation))
        keyboard.wait()  # Keep script running
        
    except KeyboardInterrupt:
        set_orientation(original_orientation)
    except Exception as e:
        print(f"Error: {str(e)}")
    finally:
        keyboard.unhook_all()

if __name__ == "__main__":
    main()
    sys.exit(0)