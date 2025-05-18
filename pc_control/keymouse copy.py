from pynput import keyboard, mouse
import time

# Initialize mouse controller
mouse_controller = mouse.Controller()

# Settings
move_speed = 12         # Base movement speed (pixels per key press)
fast_move_multiplier = 7 # Hold Shift to move faster

# Track pressed keys and states
pressed_keys = set()
left_click_held = False
right_click_held = False

def on_press(key):
    """Handle presses"""
    try:
        pressed_keys.add(key.char.lower())
    except AttributeError:
        pressed_keys.add(key)

def on_release(key):
    """Handle key releases"""
    try:
        pressed_keys.discard(key.char.lower())
    except AttributeError:
        pressed_keys.discard(key)
    
    # Exit on Escape key
    if key == keyboard.Key.esc:
        return False

def process_movement():
    """Update mouse position based on pressed keys"""
    global mouse_controller, move_speed
    
    # Calculate movement speed
    speed = move_speed
    if keyboard.Key.shift in pressed_keys:
        speed *= fast_move_multiplier
    
    # Calculate movement delta
    dx, dy = 0, 0
    if 'w' in pressed_keys: dx -= speed #a
    if 's' in pressed_keys: dx += speed #d
    if 'd' in pressed_keys: dy -= speed #w
    if 'a' in pressed_keys: dy += speed #s
    
    # Move mouse
    if dx != 0 or dy != 0:
        mouse_controller.move(dx, dy)

def process_clicks():
    """Handle mouse clicks"""
    global left_click_held, right_click_held
    
    # Left click (L key)
    if 'e' in pressed_keys and not left_click_held:
        mouse_controller.press(mouse.Button.left)
        left_click_held = True
    elif 'e' not in pressed_keys and left_click_held:
        mouse_controller.release(mouse.Button.left)
        left_click_held = False
    
    # Right click (R key)
    if 'r' in pressed_keys and not right_click_held:
        mouse_controller.press(mouse.Button.right)
        right_click_held = True
    elif 'r' not in pressed_keys and right_click_held:
        mouse_controller.release(mouse.Button.right)
        right_click_held = False

# Start keyboard listener
with keyboard.Listener(
    on_press=on_press,
    on_release=on_release) as listener:

    print("Mouse Keyboard Controller Active!")
    print("WASD: Move | L: Left Click | R: Right Click")
    print("Shift: Speed Boost | Esc: Exit")

    while listener.running:
        process_movement()
        process_clicks()
        time.sleep(0.01)  # Reduce CPU usage

print("Exiting...")