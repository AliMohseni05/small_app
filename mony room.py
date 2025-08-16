import random
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from IPython.display import HTML

# Initialize parameters
num_people = 50
initial_money = 100
max_turns = 100  # Maximum number of turns to simulate

# Initialize people and money
people = list(range(num_people))
money = [initial_money] * num_people
active = [True] * num_people  # Track active participants

# Store history for visualization
history = [money.copy()]

# Simulation loop
for turn in range(max_turns):
    # Create a list of active participants with money
    active_people = [i for i in people if active[i] and money[i] > 0]
    
    # If less than 2 active people, stop simulation
    if len(active_people) < 2:
        break
    
    # Each active person gives money to a random other active person
    for i in active_people:
        # Find a valid recipient (different person, also active)
        recipient = random.choice(active_people)
        while recipient == i:
            recipient = random.choice(active_people)
        
        # Exchange money (giver loses 1, recipient gains 1)
        money[i] -= 1
        money[recipient] += 1
    
    # Deactivate people with no money
    for i in active_people:
        if money[i] <= 0:
            active[i] = False
    
    # Record current state
    history.append(money.copy())

# Visualization setup
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [3, 1]})
plt.subplots_adjust(hspace=0.4)

# 1. Animated histogram of wealth distribution
def update_hist(frame):
    ax1.clear()
    current_money = history[frame]
    active_money = [m for m in current_money if m > 0]
    
    ax1.hist(active_money, bins=30, color='skyblue', edgecolor='black')
    ax1.set_title(f'Wealth Distribution (Turn {frame}/{len(history)-1})')
    ax1.set_xlabel('Money')
    ax1.set_ylabel('Number of People')
    ax1.set_xlim(0, max(max(history[0]), max(history[-1])))
    ax1.grid(alpha=0.3)
    
    # Add statistics
    avg = np.mean(active_money)
    std = np.std(active_money)
    ax1.axvline(avg, color='r', linestyle='dashed', linewidth=1)
    ax1.text(0.98, 0.95, f'Active: {len(active_money)}/{num_people}\nAvg: {avg:.1f} Std: {std:.1f}',
             transform=ax1.transAxes, ha='right', va='top',
             bbox=dict(facecolor='white', alpha=0.8))

# 2. Time series of wealth metrics
def plot_metrics():
    ax2.clear()
    active_counts = [len([m for m in state if m > 0]) for state in history]
    max_wealth = [max(state) for state in history]
    min_wealth = [min([m for m in state if m > 0] or [0]) for state in history]
    avg_wealth = [np.mean([m for m in state if m > 0]) for state in history]
    
    turns = list(range(len(history)))
    ax2.plot(turns, active_counts, 'b-', label='Active People')
    ax2.plot(turns, max_wealth, 'r-', label='Max Wealth')
    ax2.plot(turns, min_wealth, 'g-', label='Min Wealth (active)')
    ax2.plot(turns, avg_wealth, 'm--', label='Avg Wealth')
    
    ax2.set_title('Wealth Metrics Over Time')
    ax2.set_xlabel('Turn')
    ax2.set_ylabel('Value')
    ax2.legend(loc='upper right')
    ax2.grid(alpha=0.3)
    ax2.set_xlim(0, len(history)-1)

# Create animation
ani = FuncAnimation(fig, lambda frame: (update_hist(frame), plot_metrics()), 
                    frames=len(history), interval=200, repeat=False)

# Display the animation
plt.close()
HTML(ani.to_jshtml())