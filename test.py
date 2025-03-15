def display_tasks(tasks):
    print("\nYour To-Do List:")
    for idx, task in enumerate(tasks, start=1):
        print(f"{idx}. {task['name']} {'✓' if task['completed'] else ''}")

def main():
    tasks = []
    while True:
        action = input("\nWould you like to [add], [view], [complete], or [remove] a task? (Type 'quit' to exit): ").strip().lower()
        
        if action == "add":
            task_name = input("Enter the task name: ")
            tasks.append({"name": task_name, "completed": False})
            print(f"Added task: {task_name}")

        elif action == "view":
            display_tasks(tasks)

        elif action == "complete":
            task_num = int(input("Enter the task number to mark as complete: ")) - 1
            if 0 <= task_num < len(tasks):
                tasks[task_num]['completed'] = True
                print(f"Task {task_num + 1} marked as complete.")
            else:
                print("Invalid task number!")

        elif action == "remove":
            task_num = int(input("Enter the task number to remove: ")) - 1
            if 0 <= task_num < len(tasks):
                removed_task = tasks.pop(task_num)
                print(f"Removed task: {removed_task['name']}")
            else:
                print("Invalid task number!")

        elif action == "quit":
            print("Goodbye! Keep swimming!")
            break

        else:
            print("Sorry, I didn't catch that! Try again.")

if __name__ == "__main__":
    main()