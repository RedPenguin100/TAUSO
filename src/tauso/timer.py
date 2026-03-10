import time


class Timer:
    def __init__(self, name="Task"):
        """Initialize with a task name for clear console output."""
        self.name = name

    def __enter__(self):
        self.start_time = time.time()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end_time = time.time()
        self.elapsed_time = self.end_time - self.start_time
        # Automatically print the duration when the block finishes
        print(f"[{self.name}] finished in {self.elapsed_time:.4f}s")
