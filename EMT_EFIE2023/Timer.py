from datetime import datetime as dt

class Timer:
    def __init__(self, name=None):
        self.name = name
        self.start_time = None
        self.end_time = None

    def __enter__(self):
        self.start_time = dt.now()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end_time = dt.now()
        self.elapsed_time = self.end_time - self.start_time
        if self.name:
            print(f"Timer '{self.name}' elapsed time: {self.elapsed_time}")
        else:
            print(f"Elapsed time: {self.elapsed_time}")