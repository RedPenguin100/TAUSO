import logging
import time

logger = logging.getLogger(__name__)


class Timer:
    def __init__(self, name="Task"):
        self.name = name

    def __enter__(self):
        self.start_time = time.time()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end_time = time.time()
        self.elapsed_time = self.end_time - self.start_time
        logger.info("[%s] finished in %.4fs", self.name, self.elapsed_time)
