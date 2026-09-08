import logging
import time

logger = logging.getLogger(__name__)


class Timer:
    def __init__(self, name="Task", log=True):
        self.name = name
        self.log = log

    def __enter__(self):
        self.start_time = time.time()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end_time = time.time()
        self.elapsed_time = self.end_time - self.start_time
        if self.log:  # a caller that reports the time itself asks for log=False
            logger.info("[%s] finished in %.4fs", self.name, self.elapsed_time)
