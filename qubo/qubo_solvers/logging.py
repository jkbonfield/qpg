import logging
import sys


class InfoFilter(logging.Filter):
    def filter(self, record):
        return record.levelno in [logging.INFO]
    

class WarnErrorFilter(logging.Filter):
    def filter(self, record):
        return record.levelno in [logging.WARN, logging.ERROR]


def get_logger(name: str):
    root = logging.getLogger(name)
    root.setLevel(logging.INFO)
    if (root.hasHandlers()):
        root.handlers.clear()

    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s","%H:%M:%S")

    out_handler = logging.StreamHandler(sys.stdout)
    out_handler.setLevel(logging.INFO)
    out_handler.addFilter(InfoFilter())
    out_handler.setFormatter(formatter)

    err_handler = logging.StreamHandler(sys.stderr)
    err_handler.setLevel(logging.WARN)
    err_handler.addFilter(WarnErrorFilter())
    err_handler.setFormatter(formatter)

    root.addHandler(out_handler)
    root.addHandler(err_handler)
    return root