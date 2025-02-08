from datetime import datetime


def log(message):
    timestamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    print(
        f"{timestamp} {message}"
    )  # You can also write to a log file if needed
