import shlex
import subprocess

from .logger import log


class CommandExecutor:
    def execute(self, command):
        try:
            log(f"Executing command: {command}")
            result = subprocess.run(
                shlex.split(command), capture_output=True, text=True
            )
            return result.stdout if result.returncode == 0 else result.stderr
        except Exception as e:
            return f"Error executing command: {e}"
