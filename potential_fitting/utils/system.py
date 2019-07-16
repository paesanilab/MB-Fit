# external package imports
import subprocess
from enum import Enum

# absolute module imports
from potential_fitting.exceptions import CommandNotFoundError, CommandExecutionError

def call(command, *args, in_file = None, out_file = None):
    """
    Performs a system call with the given command and arguments.

    A CommandNotFoundError will be raised if the command is not found (exit code 127).

    A CommandExecutionError will be raised if the exit code is not 0 or 127.

    Args:
        command             - The command to be executed.
        args                - Any arguments to follow the command.
        in_file             - File handle to redirect input from.
        out_file            - File handle to redirect output to.

    Returns:
        None.
    """

    try:
        subprocess.run([command, *args], stdin = in_file, stdout = out_file, stderr = subprocess.PIPE, check = True)
    except subprocess.CalledProcessError as e:
        if e.returncode == 127:
            raise CommandNotFoundError(command) from None
        raise CommandExecutionError(command, e.cmd, e.returncode, e.stderr) from None
    except FileNotFoundError:
        raise CommandNotFoundError(command) from None

def format_print(string, bold = False, italics = False, color = None):
    if bold:
        string = '\33[1m' + string + '\33[0m'
    if italics:
        string = '\33[3m' + string + '\33[0m'
    if color is not None:
        string = '\33[{}m'.format(color.value) + string + '\33[0m'

    print(string)

class Color(Enum):
    RED = 31
    GREEN = 92
    BLUE = 34
    NORMAL = 0
    YELLOW = 33


