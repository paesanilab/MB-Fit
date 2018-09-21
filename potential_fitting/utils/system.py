import subprocess

from potential_fitting.exceptions import CommandNotFoundError, CommandExecutionError

def call(command, *args):
    """
    Performs a system call with the given command and arguments

    A CommandNotFoundError will be raised if the command is not found (exit code 127)

    A CommandExecutionError will be raised if the exit code is not 0

    Args:
        command - the command to be executed
        args    - any arguments to follow the command

    Returns:
        None
    """

    try:
        subprocess.run([command, *args], stderr = subprocess.PIPE, check = True)
    except subprocess.CalledProcessError as e:
        if e.returncode == 127:
            raise CommandNotFoundError(command)
        raise CommandExecutionError(command, e.cmd, e.returncode, e.stderr)
    except FileNotFoundError:
        raise CommandNotFoundError(command)
