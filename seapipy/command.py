import subprocess


def run_evolver(se_input_filepath: str, evolver_filepath: str = "evolver ") -> int:
    """
    Run the Surface Evolver file using the Surface Evolver software

    :param se_input_filepath: Path to Surface Evolver file
    :rtype se_input_filepath: str
    :param evolver_filepath: Path to the Surface Evolver interpreter
    :rtype evolver_filepath: str
    :return: Return code of the subprocess.run() function
    :rtype: int
    """
    p = subprocess.run([evolver_filepath, se_input_filepath], stdout=subprocess.DEVNULL)
    return p.returncode == 0
