import subprocess


def run_evolver(se_input_filepath: str, evolver_filepath: str = "evolver "):
    p = subprocess.run([evolver_filepath, se_input_filepath], stdout=subprocess.DEVNULL)
    return p.returncode == 0
