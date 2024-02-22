import os
import subprocess

base_dir = r"p:\11209905-dca-sfincs-river\01_models\KoldingA_PAK_sbg20"

# Specify the path to your .bat file
bat_file_path = os.path.join(base_dir, 'run.bat') 


# # Run the updated .bat file
os.chdir(base_dir)

# Open the .bat file in a subprocess, capturing both standard output and standard error
process = subprocess.Popen(bat_file_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

# Print the output in real-time
while True:
    output = process.stdout.readline()
    if not output and process.poll() is not None:
        break
    print(output.strip())

# Wait for the process to complete
process.wait()