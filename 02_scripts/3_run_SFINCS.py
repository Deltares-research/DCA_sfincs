import os
import subprocess

base_dir = r"p:\11209905-dca-sfincs-river\01_models"

modelnames = ["KoldingA_PAK_res25_sub5", "KoldingA_PAK_res25_sub5_riv1m", "KoldingA_PAK_res25_sub5_riv3m" ]

for modelname in modelnames:

    print(f"Running: {modelname}")
    run_dir = os.path.join(base_dir, modelname)

    # Specify the path to your .bat file
    bat_file_path = os.path.join(run_dir, 'run.bat') 


    # # Run the updated .bat file
    os.chdir(run_dir)

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