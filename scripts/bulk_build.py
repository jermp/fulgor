#!/usr/bin/env python3
import pathlib
import subprocess
import time

# Define datasets as a list of dictionaries
DATASETS = [
    {
        "directory": "~/allthebacteria",
        "basename": "salmonella_enterica",
        "outputfilename": "se",
        "amounts": ["100k", "300k", "500k"],
    }, {
        "directory": "~/allthebacteria",
        "basename": "escherichia_coli",
        "outputfilename": "ec",
        "amounts": ["100k", "200k", "300k"],
    }, {
        "directory": "~/allthebacteria",
        "basename": "",
        "outputfilename": "atb",
        "amounts": ["top10"],
    }, {
        "directory": "~/hprc",
        "basename": "",
        "outputfilename": "hprc",
        "amounts": ["all"],
    }, {
        "directory": "~/gutbacteria",
        "basename": "",
        "outputfilename": "gb",
        "amounts": ["all"],
    },
]


def run_fulgor_build():
    now = time.strftime('%y%m%d_%H%M', time.localtime(int(time.time())))
    logs_dir = f"~/AMB/logs/build-{now}"

    pathlib.Path(logs_dir).expanduser().mkdir(parents=True, exist_ok=True)

    for ds in DATASETS:
        directory = ds["directory"]
        basename = ds["basename"]
        output = ds["outputfilename"]

        pathlib.Path(f"{directory}/AMB/indexes").expanduser().mkdir(parents=True, exist_ok=True)
        pathlib.Path("tmp").mkdir(parents=True, exist_ok=True)

        for amount in ds["amounts"]:
            # Construct relative paths
            list_file = f"{directory}/filenames/{basename}-{amount}.txt"
            if basename == "":
                list_file = f"{directory}/filenames/{amount}.txt"
            out_file = f"{directory}/AMB/indexes/{output}-{amount}"
            hfur_log_file = f"{logs_dir}/{output}-{amount}.hfur.build.log"
            mfur_log_file = f"{logs_dir}/{output}-{amount}.mfur.build.log"

            # Construct the target command
            hfur_command = (
                f"/usr/bin/time -v ./fulgor build "
                f"-l {list_file} "
                f"-o {out_file} "
                f"-k 31 -m 19 -g 32 "
                f"-d tmp "
                f"-t 48 "
                f"--verbose 2>&1 | tee {hfur_log_file}"
            )

            print(f"Command: {hfur_command}\n")

            # Run synchronously and wait for process completion
            subprocess.run(
                hfur_command,
                shell=True,
                executable="/bin/bash"
            )

            '''
            # Construct the target command
            mfur_command = (
                f"/usr/bin/time -v ./fulgor build "
                f"-l {list_file} "
                f"-o {out_file} "
                f"-k 31 -m 19 -g 32 "
                f"-d tmp "
                f"-t 48 "
                f"--meta "
                f"--verbose 2>&1 | tee {mfur_log_file}"
            )

            print(f"Command: {mfur_command}\n")

            # Run synchronously and wait for process completion
            subprocess.run(
                mfur_command,
                shell=True,
                executable="/bin/bash"
            )
            '''


if __name__ == "__main__":
    run_fulgor_build()
