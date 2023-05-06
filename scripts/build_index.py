#!/usr/bin/env python3

import argparse
import subprocess
import os
import pprint
from pathlib import Path
from timeit import default_timer as timer


def run_cmd(args):
    start = timer()
    proc = subprocess.run(args)
    if proc.returncode != 0:
        print(f"{args[0]} returned exit code {proc.returncode}")
    end = timer()
    print(f"{args[0]} {args[1]} took {(end - start)}")


VALID_STEPS = ["invert", "sort_unique", "permute_unitigs", "build"]


def build(args):

    cf_input_prefix = args.cuttlefish_input
    step_list = []

    invert_cmd = [
        args.bin_dir + "/fulgor",
        "invert",
        "-i",
        cf_input_prefix,
        "-g",
        str(args.working_mem),
        "-d",
        args.tmp_dir,
        "--verbose",
    ]
    step_list.append(invert_cmd)

    sort_unique_cmd = [
        args.bin_dir + "/fulgor",
        "sort-unique",
        "-i",
        cf_input_prefix,
        "-g",
        str(args.working_mem),
        "-d",
        args.tmp_dir,
        "--verbose",
    ]
    step_list.append(sort_unique_cmd)

    permute_unitigs_cmd = [
        args.bin_dir + "/fulgor",
        "permute-unitigs",
        "-i",
        cf_input_prefix,
        "-g",
        str(args.working_mem),
        "-d",
        args.tmp_dir,
        "--verbose",
    ]
    step_list.append(permute_unitigs_cmd)

    build_cmd = [
        args.bin_dir + "/fulgor",
        "build",
        "-i",
        cf_input_prefix,
        "-k",
        str(args.k),
        "-m",
        str(args.m),
        "-d",
        args.tmp_dir,
        "--verbose",
    ]
    step_list.append(build_cmd)

    step_idx = VALID_STEPS.index(args.s)
    ns = len(VALID_STEPS)
    for i, (step_name, step_cmd) in enumerate(
        zip(VALID_STEPS[step_idx:], step_list[step_idx:])
    ):
        print(f"executing step {i+1}/{ns} : {step_name}")
        print(f"cmd : {' '.join(step_cmd)}")
        if not args.dry_run:
            run_cmd(step_cmd)


def build_step(s):
    if s not in VALID_STEPS:
        raise ValueError
    return s


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="build_color_index", description="Builds the fulgor index"
    )
    parser.add_argument("cuttlefish_input")
    parser.add_argument(
        "--bin-dir", default=".", help="directory where executables exist"
    )
    parser.add_argument("-k", default=31)
    parser.add_argument("-m", default=17)
    parser.add_argument("--tmp-dir", default=".")
    parser.add_argument("-d", "--dry-run", action="store_true")
    parser.add_argument("-g", "--working-mem", default=8)
    parser.add_argument(
        "-s",
        type=build_step,
        help="which step to start from; these are, in order, [invert, sort_unique, permute_unitigs, build]",
        default="invert",
    )

    args = parser.parse_args()

    if args.bin_dir == ".":
        args.bin_dir = os.getcwd()

    build(args)
