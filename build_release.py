#!/usr/bin/env python3
import subprocess
import fileinput
import os
import sys

def get_ver():
    with open("kmc_core/defs.h") as f:
        for line in f.readlines():
            line = line.strip()
            if "Version:" in line:
                return line.split("Version:")[-1].strip()
    print("Error: cannot read R_NOMAD_VERSION")
    sys.exit(1)

def get_os():
    if os.name == 'nt':
        return 'windows'
    elif os.name == 'posix':
        if os.uname()[0] == 'Linux':
            return 'linux'
        elif os.uname()[0] == 'Darwin':
            return 'mac'
        else:
            print("Error: unknown os", os.uname()[0])
            sys.exit(1)
    else:
        print("Error: unknown os.name", os.name)
        sys.exit(1)

def get_hardware():
    if os.name == 'nt':
        print("Error: windows is currently unsupported in this script")
        sys.exit(1)
    elif os.name == 'posix':
        if os.uname()[4] == 'x86_64':
            return 'x64'
        elif os.uname()[4] == 'aarch64' or os.uname()[4] == 'arm64':
            return 'arm64'
        else:
            print("Error: unknown hardware", os.uname()[4])
            sys.exit(1)
    else:
        print("Error: unknown os.name", os.name)
        sys.exit(1)


def run_cmd(cmd):    
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()

system = get_os()
hardware = get_hardware()

if system == 'windows':
    print("Error: this script does not currently support windows build")
    sys.exit(1)

ver = get_ver()

print(f"building\n\tVersion: {ver}\n\tOperating system: {system}\n\tHardware: {hardware}")

run_cmd("git submodule init")
run_cmd("git submodule update")
run_cmd("make clean")
run_cmd("make -j32 kmc kmc_dump kmc_tools")

run_cmd(f"tar -c bin/kmc bin/kmc_tools bin/kmc_dump bin/libkmc_core.a include/kmc_runner.h | pigz > KMC{ver}.{system}.{hardware}.tar.gz")
