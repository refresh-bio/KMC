#!/usr/bin/env python3
import sys
import os
import platform
import os
def is_windows():
    return 'win' in platform.system().lower()
def is_linux():
    return 'linux' in platform.system().lower()
def is_mac():
    return 'darwin' in platform.system().lower()

if is_windows():
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'x64', 'Release'))
elif is_linux() or is_mac():
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
    pass
else:
    raise RuntimeError("system other than linux, mac os and windows was detected which is now not supported.")    
