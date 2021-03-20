#!/usr/bin/env python3

import os

install_req = "pip3 install -r requirements.txt"
os.system(install_req)

install_setup = "pip3 install . --upgrade"
os.system(install_setup)
