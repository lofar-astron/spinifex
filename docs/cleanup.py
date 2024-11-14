#!/usr/bin/env python3

#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

import os

file_dir = os.path.dirname(os.path.realpath(__file__))

clean_dir = os.path.join(file_dir, "source", "source_documentation")
print(f"Cleaning.. {clean_dir}/*")

if not os.path.exists(clean_dir):
    exit()

for file_name in os.listdir(clean_dir):
    file = os.path.join(clean_dir, file_name)
    
    if file_name == "index.rst":
        continue

    print(f"Removing.. {file}")
    os.remove(file)
