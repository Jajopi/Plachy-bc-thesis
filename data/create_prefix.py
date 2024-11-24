#!/usr/bin/env python3

import sys
import os

def create_prefix(args):
    if len(args) != 3:
        print("specify FILE_NAME and PREFIX_LENGTH")
        return
    
    input_file_name = args[1]
    length = int(args[2])
    print(input_file_name, length)

    with open(input_file_name) as file:
        header = file.readline()
        content = file.read(length)
    
    file_name_parts = os.path.splitext(input_file_name)
    output_file_name = file_name_parts[0] + f"_prefix_{length}" + file_name_parts[1]
    with open(output_file_name, "w") as file:
        file.writelines((header, content, "\n"))

if __name__ == "__main__":
    create_prefix(sys.argv)
