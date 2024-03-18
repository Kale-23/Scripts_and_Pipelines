#! /usr/bin/env python3

import os
import argparse
import shutil

def delete_files_not_sf(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if not file.endswith('.sf'):
                os.remove(os.path.join(root, file))

def remove_empty_directories(directory):
    for root, dirs, _ in os.walk(directory, topdown=False):
        for d in dirs:
            dir_path = os.path.join(root, d)
            try:
                if not os.listdir(dir_path):
                    os.rmdir(dir_path)
            except OSError as e:
                print(f"Error: {e} - Skipping removal of '{dir_path}'")

def create_tsv(directory):
    outlist = []
    outlist2 = []
    with open(f'{directory}/libraries.tsv', 'w') as tsv_file:
        with open(f'{directory}/libraries_with_context.tsv', 'w') as tsv_file2:
            #tsv_file.write('Directory Name\tType\n')
            for root, dirs, _ in os.walk(directory):
                for d in dirs:
                    if "SG" in d:
                    #if d.startswith('SG'):
                        outlist2.append([0, f'{d}\tslime\t{os.path.basename(root)}\n'])
                        outlist.append([0, f'{d}\tslime\n'])
                    elif "SK" in d:
                    #elif d.startswith('SK'):
                        outlist2.append([1, f'{d}\tskin\t{os.path.basename(root)}\n'])
                        outlist.append([0, f'{d}\tskin\n'])
            outlist.sort()
            outlist2.sort()
            for item in outlist:
                tsv_file.write(item[1])
            for item in outlist2:
                tsv_file2.write(item[1])

def move_directories_to_mapping(directory):
    mapping_dir = os.path.join(directory, 'mapping')
    os.makedirs(mapping_dir, exist_ok=True)
    for item in os.listdir(directory):
        if item == "mapping":
            continue
        item_path = os.path.join(directory, item)
        if os.path.isdir(item_path):
            shutil.move(item_path, os.path.join(mapping_dir, item))

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str, required=True, help="Directory of Salmon quant outputs")
args = parser.parse_args()

delete_files_not_sf(args.d)
remove_empty_directories(args.d)
create_tsv(args.d)
move_directories_to_mapping(args.d)
