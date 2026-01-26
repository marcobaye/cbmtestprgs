#!/usr/bin/env python3
import argparse
import d64util
import sys  # for sys.stderr and sys.exit

def show_help():
    # TODO:
    #d64.py create IMAGE             create new empty disc image
    #d64.py addfile IMAGE FILE       add file to image
    #d64.py extract [--full] [--p00] IMAGE [FILENUM]     extract files from disc image
    #d64.py create --typeswitches --force IMG.D64
    #d64.py addfile --forcetype --relsize IMG.D64 FILE
    #d64.py check IMG.D64
    print("""
Usage:
    d64.py [-h] [help]              show this help
    d64.py help [MODE]              show help about mode
    d64.py list                     list supported image formats
    d64.py [dir] IMAGE              show directory
    d64.py checkfile IMAGE FILENUM  check blocks of file
""")

def mode_checkfile():
    parser = argparse.ArgumentParser(allow_abbrev = False, description =
"""
This mode checks the block allocation of a single file.
The file must be specified by its directory index.
""")
    parser.add_argument("image", metavar="IMAGE.D64", help="Disk image file.")
    parser.add_argument("file_index", metavar="FILENUM", type=int, help="Dir index of file to check.")
    args = parser.parse_args(sys.argv[2:])
    image = d64util.DiskImage(args.image)
    image.check_file(args.file_index)

def one_arg(arg):
    if arg in ("help", "-h"):
        show_help()
    elif arg == "list":
        d64util.list_formats()
    else:
        d64util._process_file(arg)

def not_yet():
    print("not yet implemented")

def _main():
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    if len(sys.argv) == 2:
        one_arg(sys.argv[1])
    else:
        mode = sys.argv[1]
        if mode == "help":
            not_yet()
        elif mode == "dir":
            not_yet()
        elif mode == "create":
            not_yet()
        elif mode == "addfile":
            not_yet()
        elif mode == "extract":
            not_yet()
        elif mode == "checkfile":
            mode_checkfile()
        else:
            d64util._main()

if __name__ == "__main__":
    _main()
