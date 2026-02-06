#!/usr/bin/env python3
import argparse
import d64util
import sys  # for sys.stderr and sys.exit

def show_help():
    # TODO:
    #d64.py create IMAGE             create new empty disc image
    #d64.py addfile IMAGE FILE       add file to image
    #d64.py extract [--full] [--p00] IMAGE [FILENUM]     extract files from disc image
    #d64.py check IMG.D64
    print("""
Usage:
    d64.py [-h] [help]              show this help
    d64.py help [MODE]              show help about mode
    d64.py list                     list supported image formats
    d64.py [dir] IMAGE              display directory
    d64.py create IMAGE             create new image file
    d64.py checkfile IMAGE FILENUM  check blocks of file
    d64.py add IMAGE FILE           add file to image
    d64.py bam IMAGE                display block availability map
    d64.py errors IMAGE             display optional error chunk
""")

def mode_create():
    #sys.argv[1] += " create"
    parser = argparse.ArgumentParser(allow_abbrev = False, description =
"""
This mode creates a new empty disk image file.
""")
    # TODO: add "create" string so it is included in help message.
    #parser.add_argument("create", action=None, help="mode")
    parser.add_argument("image", metavar="IMAGEFILE", help="disk image file to create")
    parser.add_argument("-t", "--type", metavar="TYPE", type=str, default=None, help="type of image")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite existing file")
    parser.add_argument("-e", "--error", action="store_true", help="create with error chunk")
    parser.add_argument("--tracks", metavar="TRACKS", help="number of tracks for CMD native partition")
    args = parser.parse_args(sys.argv[2:])
    #print("image:", args.image)
    #print("type:", args.type)
    #print("force:", args.force)
    print("not yet")

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

def mode_delete():
    parser = argparse.ArgumentParser(allow_abbrev = False, description =
"""
This mode deletes a single file from the image.
The file must be specified by its directory index.
""")
    parser.add_argument("image", metavar="IMAGE.D64", help="Disk image file.")
    parser.add_argument("file_index", metavar="FILENUM", type=int, help="Dir index of file to delete.")
    args = parser.parse_args(sys.argv[2:])
    image = d64util.DiskImage(args.image, d64util.ImgMode.WRITEBACK)
    image.delete_file(args.file_index)
    image.writeback()   # flush to file

def mode_add():
    parser = argparse.ArgumentParser(allow_abbrev = False, description =
"""
This mode adds a file to the image.
""")
    parser.add_argument("image", metavar="IMAGE.D64", help="Disk image file.")
    parser.add_argument("file", metavar="FILE", help="File to add to image.")
    # TODO: add "--forcefiletype"
    # TODO: add "--relsize"
    args = parser.parse_args(sys.argv[2:])
    image = d64util.DiskImage(args.image, d64util.ImgMode.WRITEBACK)
    file = open(args.file, "rb")
    filebody = file.read(16 * 1024 * 1024)  # 16 MiB ought to be enough for anybody
    image.add_file(args.file, filebody)
    image.writeback()   # flush to file

def mode_bam():
    parser = argparse.ArgumentParser(allow_abbrev = False, description =
"""
This mode displays the image's "block availability map".
CAUTION, the format of the output varies between image types!
""")
    parser.add_argument("image", metavar="IMAGE.D64", help="Disk image file.")
    args = parser.parse_args(sys.argv[2:])
    image = d64util.DiskImage(args.image)
    image.bam_display()

def mode_errors():
    parser = argparse.ArgumentParser(allow_abbrev = False, description =
"""
This mode displays the image's error chunk, if there is one.
""")
    parser.add_argument("image", metavar="IMAGE.D64", help="Disk image file.")
    args = parser.parse_args(sys.argv[2:])
    image = d64util.DiskImage(args.image)
    image.errorchunk_display()

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
            mode_create()
        elif mode == "checkfile":
            mode_checkfile()
        elif mode == "delete":
            mode_delete()
        elif mode == "add":
            mode_add()
        elif mode == "extract":
            not_yet()
        elif mode == "bam":
            mode_bam()
        elif mode == "errors":
            mode_errors()
        else:
            d64util._main()

if __name__ == "__main__":
    _main()
