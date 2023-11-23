#!/usr/bin/env python3

import sys
import d64util
import argparse

def get_dir_art(diskimage):
    """
    Read dir art from image and return as list of byte sequences.
    """
    artdisk = d64util.DiskImage(diskimage)
    art = []
    for entry in artdisk.read_directory_entries(include_invisible=False):
        bin30 = entry[1]    # [0] would be index, we only want the data
        art.append(bin30[3:19]) # and we only keep the name field
    print("Found %d entries of dir art." % len(art))
    return art

def apply_dir_art(art, targetdisk, targetskip, use_all_art, add_spacer):
    """
    Do the actual merge of original directory and dir art.
    """
    # open image in writeback mode so it does not get clobbered on error:
    img = d64util.DiskImage(targetdisk, d64util.ImgMode.WRITEBACK)
    # get current directory contents
    olddir = []
    for entry in img.read_directory_entries(include_invisible=False):
        bin30 = entry[1]    # [0] would be index, we only want the data
        olddir.append(bin30)
    print("Original directory has %d entries." % len(olddir))
    # now build new directory
    newdir = []
    applied_art = 0 # counter for debug output
    # part 1: if wanted, keep a number of entries unchanged
    if targetskip:
        print("Skipping %d entries of target dir." % targetskip)
        newdir.extend(olddir[:targetskip])
        olddir = olddir[targetskip:]
    # part 2: handle as many two-char file names as there are
    while len(olddir):
        # offsets 0/1/2 hold type and blockcount, shift-space-padded file name starts at offset 3:
        # if name has less than two chars, stop
        if olddir[0][4] == 0xa0:
            break
        # if name has more than two chars, stop
        if olddir[0][5] != 0xa0:
            break
        # file name is two chars long, so process:
        dir_entry = olddir.pop(0)
        if len(art):
            # we still have art, so use it:
            art_entry = art.pop(0)
            dir_entry[5:19] = art_entry[2:16]   # overwrite last 14 chars of filename
            applied_art += 1
        newdir.append(dir_entry)
    # part 3: if there is still some art left, decide what to do with it
    if use_all_art:
        while len(art):
            # create dummy entries for dir art
            art_entry = art.pop(0)
            dir_entry = img.build_dummy_dir_entry(art_entry[:16])   # all 16 chars of art
            applied_art += 1
            newdir.append(dir_entry)
        print("Applied all %d entries of dir art." % applied_art)
    else:
        print("Applied %d of %d entries of dir art." % (applied_art, applied_art + len(art)))
    # now we are done with dir art!
    # part 4: if we still have original directory entries left, add spacer:
    if len(olddir) and add_spacer:
        newdir.append(img.build_dummy_dir_entry(b"----------------"))
        #newdir.append(img.build_dummy_dir_entry(b"----------------"))
        #newdir.append(img.build_dummy_dir_entry(b"----------------"))
    # part 5: append remaining original directory entries unchanged
    newdir.extend(olddir)
    # part 6: write new directory
    img.write_directory(newdir)
    img.writeback() # flush to file
    # ...aaand we're done!

def main():
    if len(sys.argv) == 1:
        sys.argv.append("-h")   # if run without arguments, show help instead of complaining!
    parser = argparse.ArgumentParser(allow_abbrev = False, description =
"""
This program copies directory art from one disk image to another.
It works with 1541, 1571 and 1581 images.
Only two-char file names will be changed.
""")
    parser.add_argument("--artskip", metavar="NUM", type=int, default=0, help="Number of entries to skip at start of dir art.")
    parser.add_argument("--targetskip", metavar="NUM", type=int, default=0, help="Number of entries to skip at start of target dir.")
    parser.add_argument("--keepsize", action="store_true", help="Stop applying dir art at end of target dir.")
    parser.add_argument("--spacer", action="store_true", help="Add '----------------' spacer.")
    parser.add_argument("-d", "--debug", action="count", default=0, help="Increase debugging output.")
    parser.add_argument("artdisk", metavar="ARTDISK.D64", help="Disk image file to read directory art from.")
    parser.add_argument("targetdisk", metavar="TARGETDISK.D64", help="Disk image file to apply directory art to.")
    args = parser.parse_args()
    d64util._debuglevel += args.debug

    # step 1: read dir art from disk image
    dir_art = get_dir_art(args.artdisk)
    # step 2: if wanted, skip over the first few entries 
    if args.artskip:
        print("Skipping first %d entries of dir art." % args.artskip)
        dir_art = dir_art[args.artskip:]
        print("Now %d entries of dir art left." % len(dir_art))
    # step 3: put directory art on target disk image
    apply_dir_art(art=dir_art, targetdisk=args.targetdisk, targetskip=args.targetskip, use_all_art = not args.keepsize, add_spacer = args.spacer)

if __name__ == "__main__":
    main()
