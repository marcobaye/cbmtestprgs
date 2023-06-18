#!/usr/bin/env python3

import argparse
import sys  # for sys.stderr and sys.exit

# this library accepts disk images with error info, but atm does not honor it!

debuglevel = 1

def debug(minlevel, *unnamed, **named):
    """Helper function for debugging output."""
    if debuglevel >= minlevel:
        print("debug%d:" % minlevel, *unnamed, **named)

def popcount(integer):
    """Helper function to check bit fields in BAM."""
    counted = 0
    while integer:
        counted += 1
        integer = integer & (integer - 1)
    return counted

def _bam_offset_and_bit(sector):
    """
    Convert sector number to byte offset and bit value for accessing BAM bitmap.
    """
    # bitmap bytes are little-endian, lsb first:
    return sector >> 3, 1 << (sector & 7)

# TODO: use a ts tuple instead of track/sector, so ts and lba can be used interchangeably

################################################################################
# virtual base class

class d64(object):
    """
    This class describes a cbm disc image. There are subclasses for 1541/40track/1571/1581.
    """
    def __init__(self):
        raise Exception("Only subclasses (d41, d71, d81) can be instantiated!")

    def _populate(self, fh, body, readonly, writeback, writethrough):
        """
        Called after correct subclass has been instantiated.
        """
        self.fh = fh
        self.body = body
        self.readonly = readonly    # we don't really need to store this info
        self.writebackmode = writeback  # in three vars, but it might make
        self.writethroughmode = writethrough    # the code more readable.
        self.dirty = False  # no need to writeback yet

    def _check_track_num(self, track):
        """_check_track_num(int)

        Throw exception if track number is invalid for this image type.
        """
        if track < 1:
            raise Exception("Given track number was lower than 1.")
        if track > self.maxtrack:
            raise Exception("Exceeded maximum track number.")

    def _check_ts(self, track, sector):
        """_check_ts(int, int)

        Throw exception if track/sector address is invalid for this image type.
        """
        self._check_track_num(track)
        # now check sector number
        if sector < 0:
            raise Exception("Given sector number was negative.")
        if sector >= self.sectors_of_track(track):
            raise Exception("Exceeded maximum sector number of track %d." % track)

    def _virtualfn(self):
        raise Exception("BUG: A virtual function was called.")

    def sectors_of_track(self, track):
        """sectors_of_track(int) -> int

        Return number of sectors in given track.
        Sector numbers start at 0, so the maximum sector number is one less than this.
        """
        self._virtualfn()

    def ts_to_lba(self, track, sector):
        """ts_to_lba(int, int) -> int

        Convert track and sector to logical block addressing (0-based block number).
        """
        self._virtualfn()

    def read_ts(self, track, sector):
        """
        Read block given via track/sector number and return as bytes or bytearray.
        """
        debug(2, "Reading t%d s%d." % (track, sector))
        self._check_ts(track, sector)
        offset = 256 * self.ts_to_lba(track, sector)
        block = self.body[offset:offset+256]
        return block

    def write_ts(self, track, sector, data):
        """
        Write data to block given via track/sector number.
        """
        if len(data) != 256:
            raise Exception("block should be 256 bytes")
        debug(2, "Writing t%d s%d." % (track, sector))
        self._check_ts(track, sector)
        offset = 256 * self.ts_to_lba(track, sector)
        self.body[offset:offset+256] = data # will crash in readonly mode!
        if self.writebackmode:
            self.dirty = True
        else:
            # writethrough mode
            self.fh.seek(offset)
            self.fh.write(data)
            self.fh.flush()

    def writeback(self):
        """
        If data has been changed, write to file.
        """
        if not self.writebackmode:
            raise Exception("writeback was called outside of writeback mode")
        if self.dirty:
            debug(1, "Flushing to disk.")
            self.fh.seek(0)
            self.fh.write(self.body)
            self.fh.flush()
            self.dirty = False
        else:
            debug(1, "Nothing changed, no need to flush to disk.")

    def _check_totals(self, firsttrack, howmanytracks, total_offset_step, total_ts, bitsbytes, bits_offset_step, bits_ts=None):
        """
        Helper function to compare free blocks totals to free blocks bitmaps.
        firsttrack: number to add to zero-based loop counter to get meaningful debug/error messages
        howmanytracks: number of entries to check
        total_offset_step: offset of first total in block, and step size to get to the next total
        total_ts: track and sector where totals are stored
        bitsbytes: length of bitmap in bytes
        bits_offset_step: offset of first bitmap in block, and step size to get to the next bitmap
        bits_ts: track and sector where bitmaps are stored (if different from total_ts)
        """
        total_block = self.read_ts(*total_ts)
        if bits_ts == None:
            bits_block = total_block
        else:
            bits_block = self.read_ts(*bits_ts)
        total_offset, total_step = total_offset_step
        bits_offset, bits_step = bits_offset_step
        for entry in range(howmanytracks):
            track = entry + firsttrack
            total = total_block[total_offset + entry * total_step]
            sum = 0
            for i in range(bitsbytes):
                sum += popcount(bits_block[bits_offset + entry * bits_step + i])
            if total == sum:
                debug(9, "track %d: %d free blocks" % (track, sum))
            else:
                print("BAM error for track %d: counter says %d, bitfield says %d!" % (track, total, sum), file=sys.stderr)

    def _fill_free_blocks_dict(self, d, firsttrack, howmanytracks, offset_step, ts):
        """
        Helper function to read "free blocks" numbers from BAM.
        d: dictionary to put results in
        firsttrack: number to add to zero-based loop counter to get track number
        howmanytracks: number of entries to process
        offset_step: offset of first number in block, and step size to get to the next one
        ts: track and sector where numbers are stored
        The directory track is included in "all", but not in "shown"!
        """
        if "all" not in d:
            d["all"] = 0
        if "shown" not in d:
            d["shown"] = 0
        offset, step = offset_step
        block = self.read_ts(*ts)
        for entry in range(howmanytracks):
            track = entry + firsttrack
            freeblocks = block[offset + entry * step]
            d[track] = freeblocks
            d["all"] += freeblocks
            if track == self.directory_ts[0]:
                d["dir"] = freeblocks
            else:
                d["shown"] += freeblocks

    def bam_check(self):
        """
        Check BAM for internal consistency: Do counters match bit fields?
        This does not check allocation of files!
        """
        self._virtualfn()

    def read_header_fields(self):
        """
        Return disk name and five-byte "pseudo id".
        """
        track, sector = self.header_ts
        of = self.header_name_offset
        block = self.read_ts(track, sector)
        return block[of:of+16], block[of+18:of+23]

    def read_free_blocks(self):
        """
        Return dictionary of free blocks per track. There are three additional keys,
        "shown", "dir" and "all", where "shown"+"dir"="all"
        """
        self._virtualfn()

    def _release_block(self, ts, entry, totals_offset_step, totals_block, bitmaps_offset_step, bitmaps_block=None):
        """
        Helper function to release a single block in BAM.
        ts: track (for debugging output) and sector (to determine bit position)
        entry: which one of totals/bitmaps to use
        totals_offset_step: where to find totals
        totals_block: block with totals
        bitmaps_offset_step: where to find bitmaps
        bitmaps_block: block with bitmaps (if different from totals_block)
        """
        # process args
        track, sector = ts
        totals_offset, totals_step = totals_offset_step
        bitmaps_offset, bitmaps_step = bitmaps_offset_step
        if bitmaps_block == None:
            bitmaps_block = totals_block
        # calculate offsets
        byte_offset, bit_value = _bam_offset_and_bit(sector)
        bitmaps_offset += entry * bitmaps_step + byte_offset
        totals_offset += entry * totals_step
        if bitmaps_block[bitmaps_offset] & bit_value:
            raise Exception("Attempted to free a block (t%ds%d) that is already free." % (track, sector))
        else:
            bitmaps_block[bitmaps_offset] |= bit_value
            totals_block[totals_offset] += 1
        # TODO - compare totals to maximum for this track?
        #raise Exception("BAM is corrupt, totals do not match bitmap.")

    def _try_to_allocate(self, wanted_ts, entry, totals_offset_step, totals_ts, bitmaps_offset_step, bitmaps_ts=None, exact=True):
        """
        Helper function to allocate a single block in BAM.
        If block is available, allocate it and return t/s.
        If block is not available, return None.
        wanted_ts: track and sector to allocate
        entry: which one of totals/bitmaps to use
        totals_offset_step: where to find totals
        totals_ts: track and sector where totals are stored
        bitmaps_offset_step: where to find bitmaps
        bitmaps_ts: track and sector where bitmaps are stored (if different from totals_ts)
        exact: if False, function may allocate and return a different sector from this track
        """
        # process args
        track, wanted_sector = wanted_ts
        totals_offset, totals_step = totals_offset_step
        totals_block = self.read_ts(*totals_ts)
        bitmaps_offset, bitmaps_step = bitmaps_offset_step
        if bitmaps_ts == None:
            bitmaps_block = totals_block
        else:
            bitmaps_block = self.read_ts(*bitmaps_ts)
        # calculate offsets
        bitmaps_offset += entry * bitmaps_step
        totals_offset += entry * totals_step
        cand_sector = wanted_sector # we start the search with the wanted sector...
        num_sectors = self.sectors_of_track(track)  # ...and this is where we wrap around
        while True:
            byte_offset, bit_value = _bam_offset_and_bit(cand_sector)
            bitmap_byte_offset = bitmaps_offset + byte_offset
            # bit clear: sector is not available, i.e. is allocated or does not even exist
            # bit set: sector is available, i.e. can be allocated
            available = bool(bitmaps_block[bitmap_byte_offset] & bit_value)
            if available:
                # allocate block
                bitmaps_block[bitmap_byte_offset] &= (255 - bit_value)
                if totals_block[totals_offset]:
                    totals_block[totals_offset] -= 1
                else:
                    raise Exception("BAM is corrupt, totals do not match bitmap.")
                # write bam block(s)
                self.write_ts(*totals_ts, totals_block)
                if bitmaps_block != totals_block:
                    self.write_ts(*bitmaps_ts, bitmaps_block)
                debug(2, "Allocated t/s", track, cand_sector)
                return track, cand_sector   # block has been allocated
            debug(3, "t/s", track, cand_sector, "is not available")
            if exact:
                break   # fail (the wanted sector is not available)
            # try the next sector on this track:
            cand_sector += 1
            # if out of range, wrap around:
            if cand_sector >= num_sectors:
                cand_sector -= num_sectors
            # all done?
            if cand_sector == wanted_sector:
                break   # fail (no sectors left on this track)
        # if we arrive here, the request could not be satisfied
        return None # failure

    def release_blocks(self, set_of_ts):
        """
        Free all blocks given as t/s tuples.
        """
        self._virtualfn()

    def read_directory_entries(self, include_invisible=False):
        """
        Return directory entries as tuples:
        (index, raw entry (30 bytes), tuple might grow in future...)
        Setting include_invisible to True yields even the empty entries.
        """
        track, sector = self.directory_ts
        all_used_ts = set() # for sanity check
        entry_number = 0    # index, so caller can unambiguously reference each entry
        while track:
            # sanity check
            if (track, sector) in all_used_ts:
                raise Exception("Directory loops back to itself, please check disk image!")
            all_used_ts.add((track, sector))
            # read a directory sector
            block = self.read_ts(track, sector)
            track = block[0]
            sector = block[1]
            if track == 0:
                debug(2, "Last! Link is t%d s%d." % (track, sector))
            readidx = 2
            for i in range(8):
                bin30 = block[readidx:readidx+30]
                readidx += 32
                if bin30[0] or include_invisible:
                    yield entry_number, bin30   # CAUTION, maybe more fields will get added in future!
                entry_number += 1
        # "track" is zero, so there is no next block
        if sector != 255:
            print("WARNING: sector value of final dir block is %d instead of 255!" % sector, file=sys.stderr)

    def write_directory(self, new_dir): # TODO: add flag for "accept oversized dirs and use other tracks"
        """
        Overwrite directory with data given as list of 30-byte entries.
        """
        if len(new_dir) > self.std_max_dir_entries:
            raise Exception("New dir has too many entries (%d > %d)." % (len(new_dir), self.std_max_dir_entries))
        track, sector = self.directory_ts
        all_used_ts = set() # for sanity check
        init_link_ptrs = False  # flag needed when growing directory (to init fresh block)
        while track:
            # sanity check
            if (track, sector) in all_used_ts:
                raise Exception("Directory loops back to itself, please check disk image!")
            all_used_ts.add((track, sector))
            # read a directory sector
            block = self.read_ts(track, sector)
            if block[0] == 0:
                debug(2, "Last! Link is t%d s%d." % (block[0], block[1]))
            if init_link_ptrs:
                block[0:2] = 0x00, 0xff # init fresh dir block
            writeidx = 2
            # overwrite all eight entries
            for i in range(8):
                if len(new_dir):
                    entry = new_dir.pop(0)
                else:
                    entry = bytes(30)
                block[writeidx:writeidx+30] = entry
                writeidx += 32
            # check current link ptr (may get modified)
            next_track, next_sector = block[0:2]
            # now for the interesting part, enlarging/shrinking directory:
            if len(new_dir) and next_track == 0:
                # we have more data but no block to put it
                block[0:2] = self.get_new_block(track, sector) # get a new block and let current block point to it
                init_link_ptrs = True   # make sure link pointers get overwritten with 00/ff from now on
            elif len(new_dir) == 0 and next_track:
                # we have more block(s) but no data for them
                self.free_block_chain(next_track, next_sector)  # free all following blocks
                block[0:2] = 0x00, 0xff # mark this block as final
            # write back
            self.write_ts(track, sector, block)
            # get (potentially modified) link ptr
            track, sector = block[0:2]
        # "track" is zero, so there is no next block
        if sector != 255:
            print("WARNING: sector value of final dir block is %d instead of 255!" % sector, file=sys.stderr)
        # sanity check
        if len(new_dir):
            sys.exit("BUG: still data left after writing dir!")

    def build_dummy_dir_entry(self, namebytes=b""):
        """
        Return a 30-byte dummy directory entry (DEL).
        A name can be given, but must be bytes or bytearray!
        """
        # pad name with shift-space:
        name = namebytes + b"\xa0\xa0\xa0\xa0\xa0\xa0\xa0\xa0\xa0\xa0\xa0\xa0\xa0\xa0\xa0\xa0"
        entry = b"\x80" + bytes(self.header_ts) + name[:16] + bytes(11)
        # 0x80 is for a visible DEL entry,
        # tt/ss point to header block (so VALIDATE does not complain?),
        # name is 16 bytes of filename, padded with shift-space,
        # and then there are 11 bytes for various (DOS/GEOS) purposes, all 0 here
        return entry

    def get_new_block(self, prev_track, prev_sector):
        """
        Find and allocate a new block after given block, return t/s.
        """
        cand_track = prev_track
        # to get a candidate sector, add interleave:
        if cand_track == self.directory_ts[0]:
            cand_sector = prev_sector + self.std_directory_interleave   # may be out of range, see check below!
        else:
            cand_sector = prev_sector + self.std_file_interleave    # may be out of range, see check below!
            raise Exception("Block allocation is only supported for directory track at the moment!")
            # ...and to fully support files, we'd also need to iterate over all tracks of disk!
        # if out of range, wrap around:
        num_sectors = self.sectors_of_track(cand_track)
        if cand_sector >= num_sectors:
            cand_sector -= num_sectors
        # just do what the DOS does, allocate the next possible sector:
        ts = self.try_to_allocate(cand_track, cand_sector, exact=False)
        if ts == None:
            raise Exception("Directory track is full!")
        return ts   # return track and sector

    def free_block_chain(self, track, sector):
        """
        Free all connected blocks, following the chain of link pointers.
        """
        all_used_ts = set()
        while track:
            # sanity check
            if (track, sector) in all_used_ts:
                raise Exception("Block chain loops back to itself, please check disk image!")
            all_used_ts.add((track, sector))
            # read block, free it, go on with link
            block = self.read_ts(track, sector)
            track, sector = block[0:2]
        self.release_blocks(all_used_ts)

################################################################################
# 1581 class

class _d81(d64):
    name = "1581"
    blocks_total = 3200
    maxtrack = 80
    header_ts = (40, 0)
    header_name_offset = 4  # diskname and five-byte "pseudo id"
    directory_ts = (40, 3)
    std_max_dir_entries = 296   # for writing directory
    std_directory_interleave = 1    # 1581 uses interleave 1 because of track cache
    std_file_interleave = 1 # 1581 uses interleave 1 because of track cache

    def __init__(self):
        # this inhibits the base class's exception...
        pass    # ...but there is nothing to do!

    def sectors_of_track(self, track):
        self._check_track_num(track)
        return 40

    def ts_to_lba(self, track, sector):
        self._check_ts(track, sector)
        return (track - 1) * 40 + sector

    def bam_check(self):
        debug(1, "Checking 1581 BAM")
        self._check_totals(1, 40, (16, 6), (40, 1), 5, (17, 6))
        self._check_totals(41, 40, (16, 6), (40, 2), 5, (17, 6))

    def read_free_blocks(self):
        d = dict()
        self._fill_free_blocks_dict(d, 1, 40, (16, 6), (40, 1))
        self._fill_free_blocks_dict(d, 41, 40, (16, 6), (40, 2))
        return d

    def release_blocks(self, set_of_ts):
        bamblock1 = self.read_ts(40, 1)
        bamblock2 = self.read_ts(40, 2)
        dirty1 = False
        dirty2 = False
        for track, sector in set_of_ts:
            if track <= 40:
                self._release_block((track, sector), track - 1, (16, 6), bamblock1, (17, 6))
                dirty1 = True
            else:
                self._release_block((track, sector), track - 41, (16, 6), bamblock2, (17, 6))
                dirty2 = True
        if dirty1:
            self.write_ts(40, 1, bamblock1)
        if dirty2:
            self.write_ts(40, 2, bamblock2)

    def try_to_allocate(self, track, sector, exact):
        if track <= 40:
            ts = self._try_to_allocate((track, sector), track - 1, (16, 6), (40, 1), (17, 6), exact=exact)
        else:
            ts = self._try_to_allocate((track, sector), track - 41, (16, 6), (40, 2), (17, 6), exact=exact)
        return ts

################################################################################
# 1541 class

class _d41(d64):
    name = "1541"
    blocks_total = 683
    maxtrack = 35
    header_ts = (18, 0)
    header_name_offset = 144    # diskname and five-byte "pseudo id"
    directory_ts = (18, 1)
    std_max_dir_entries = 144   # for writing directory
    std_directory_interleave = 3
    std_file_interleave = 10

    def __init__(self):
        # build lookup tables for "sectors per track" and "blocks before track"
        self._sectors_of_track = {}
        self._blocks_before_track = {}
        lba = 0
        for track in range(1, self.maxtrack + 1):
            if track <= 17:
                sectors = 21
            elif track <= 24:
                sectors = 19
            elif track <= 30:
                sectors = 18
            else:
                sectors = 17
            self._sectors_of_track[track] = sectors
            self._blocks_before_track[track] = lba
            lba += sectors

    def sectors_of_track(self, track):
        self._check_track_num(track)
        return self._sectors_of_track[track]

    def ts_to_lba(self, track, sector):
        self._check_ts(track, sector)
        return self._blocks_before_track[track] + sector

    def bam_check(self):
        debug(1, "Checking 1541 BAM")
        self._check_totals(1, 35, (4, 4), (18, 0), 3, (5, 4))

    def read_free_blocks(self):
        d = dict()
        self._fill_free_blocks_dict(d, 1, 35, (4, 4), (18, 0))
        return d

    def release_blocks(self, set_of_ts):
        bamblock = self.read_ts(18, 0)
        dirty = False
        for track, sector in set_of_ts:
            self._release_block((track, sector), track - 1, (4, 4), bamblock, (5, 4))
            dirty = True
        if dirty:
            self.write_ts(18, 0, bamblock)

    def try_to_allocate(self, track, sector, exact):
        return self._try_to_allocate((track, sector), track - 1, (4, 4), (18, 0), (5, 4), exact=exact)

################################################################################
# 40-track-1541 class

class _d41_40(_d41):
    name = "40-track 1541"
    blocks_total = 683 + 5 * 17
    maxtrack = 40

    def try_to_allocate(self, track, sector, exact):
        if track <= 35:
            return super().try_to_allocate(track, sector, exact)   # just call 1541 method
        # AFAIK there are two different ways to store the additional 5*4 bytes,
        # either at the start or at the end of the "unused" part of t18s0.
        raise Exception("Allocation of tracks 36..40 is not yet supported!")

################################################################################
# 1571 class

class _d71(_d41):
    name = "1571"
    blocks_total = 2 * 683
    maxtrack = 70
    std_file_interleave = 6

    def __init__(self):
        super().__init__()  # let 1541 class build the lookup tables...
        # ...and now fix them for tracks 36..70:
        lba = self._blocks_before_track[36] # lba of "second half"
        for track in range(36, self.maxtrack + 1):
            sectors = self._sectors_of_track[track - 35]
            self._sectors_of_track[track] = sectors # use the same values for tracks 36..70 as for tracks 1..35
            self._blocks_before_track[track] = lba
            lba += sectors

    def bam_check(self):
        super().bam_check() # let 1541 class check the first side...
        # ...and now check second side:
        debug(1, "Checking 1571 BAM (second side)")
        self._check_totals(36, 35, (221, 1), (18, 0), 3, (0, 3), (53, 0))

    def read_free_blocks(self):
        d = super().read_free_blocks()  # let 1541 class do the first side...
        # and now add second side:
        self._fill_free_blocks_dict(d, 36, 35, (221, 1), (18, 0))
        # TODO - find a way to include "would show XYZ in a 1541 drive" info!
        return d

    def release_blocks(self, set_of_ts):
        bamblock180 = self.read_ts(18, 0)
        bamblock530 = self.read_ts(53, 0)
        dirty180 = False
        dirty530 = False
        for track, sector in set_of_ts:
            if track <= 35:
                self._release_block((track, sector), track - 1, (4, 4), bamblock180, (5, 4))
                dirty180 = True
            else:
                self._release_block((track, sector), track - 36, (221, 1), bamblock180, (0, 3), bamblock530)
                dirty180 = True
                dirty530 = True
        if dirty180:
            self.write_ts(18, 0, bamblock180)
        if dirty530:
            self.write_ts(53, 0, bamblock530)

    def try_to_allocate(self, track, sector, exact):
        if track <= 35:
            ts = self._try_to_allocate((track, sector), track - 1, (4, 4), (18, 0), (5, 4), exact=exact)    # all at t18s0
        else:
            ts = self._try_to_allocate((track, sector), track - 36, (221, 1), (18, 0), (0, 3), (53, 0), exact=exact)    # totals at t18s0, bitmaps at t53s0
        return ts

################################################################################
# wrapper stuff

# collect supported sizes and largest of them so files can be identified
_type_of_size = dict()
for imgtype in (_d41, _d41_40, _d71, _d81):
    _type_of_size[imgtype.blocks_total * 256] = (imgtype, False)    # no error info
    _type_of_size[imgtype.blocks_total * 257] = (imgtype, True) # with error info
_largest_size = max(_type_of_size.keys())
del imgtype # we do not want this to pop up in the online help...

def DiskImage(filename, writeback=False, writethrough=False):
    """
    Identify disk image file and return as d64 object.
    'writeback' and 'writethrough' are mutually exclusive.
    If neither is true, data is read-only.
    If 'writeback' mode is selected, you need to call writeback() to flush changes to file.
    """
    if writeback and writethrough:
        raise Exception("writeback and writethrough are mutually exclusive")
    readonly = not (writeback or writethrough)
    # open file and read data
    mode = "rb" if readonly else "r+b"
    fh = open(filename, mode)
    body = fh.read(_largest_size + 1)   # 1 more than largest supported type
    # FIXME - old version intercepted OSError! fix callers so they do that!
    if not readonly:
        body = bytearray(body)  # allow changes
    filesize = len(body)
    if filesize not in _type_of_size:
        raise Exception("Could not process " + filename + ": Image type not recognised")
    img_type, error_info = _type_of_size[filesize]
    obj = img_type()
    debug(1, filename, "is a", obj.name, "disk image" + (" with error info." if error_info else "."))
    obj._populate(fh=fh, body=body, readonly=readonly, writeback=writeback, writethrough=writethrough)
    return obj

################################################################################
# "main program"

petscii_graphics = (
    "        (control  codes)        " +
    " !\"#$%&'()*+,-./0123456789:;<=>?" +
    "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[£]↑←" +
    "━♠┃━⎻⎺⎼⎢⎥╮╰╯ᒪ╲╱ᒥᒣ●⎽♥▏╭╳○♣▕♦╋⡇┃π◥" +    # copy of c0..df, for output only
    "        (control  codes)        " +
    " ▌▄▔▁▏▒▕⣤◤▊┣▗┗┓▂┏┻┳┫▎▍▋▆▅▃ᒧ▖▝┛▘▚" +
    "━♠┃━⎻⎺⎼⎢⎥╮╰╯ᒪ╲╱ᒥᒣ●⎽♥▏╭╳○♣▕♦╋⡇┃π◥" +
    " ▌▄▔▁▏▒▕⣤◤▊┣▗┗┓▂┏┻┳┫▎▍▋▆▅▃ᒧ▖▝┛▘π")     # copy of a0..be, for output only
petscii_lowercase = (
    "        (control  codes)        " +
    " !\"#$%&'()*+,-./0123456789:;<=>?" +
    "@abcdefghijklmnopqrstuvwxyz[£]↑←" +
    "━ABCDEFGHIJKLMNOPQRSTUVWXYZ╋⡇┃▒▧" +    # copy of c0..df, for output only
    "        (control  codes)        " +
    " ▌▄▔▁▏▒▕⣤▨▊┣▗┗┓▂┏┻┳┫▎▍▋▆▅▃✓▖▝┛▘▚" +
    "━ABCDEFGHIJKLMNOPQRSTUVWXYZ╋⡇┃▒▧" +
    " ▌▄▔▁▏▒▕⣤▨▊┣▗┗┓▂┏┻┳┫▎▍▋▆▅▃✓▖▝┛▘▒")     # copy of a0..be, for output only
# CAUTION, these symbols need to be shown in reverse: "▊▋▆▅"

def from_petscii(bindata, second_charset):
    """Helper fn to convert petscii bytes to something printable."""
    charset = petscii_lowercase if second_charset else petscii_graphics
    ret = ""
    for b in bindata:
        # first check if we need to display in reverse
        revs = b in [0xaa, 0xb6, 0xb7, 0xb8, 0xea, 0xf6, 0xf7, 0xf8]
        # then check for control codes, which are shown as other symbols in reverse
        if (b & 0x60) == 0: # control code?
            b |= 64 # then replace
            revs = True # and reverse
        if revs:
            ret += "\033[7m"    # ANSI reverse
        ret += charset[b]
        if revs:
            ret += "\033[27m"   # ANSI revs off
    return ret

def filetypes(filetype):
    """
    Convert CBM file type to text representation.
    """
    filetype &= 15
    if filetype == 0:
        return b"DEL"
    if filetype == 1:
        return b"SEQ"
    if filetype == 2:
        return b"PRG"
    if filetype == 3:
        return b"USR"
    if filetype == 4:
        return b"REL"
    if filetype == 5:
        return b"CBM"   # introduced in 1581
    if filetype == 6:
        return b"DIR"   # introduced by CMD
    return b"0X%X" % filetype

def quote(name16):
    """
    Helper function to add opening and closing quotes at correct positions.
    """
    name = b"\xa0" + name16 + b"\xa0"   # add shift-spaces before and after, then...
    return name.replace(b"\xa0", b'"', 2)   # ...replacing the first two should do the trick.
    # FIXME - I think the result may differ from CBM DOS if the initial name actually contains quotes.

def show_directory(img, second_charset, full=False):
    """
    Show directory of image file.
    """
    name, id5 = img.read_header_fields()
    qname = quote(name)
    qname = from_petscii(qname, second_charset)
    id5 = from_petscii(id5, second_charset)
    print('header:', qname, id5)
    nonempty = 0
    empty = 0
    for entry in img.read_directory_entries(include_invisible=True):
        index = entry[0]
        bin30 = entry[1]
        filetype = bin30[0]
        if filetype:
            nonempty += 1
        else:
            empty += 1
            if full == False:
                continue
        ts = bin30[1:3]
        qname = quote(bin30[3:19])
        misc = bin30[19:28]
        blocks = bin30[28:30]
        line = '%3d: %5d ' % (index, int.from_bytes(blocks, "little"))
        line += from_petscii(qname, second_charset)
        line += " " if filetype & 128 else "*"
        line += from_petscii(filetypes(filetype), second_charset)
        line += "<" if filetype & 64 else " "
        if full:
            line += " : "
            line += ts.hex(" ")
            line += " : "
            line += misc.hex(" ")
        print(line)
    freeblocks = img.read_free_blocks()
    print(freeblocks["shown"], "blocks free (+%d in dir track)" % freeblocks["dir"])
    print("(%d directory entries, +%d empty)" % (nonempty, empty))

def process_file(file, second_charset, full=False):
    try:
        img = DiskImage(file)
    except OSError as e:
        print("Error: Could not open " + file + ": " + e.strerror, file=sys.stderr)
        return
    except Exception as e:
        print(e, "\n", file=sys.stderr)
        return
    img.bam_check()
    show_directory(img, second_charset, full=full)

def main():
    global debuglevel
    if len(sys.argv) == 1:
        sys.argv.append("-h")   # if run without arguments, show help instead of complaining!
    parser = argparse.ArgumentParser(allow_abbrev = False, description =
"""
This is a library for accessing d64 disk image files.
It works with 1541, 1571 and 1581 images.
If run directly, the directory of the given file(s) is displayed.
""")
    parser.add_argument("-c", "--charset", action="store_true", help="Use the other charset.")
    parser.add_argument("-d", "--debug", action="count", default=0, help="Increase debugging output.")
    parser.add_argument("-f", "--full", action="store_true", help="Show hidden data as well.")
    parser.add_argument("files", metavar="IMAGEFILE.D64", nargs='+', help="Disk image file.")
    args = parser.parse_args()
    debuglevel += args.debug
    for file in args.files:
        process_file(file, args.charset, full=args.full)

if __name__ == "__main__":
    main()
