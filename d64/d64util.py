#!/usr/bin/env python3
"""
for accessing commodore disc image files
"""
import argparse
import sys  # for sys.stderr and sys.exit

# this library accepts disk images with error info, but atm does not honor it!
# TODO: honor error info
# TODO: allow "ts tuple" and "lba number" to be used interchangeably

_debuglevel = 1

def _debug(minlevel, *unnamed, **named):
    """Helper function for debugging output."""
    if _debuglevel >= minlevel:
        print("debug%d:" % minlevel, *unnamed, **named)

def _popcount(integer):
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

# error codes
_errorcodes_map = {
    1: 1,   # ok
    2: 2,   # header prefix not found
    3: 3,   # sync not found
    4: 4,   # data prefix not found
    5: 5,   # data checksum error
    6: 6,   # decoding error
    9: 9,   # header checksum error
    11: 11, # id mismatch
    # mapping for wrong values (from status):
    0: 1,
    20: 2,
    21: 3,
    22: 4,
    23: 5,
    24: 6,
    27: 9,
    29: 11,
}
_errorcodes_chars = "X.2s4c6789hi"  # characters for display (X=illegal, .=ok)

################################################################################
# virtual base class

class d64(object):
    """
    This class describes a cbm disc image. There are subclasses for 1541,
        40track, 1571, 1581, ...

    Constants:
        name: the name of the format, for example "1581"
        blocks_total: the number of 256-byte blocks per image file
        maxtrack: highest track number (lowest is always 1)
        track_length_changes: dict with "new" sectors-per-track value
        header_ts_and_offset: where to find disc name and five-byte
            "pseudo id"
        bam_blocks: list of ts values where block availability map is held
        bam_start_size_maxperblock: start offset in bam blocks, size of
            entries, max number of entries per block
        directory_ts: track and sector of first directory block
        std_max_dir_entries: maximum number of directory entries
        std_directory_interleave: block interleave for directory
        std_file_interleave: block interleave for files
    """

    def __init__(self):
        raise Exception("Only subclasses (1541, 1571, 1581, ...) can be instantiated!")

    def _populate(self, fh, body, readonly, writeback, writethrough):
        """
        Called after correct subclass has been instantiated.
        """
        self.fh = fh
        self.body = body
        self.readonly = readonly    # we don't really need to store this info
        self.writebackmode = writeback  # in three vars, but it might make
        self.writethroughmode = writethrough    # the code more readable.
        self.need_writeback = False
        # build lookup tables for "sectors per track" and "blocks before track":
        self._sectors_of_track = {}
        self._blocks_before_track = {}
        lba = 0
        sectors = None  # trigger exception if track_length_changes has no "1" key!
        for track in range(1, self.maxtrack + 1):
            sectors = self.track_length_changes.get(track, sectors) # get new length or keep old one
            self._sectors_of_track[track] = sectors
            self._blocks_before_track[track] = lba
            lba += sectors
        # sanity check:
        if self.blocks_total != lba:
            sys.exit("BUG: total number of blocks is inconsistent!")
        # handle error codes separately
        self.errorblock = self.body[lba * 256:]
        self.body = self.body[:lba * 256]
        # display error block
        if _debuglevel >= 2:
            self._print_error_block()

    def _check_track_num(self, track):
        """_check_track_num(int)

        Throw exception if track number is invalid for this image type.
        """
        if track < 1:
            raise Exception("Given track number was lower than 1.")
        if track > self.maxtrack:
            raise Exception("Exceeded maximum track number.")

    def _check_ts(self, ts):
        """_check_ts((int, int))

        Throw exception if track/sector address is invalid for this image type.
        """
        track, sector = ts
        self._check_track_num(track)
        # now check sector number
        if sector < 0:
            raise Exception("Given sector number was negative.")
        if sector >= self.sectors_of_track(track):
            raise Exception("Exceeded maximum sector number of track %d." % track)

    def sectors_of_track(self, track):
        """sectors_of_track(int) -> int

        Return number of sectors in given track. Sector numbers start
        at 0, so the maximum sector number is one less than this.
        """
        self._check_track_num(track)
        return self._sectors_of_track[track]

    def ts_to_lba(self, ts):
        """ts_to_lba((int, int)) -> int

        Convert track and sector to logical block address
        (0-based block number).
        """
        self._check_ts(ts)
        track, sector = ts
        return self._blocks_before_track[track] + sector

    def read_ts(self, ts):
        """
        Read block indicated via track/sector tuple and return as bytes or
        bytearray.
        """
        _debug(2, "Reading t%d s%d." % ts)
        self._check_ts(ts)
        offset = 256 * self.ts_to_lba(ts)
        block = self.body[offset:offset+256]
        return block

    def write_ts(self, ts, data):
        """
        Write data to block indicated via track/sector tuple.
        """
        if len(data) != 256:
            raise Exception("block should be 256 bytes")
        _debug(2, "Writing t%d s%d." % ts)
        self._check_ts(ts)
        offset = 256 * self.ts_to_lba(ts)
        self.body[offset:offset+256] = data # will crash in readonly mode!
        if self.writebackmode:
            self.need_writeback = True
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
        if self.need_writeback:
            _debug(1, "Flushing to disk.")
            self.fh.seek(0)
            self.fh.write(self.body)
            self.fh.write(self.errorblock)
            self.fh.flush()
            self.need_writeback = False
        else:
            _debug(1, "Nothing changed, no need to flush to disk.")

    def _print_error_block(self):
        """
        Show contents of error block
        """
        if self.errorblock:
            print("Error block:")
            offset = 0
            for track in range(1, self.maxtrack + 1):
                out = "%3d: " % track
                for sector in range(self.sectors_of_track(track)):
                    code = self.errorblock[offset]
                    offset += 1
                    code = _errorcodes_map.get(code, 0) # 0 is invalid
                    char = _errorcodes_chars[code]
                    out += char
                print(out)
        else:
            print("Image does not have an error block.")

    def _check_totals(self, ts, first_byte_offset, size, howmanytracks, firsttrack):
        """
        Helper function to compare free blocks totals to free blocks bitmaps.

        This fn is not used for second side of 1571, where data is split over two blocks.

        ts: track and sector where data is stored
        first_byte_offset: offset of data in block
        size: bytes per entry
        howmanytracks: number of entries to process
        firsttrack: number to add to zero-based loop counter to get meaningful debug/error messages
        """
        block = self.read_ts(ts)
        for entry in range(howmanytracks):
            track = entry + firsttrack
            total = block[first_byte_offset + entry * size]
            sum = 0
            for i in range(size - 1):
                sum += _popcount(block[first_byte_offset + entry * size + 1 + i])
            if total == sum:
                _debug(9, "track %d: %d free blocks" % (track, sum))
            else:
                print("BAM error for track %d: counter says %d, bitfield says %d!" % (track, total, sum), file=sys.stderr)

    def _fill_free_blocks_dict(self, d, ts, offset, step, howmanytracks, firsttrack):
        """
        Helper function to read "free blocks" numbers from BAM.

        d: dictionary to put results in
        ts: track and sector where numbers are stored
        offset: offset of first entry in block
        step: step size to get to the next entry
        howmanytracks: number of entries to process
        firsttrack: number to add to zero-based loop counter to get track number

        The directory track is included in "all", but not in "shown"!
        """
        if "all" not in d:
            d["all"] = 0
        if "shown" not in d:
            d["shown"] = 0
        block = self.read_ts(ts)
        for entry in range(howmanytracks):
            track = entry + firsttrack
            freeblocks = block[offset + entry * step]
            d[track] = freeblocks
            d["all"] += freeblocks
            if track == self.directory_ts[0]:
                d["dir"] = freeblocks
            else:
                d["shown"] += freeblocks

    def bam_check(self, alt_maxtrack=None):
        """
        Check BAM for internal consistency: Do counters match bit fields?
        This does not check allocation of files!
        """
        _debug(1, "Checking " + self.name + " BAM for internal consistency")
        startoffset, size, maxtracks = self.bam_start_size_maxperblock
        starttrack = 1
        # if alternative maxtrack given, use it (used by 40track and 1571,
        # because only the 1541-compatible part of the bam uses the "standard"
        # layout)
        tracks_left = alt_maxtrack if alt_maxtrack else self.maxtrack
        for bamblock_ts in self.bam_blocks:
            self._check_totals(bamblock_ts, startoffset, size, maxtracks, starttrack)
            starttrack += maxtracks
            tracks_left -= maxtracks
            maxtracks = min(maxtracks, tracks_left)
        # sanity check:
        if tracks_left:
            sys.exit("BUG: inconsistent number of tracks when checking BAM!")

    def read_header_fields(self):
        """
        Return disk name and five-byte "pseudo id".
        """
        ts, of = self.header_ts_and_offset
        block = self.read_ts(ts)
        return block[of:of+16], block[of+18:of+23]

    def read_free_blocks(self, alt_maxtrack=None):
        """
        Return dictionary of free blocks per track.
        There are three additional keys, "shown", "dir" and "all", where
        "shown" + "dir" == "all"
        """
        _debug(3, "Reading free blocks counters")
        startoffset, size, maxtracks = self.bam_start_size_maxperblock
        starttrack = 1
        # if alternative maxtrack given, use it (used by 40track and 1571,
        # because only the 1541-compatible part of the bam uses the "standard"
        # layout)
        tracks_left = alt_maxtrack if alt_maxtrack else self.maxtrack
        d = dict()
        for bamblock_ts in self.bam_blocks:
            self._fill_free_blocks_dict(d, bamblock_ts, startoffset, size, maxtracks, starttrack)
            starttrack += maxtracks
            tracks_left -= maxtracks
            maxtracks = min(maxtracks, tracks_left)
        # sanity check:
        if tracks_left:
            sys.exit("BUG: inconsistent number of tracks when reading free blocks!")
        return d

    # TODO - simplify for non-1571 and give 1571 its own version
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
        totals_block = self.read_ts(totals_ts)
        bitmaps_offset, bitmaps_step = bitmaps_offset_step
        if bitmaps_ts == None:
            bitmaps_block = totals_block
        else:
            bitmaps_block = self.read_ts(bitmaps_ts)
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
                self.write_ts(totals_ts, totals_block)
                if bitmaps_block != totals_block:
                    self.write_ts(bitmaps_ts, bitmaps_block)
                _debug(2, "Allocated t/s", track, cand_sector)
                return track, cand_sector   # block has been allocated
            _debug(3, "t/s", track, cand_sector, "is not available")
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

    def _virtualfn(self):
        raise Exception("BUG: A virtual function was called!")

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
        ts = self.directory_ts
        all_used_ts = set() # for sanity check
        entry_number = 0    # index, so caller can unambiguously reference each entry
        while ts[0]:
            # sanity check
            if ts in all_used_ts:
                raise Exception("Directory loops back to itself, please check disk image!")
            all_used_ts.add(ts)
            try:
                self._check_ts(ts)
            except Exception as e:
                print("WARNING, stopped reading dir:", e)
                return
            # read a directory sector
            block = self.read_ts(ts)
            ts = (block[0], block[1])
            if ts[0] == 0:
                _debug(2, "Last! Link is t%d s%d." % ts)
            readidx = 2
            for i in range(8):
                bin30 = block[readidx:readidx+30]
                readidx += 32
                if bin30[0] or include_invisible:
                    yield entry_number, bin30   # CAUTION, maybe more fields will get added in future!
                entry_number += 1
        # "track" is zero, so there is no next block
        if ts[1] != 255:
            print("WARNING: sector value of final dir block is %d instead of 255!" % ts[1], file=sys.stderr)

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
            block = self.read_ts((track, sector))
            if block[0] == 0:
                _debug(2, "Last! Link is t%d s%d." % (block[0], block[1]))
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
                block[0:2] = self.get_new_block((track, sector))    # get a new block and let current block point to it
                init_link_ptrs = True   # make sure link pointers get overwritten with 00/ff from now on
            elif len(new_dir) == 0 and next_track:
                # we have more block(s) but no data for them
                self.free_block_chain((next_track, next_sector))    # free all following blocks
                block[0:2] = 0x00, 0xff # mark this block as final
            # write back
            self.write_ts((track, sector), block)
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
        entry = b"\x80" + bytes(self.header_ts_and_offset[0]) + name[:16] + bytes(11)
        # 0x80 is for a visible DEL entry,
        # tt/ss point to header block (so VALIDATE does not complain?),
        # name is 16 bytes of filename, padded with shift-space,
        # and then there are 11 bytes for various (DOS/GEOS) purposes, all 0 here
        return entry

    def get_new_block(self, prev_ts):
        """
        Find and allocate a new block after given block, return t/s.
        """
        prev_track, prev_sector = prev_ts
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

    def free_block_chain(self, ts):
        """
        Free all connected blocks, following the chain of link pointers.
        """
        track, sector = ts
        all_used_ts = set()
        while track:
            # sanity check
            if (track, sector) in all_used_ts:
                raise Exception("Block chain loops back to itself, please check disk image!")
            all_used_ts.add((track, sector))
            # read block, free it, go on with link
            block = self.read_ts((track, sector))
            track, sector = block[0:2]
        self.release_blocks(all_used_ts)

################################################################################
# CBM DOS 1 disk format, as used by 2040/3040 drives
# file extension is .d67 (at least in VICE)

class _dos1(d64):
    name = "2040/3040 (DOS 1.0)"
    blocks_total = 690  # 670 free
    maxtrack = 35
    track_length_changes = {1: 21, 18: 20, 25: 18, 31: 17}
    header_ts_and_offset = (18, 0), 144 # where to find diskname (dos 1 had no pseudo id?)
    bam_blocks = [(18, 0)]
    bam_start_size_maxperblock = (4, 4, 35)
    directory_ts = (18, 1)
    std_max_dir_entries = 152   # for writing directory
    #std_directory_interleave =
    #std_file_interleave =

    def __init__(self):
        # this inhibits the base class's exception...
        pass    # ...but there is nothing to do!

################################################################################
# CBM DOS 2 disk format, as used by 4040, 2031, 4031, 1540, 1541, 1551, 1570 drives
# file extension is mostly .d64, sometimes .d41

class _1541(d64):
    name = "1541"
    blocks_total = 683
    maxtrack = 35
    track_length_changes = {1: 21, 18: 19, 25: 18, 31: 17}
    header_ts_and_offset = (18, 0), 144 # where to find diskname and five-byte "pseudo id"
    bam_blocks = [(18, 0)]
    bam_start_size_maxperblock = (4, 4, 35)
    directory_ts = (18, 1)
    std_max_dir_entries = 144   # for writing directory
    std_directory_interleave = 3
    std_file_interleave = 10

    def __init__(self):
        # this inhibits the base class's exception...
        pass    # ...but there is nothing to do!

    def release_blocks(self, set_of_ts):
        bamblock = self.read_ts((18, 0))
        dirty = False
        for track, sector in set_of_ts:
            self._release_block((track, sector), track - 1, (4, 4), bamblock, (5, 4))
            dirty = True
        if dirty:
            self.write_ts((18, 0), bamblock)

    def try_to_allocate(self, track, sector, exact):
        return self._try_to_allocate((track, sector), track - 1, (4, 4), (18, 0), (5, 4), exact=exact)

################################################################################
# disk format of 8050 drives
# file extension is .d80

class _8050(d64):
    name = "8050"
    blocks_total = 2083 # 2052 free
    maxtrack = 77
    track_length_changes = {1: 29, 40: 27, 54: 25, 65: 23}
    header_ts_and_offset = (39, 0), 6   # where to find diskname and five-byte "pseudo id"
    bam_blocks = [(38, 0), (38, 3)] # one source says sector 1 instead of 3!
    bam_start_size_maxperblock = (6, 5, 50)
    directory_ts = (39, 1)
    std_max_dir_entries = 224   # for writing directory
    #std_directory_interleave = 3   ?
    #std_file_interleave = 10       ?

    def __init__(self):
        # this inhibits the base class's exception...
        pass    # ...but there is nothing to do!

    def release_blocks(self, set_of_ts):
        bamblock380 = self.read_ts((38, 0))
        bamblock383 = self.read_ts((38, 3))
        dirty380 = False
        dirty383 = False
        for track, sector in set_of_ts:
            if track <= 50:
                self._release_block((track, sector), track - 1, (6, 5), bamblock380, (7, 5))
                dirty380 = True
            else:
                self._release_block((track, sector), track - 51, (6, 5), bamblock383, (7, 5))
                dirty383 = True
        if dirty380:
            self.write_ts((38, 0), bamblock380)
        if dirty383:
            self.write_ts((38, 3), bamblock383)

    def try_to_allocate(self, track, sector, exact):
        if track <= 50:
            ts = self._try_to_allocate((track, sector), track - 1, (6, 5), (38, 0), (7, 5), exact=exact)
        else:
            ts = self._try_to_allocate((track, sector), track - 51, (6, 5), (38, 3), (7, 5), exact=exact)
        return ts

################################################################################
# disk format of 8250 / SFD-1001 drives
# file extension is .d82

class _8250(_8050):
    name = "8250"
    blocks_total = 4166 # 4133 free
    maxtrack = 154
    track_length_changes = {1: 29, 40: 27, 54: 25, 65: 23, 77+1: 29, 77+40: 27, 77+54: 25, 77+65: 23}
    bam_blocks = [(38, 0), (38, 3), (38, 6), (38, 9)]

    def __init__(self):
        # this inhibits the base class's exception...
        pass    # ...but there is nothing to do!

    def release_blocks(self, set_of_ts):
        bamblock380 = self.read_ts((38, 0))
        bamblock383 = self.read_ts((38, 3))
        bamblock386 = self.read_ts((38, 6))
        bamblock389 = self.read_ts((38, 9))
        dirty380 = False
        dirty383 = False
        dirty386 = False
        dirty389 = False
        for track, sector in set_of_ts:
            if track <= 50:
                self._release_block((track, sector), track - 1, (6, 5), bamblock380, (7, 5))
                dirty380 = True
            elif track <= 100:
                self._release_block((track, sector), track - 51, (6, 5), bamblock383, (7, 5))
                dirty383 = True
            elif track <= 150:
                self._release_block((track, sector), track - 101, (6, 5), bamblock386, (7, 5))
                dirty386 = True
            else:
                self._release_block((track, sector), track - 151, (6, 5), bamblock389, (7, 5))
                dirty389 = True
        if dirty380:
            self.write_ts((38, 0), bamblock380)
        if dirty383:
            self.write_ts((38, 3), bamblock383)
        if dirty386:
            self.write_ts((38, 6), bamblock386)
        if dirty389:
            self.write_ts((38, 9), bamblock389)

    def try_to_allocate(self, track, sector, exact):
        if track <= 50:
            ts = self._try_to_allocate((track, sector), track - 1, (6, 5), (38, 0), (7, 5), exact=exact)
        elif track <= 100:
            ts = self._try_to_allocate((track, sector), track - 51, (6, 5), (38, 3), (7, 5), exact=exact)
        elif track <= 150:
            ts = self._try_to_allocate((track, sector), track - 101, (6, 5), (38, 6), (7, 5), exact=exact)
        else:
            ts = self._try_to_allocate((track, sector), track - 151, (6, 5), (38, 9), (7, 5), exact=exact)
        return ts

################################################################################
# disk format of 1541-with-40-tracks-support
# file extension is mostly .d64, sometimes .d41

class _40track(_1541):
    name = "40-track 1541"
    blocks_total = 768  # 749 free
    maxtrack = 40

    def bam_check(self):
        super().bam_check(alt_maxtrack=35)  # let 1541 class check the first side...
        # ...and now check extra tracks:
        # TODO
        # AFAIK there are two different ways to store the additional 5*4 bytes,
        # either at the start or at the end of the "unused" part of t18s0.

    def read_free_blocks(self):
        d = super().read_free_blocks(alt_maxtrack=35)   # let 1541 class do the first side...
        # ...and now add extra tracks:
        # TODO
        # AFAIK there are two different ways to store the additional 5*4 bytes,
        # either at the start or at the end of the "unused" part of t18s0.
        #self._fill_free_blocks_dict(d, (18, 0), FIXME, 4, 5, 36)
        # TODO - find a way to include "would show XYZ in a 1541 drive" info!
        return d

    def try_to_allocate(self, track, sector, exact):
        if track <= 35:
            return super().try_to_allocate(track, sector, exact)   # just call 1541 method
        # AFAIK there are two different ways to store the additional 5*4 bytes,
        # either at the start or at the end of the "unused" part of t18s0.
        raise Exception("Allocation of tracks 36..40 is not yet supported!")

################################################################################
# disk format of 1571 drives
# file extension is .d71 or .d64

class _1571(_1541):
    name = "1571"
    blocks_total = 1366 # 1328 free
    maxtrack = 70
    track_length_changes = {1: 21, 18: 19, 25: 18, 31: 17, 35+1: 21, 35+18: 19, 35+25: 18, 35+31: 17}
    std_file_interleave = 6

    def bam_check(self):
        super().bam_check(alt_maxtrack=35)  # let 1541 class check the first side...
        # ...and now check second side:
        #_debug(1, "Checking 1571 BAM (second side)")
        # BAM for second side is split into two parts (totals at 18,0 and bitmaps at 53,0),
        # so we cannot use the std function:
        total_block = self.read_ts((18, 0))
        bits_block = self.read_ts((53, 0))
        for entry in range(35):
            track = entry + 36
            total = total_block[221 + entry]
            sum = 0
            for i in range(3):
                sum += _popcount(bits_block[entry * 3 + i])
            if total == sum:
                _debug(9, "track %d: %d free blocks" % (track, sum))
            else:
                print("BAM error for track %d: counter says %d, bitfield says %d!" % (track, total, sum), file=sys.stderr)

    def read_free_blocks(self):
        d = super().read_free_blocks(alt_maxtrack=35)   # let 1541 class do the first side...
        # ...and now add second side:
        self._fill_free_blocks_dict(d, (18, 0), 221, 1, 35, 36)
        # TODO - find a way to include "would show XYZ in a 1541 drive" info!
        return d

    def release_blocks(self, set_of_ts):
        bamblock180 = self.read_ts((18, 0))
        bamblock530 = self.read_ts((53, 0))
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
            self.write_ts((18, 0), bamblock180)
        if dirty530:
            self.write_ts((53, 0), bamblock530)

    def try_to_allocate(self, track, sector, exact):
        if track <= 35:
            ts = self._try_to_allocate((track, sector), track - 1, (4, 4), (18, 0), (5, 4), exact=exact)    # all at t18s0
        else:
            ts = self._try_to_allocate((track, sector), track - 36, (221, 1), (18, 0), (0, 3), (53, 0), exact=exact)    # totals at t18s0, bitmaps at t53s0
        return ts

################################################################################
# disk format of 1581 drives
# file extension is .d81 or .d64

class _1581(d64):
    name = "1581"
    blocks_total = 3200 # 3160 free
    maxtrack = 80
    track_length_changes = {1: 40}  # all tracks have 40 sectors
    header_ts_and_offset = (40, 0), 4   # where to find diskname and five-byte "pseudo id"
    bam_blocks = [(40, 1), (40, 2)]
    bam_start_size_maxperblock = (16, 6, 40)
    directory_ts = (40, 3)
    std_max_dir_entries = 296   # for writing directory
    std_directory_interleave = 1    # 1581 uses interleave 1 because of track cache
    std_file_interleave = 1 # 1581 uses interleave 1 because of track cache
    # HACK for img partition on test/demo disk:
    #header_ts_and_offset = (50, 0), 4
    #bam_blocks = [(50, 1), (50, 2)]
    #directory_ts = (50, 3)

    def __init__(self):
        # this inhibits the base class's exception...
        pass    # ...but there is nothing to do!

    def release_blocks(self, set_of_ts):
        bamblock1 = self.read_ts((40, 1))
        bamblock2 = self.read_ts((40, 2))
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
            self.write_ts((40, 1), bamblock1)
        if dirty2:
            self.write_ts((40, 2), bamblock2)

    def try_to_allocate(self, track, sector, exact):
        if track <= 40:
            ts = self._try_to_allocate((track, sector), track - 1, (16, 6), (40, 1), (17, 6), exact=exact)
        else:
            ts = self._try_to_allocate((track, sector), track - 41, (16, 6), (40, 2), (17, 6), exact=exact)
        return ts

################################################################################
# wrapper stuff

# collect supported sizes and largest of them so files can be identified
_supported = (_dos1, _1541, _40track, _1571, _1581, _8050, _8250)
_type_of_size = dict()
for imgtype in _supported:
    _type_of_size[imgtype.blocks_total * 256] = (imgtype, False)    # no error info
    _type_of_size[imgtype.blocks_total * 257] = (imgtype, True) # with error info
_largest_size = max(_type_of_size.keys())
del imgtype # we do not want this to pop up in the online help...

def list_formats():
    for type in _supported:
        print(type.name)

def DiskImage(filename, writeback=False, writethrough=False):
    """
    Identify disk image file and return as d64 object.
    'writeback' and 'writethrough' are mutually exclusive.
    If neither is true, data is read-only.
    If 'writeback' mode is selected, you need to call writeback() to flush
    changes to file.
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
    _debug(1, filename, "is a", obj.name, "disk image" + (" with error info." if error_info else "."))
    obj._populate(fh=fh, body=body, readonly=readonly, writeback=writeback, writethrough=writethrough)
    return obj

################################################################################
# "main program"

_petscii_graphics = (
    "        (control  codes)        " +
    " !\"#$%&'()*+,-./0123456789:;<=>?" +
    "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[£]↑←" +
    "━♠┃━⎻⎺⎼⎢⎥╮╰╯ᒪ╲╱ᒥᒣ●⎽♥▏╭╳○♣▕♦╋⡇┃π◥" +    # copy of c0..df, for output only
    "        (control  codes)        " +
    " ▌▄▔▁▏▒▕⣤◤▊┣▗┗┓▂┏┻┳┫▎▍▋▆▅▃ᒧ▖▝┛▘▚" +
    "━♠┃━⎻⎺⎼⎢⎥╮╰╯ᒪ╲╱ᒥᒣ●⎽♥▏╭╳○♣▕♦╋⡇┃π◥" +
    " ▌▄▔▁▏▒▕⣤◤▊┣▗┗┓▂┏┻┳┫▎▍▋▆▅▃ᒧ▖▝┛▘π")     # copy of a0..be, for output only
_petscii_lowercase = (
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
    charset = _petscii_lowercase if second_charset else _petscii_graphics
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

def _quote(name16):
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
    qname = _quote(name)    # FIXME - disk name uses a different quoting rule, end quote is *after* 16 chars!
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
        qname = _quote(bin30[3:19])
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

def _process_file(file, second_charset, full=False):
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

def _main():
    global _debuglevel
    if len(sys.argv) == 1:
        sys.argv.append("-h")   # if run without arguments, show help instead of complaining!
    parser = argparse.ArgumentParser(allow_abbrev = False, description =
"""
This is a library for accessing d64 disk image files.
If run directly, the directory of the given file(s) is displayed.
""")
    parser.add_argument("-c", "--charset", action="store_true", help="use the other charset")
    parser.add_argument("-d", "--debug", action="count", default=0, help="increase debugging output")
    parser.add_argument("-f", "--full", action="store_true", help="show hidden data as well")
    parser.add_argument("-l", "--list", action="store_true", help="list supported formats")
    parser.add_argument("files", metavar="IMAGEFILE.D64", nargs='+', help="Disk image file.")
    args = parser.parse_args()
    _debuglevel += args.debug
    if args.list:
        list_formats()
    for file in args.files:
        _process_file(file, args.charset, full=args.full)

if __name__ == "__main__":
    _main()

"""
dos 1.0:
    in 2040 and 3040 drives.
    tracks 18..24 have 20 sectors -> out of spec?
    690 blocks total, 670 free, 152 dir entries.
    no REL files.
dos 2.1:
    in 4040 drives and in upgraded 2040 and 3040 drives
    tracks 18..24 have 19 sectors
    683 blocks total, 664 free, 144 dir entries.
    new: REL files and @SAVE
dos 2.5
    in 8050 drives. disk changes are detected.
    2083 blocks total, 2052 free, 224 dir entries.
    REL files are still limited to 720 data blocks (180 kB), but these disks
    may have been been used in 8250 drives, therefore REL files on a dos 2.5
    disk may use the dos 2.7 format, i.e. with a super side sector - beware!
dos 2.6
    in 2031/4031 drives (and later in 1540, 1541, 1551, ...).
    basically dos 2.5, but
        downgraded to a single drive and
        using the same disk format as dos 2.1

dos 2.7
    in 8250 drives (double-sided)
    new: super side sector, so REL files can use the whole disk. but as this is
    incompatible with dos 2.5, support for sss can be switched on/off via m-w
    command. so image files of this type may hold REL file with or without sss!
dos 3.0
    in D9060 and D9090 hard disks

dos 3.0 (based on 2.6?)
    in 1570 and 1571 drives (but format on disk is still "2A")
    no super side sectors, so REL files are limited to 720 data blocks.
dos 3.1
    in 128 DCR (but format on disk is still "2A")

dos 10 (is this based on 2.7?)
    in 1581 drives (but format on disk is "3D")
    uses super side sector (-> up to 126 groups of six side sectors),
    so REL file can use whole disk.
"""
