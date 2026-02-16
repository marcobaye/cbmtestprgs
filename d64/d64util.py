#!/usr/bin/env python3
"""
A simple library for accessing Commodore disc image files.
"""
import argparse
import sys  # for sys.stderr and sys.exit
import enum

# this library accepts disk images with error info, but atm does not honor it!

_debuglevel = 1

def _debug(minlevel: int, *unnamed, **named):
    """Helper function for debugging output."""
    if _debuglevel >= minlevel:
        print("debug%d:" % minlevel, *unnamed, **named)

def _ceil_div(a: int, b: int) -> int:
    """ Divide and round up if remainder is nonzero """
    return -(a // -b)

def _popcount(integer: int) -> int:
    """Helper function to check bit fields in BAM."""
    return bin(integer).count("1")
    # since python 3.10 there is also:
    #   return integer.bit_count()
    # here's the old C algo:
    #counted = 0
    #while integer:
    #    counted += 1
    #    integer = integer & (integer - 1)
    #return counted

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

geos_file_types = {
    1:  b" BAS ",  # basic
    2:  b" ASS ",  # assembler
    3:  b" DAT ",  # data file
    4:  b" SYS ",  # system File
    5:  b" DAC ",  # desk accessory
    6:  b" APP ",  # application
    7:  b" DOC ",  # application data
    8:  b" FON ",  # font
    9:  b" PRN ",  # printer driver
    10: b" INP ",  # input driver
    11: b" DRV ",  # disk driver
    12: b" BOO ",  # system boot file
    13: b" TMP ",  # temporary
    14: b" AUT ",  # auto-execute file
}

class d64Exception(Exception):
    pass
# TODO: define exceptions for "LoopDetected" and "DoublyAllocated"

class cbmdosException(d64Exception):
    def __init__(self, code: int, msg: str, track: int, sector: int) -> None:
        self.code = code
        self.msg = msg
        self.track = track
        self.sector = sector
        self.status = f"{code}, {msg}, {track}, {sector}"
    def __str__(self):
        return self.status

class WriteProtectOnException(cbmdosException):
    def __init__(self):
        super().__init__(26, "write protect on", 0, 0)

class IllegalTrackOrSectorException(cbmdosException):
    def __init__(self, track: int, sector: int):
        super().__init__(66, "illegal track or sector", track, sector)

class DiskFullException(cbmdosException):
    def __init__(self):
        super().__init__(72, "disk full", 0, 0)


################################################################################
# backend

class ImgMode(enum.Enum):
    """
    This class defines enum constants to set the "caching mode" of an image file:
        READONLY:
            any attempt to write raises an Exception.
        WRITEBACK:
            changes are in RAM only, calling writeback() flushes them to file.
        DETACHED:
            changes are in RAM only, calling writeback() raises an Exception.
        DRYRUN:
            changes are in RAM only, calling writeback() does nothing.
        WRITETHROUGH:
            all changes are flushed to file immediately.
    """
    READONLY = enum.auto()
    WRITEBACK = enum.auto()
    WRITETHROUGH = enum.auto()
    DETACHED = enum.auto()
    DRYRUN = enum.auto()

class imagefile(object):
    """
    This class acts as a common back-end for all the different formats. It only
    handles reading blocks from the actual file and writing them back.
    """
    def __init__(self, fh, body: bytearray, img_mode, has_error_chunk=False):
        self.fh = fh
        self.body = body
        if type(img_mode) == ImgMode:
            self.img_mode = img_mode
        else:
            raise Exception("Invalid ImgMode!")
        self.need_writeback = False
        # handle error codes separately
        if has_error_chunk:
            self.block_count, mustbezero = divmod(len(body), 257)
        else:
            self.block_count, mustbezero = divmod(len(body), 256)
        if mustbezero:
            sys.exit("BUG: illegal file size!")
        self.error_chunk = self.body[self.block_count * 256:]
        self.body = self.body[:self.block_count * 256]
        # for accessing partitions:
        self.partition_start = 0
        self.partition_size = self.block_count

    def partition_set(self, start_and_size: (int, int)):
        """
        Set start and size of fs in file
        (so user can access a 1541 partition inside a .d4m file)
        """
        self.partition_start, self.partition_size = start_and_size

    def has_error_chunk(self) -> bool:
        return bool(len(self.error_chunk))

    def errorchunk_read_entry(self, lba: int) -> int:
        return self.error_chunk[self.partition_start + lba]

    def block_read(self, lba: int) -> bytes or bytearray:
        offset = 256 * (self.partition_start + lba)
        block = self.body[offset:offset+256]
        return block

    def block_write(self, lba: int, data: bytes) -> None:
        offset = 256 * (self.partition_start + lba)
        self.body[offset:offset+256] = data
        if self.img_mode == ImgMode.READONLY:
            raise WriteProtectOnException() # "Tried to write in readonly mode."
        elif self.img_mode == ImgMode.WRITEBACK:
            self.need_writeback = True
        elif self.img_mode == ImgMode.DETACHED:
            pass
        elif self.img_mode == ImgMode.DRYRUN:
            pass
        elif self.img_mode == ImgMode.WRITETHROUGH:
            self.fh.seek(offset)
            self.fh.write(data)
            self.fh.flush()
        else:
            sys.exit("BUG: invalid image mode")

    def writeback(self):
        """
        If data has been changed, write to file.
        """
        if self.img_mode == ImgMode.READONLY:
            raise Exception("Called writeback() in readonly mode.")
        elif self.img_mode == ImgMode.WRITEBACK:
            if self.need_writeback:
                _debug(1, "Flushing to disk.")
                self.fh.seek(0)
                self.fh.write(self.body)
                self.fh.write(self.error_chunk)
                self.fh.flush()
                self.need_writeback = False
            else:
                _debug(1, "Nothing changed, no need to flush to disk.")
        elif self.img_mode == ImgMode.DETACHED:
            raise Exception("Called writeback() in detached mode.")
        elif self.img_mode == ImgMode.DRYRUN:
            pass
        elif self.img_mode == ImgMode.WRITETHROUGH:
            raise Exception("Called writeback() in writethrough mode.")
        else:
            sys.exit("BUG: invalid image mode")

################################################################################
# virtual base class

class d64(object):
    """
    This class describes a cbm disc image. There are subclasses for 1541,
    40track, 1571, 1581, ...
    """
    # Constants:
    #    name: the name of the format, for example "1581"
    #    blocks_total: the number of 256-byte blocks per image file
    #    mintrack: lowest track number (1)
    #    maxtrack: highest track number
    #    track_length_changes: dict with "new" sectors-per-track value
    #    header_ts: where to find disc name and five-byte "pseudo id"
    #    header_offset: byte offset in block
    #    bam_blocks: list of ts values where block availability map is held
    #    bam_start_size_maxperblock: start offset in bam blocks, size of
    #        entries, max number of entries per block
    #    directory_ts: track and sector of first directory block
    #    std_max_dir_entries: maximum number of directory entries
    #    std_directory_interleave: block interleave for directory
    #    std_file_interleave: block interleave for files
    #    filetypes: dict of supported file types
    #    use_super_side_sector: REL files can be larger than 720+6 blocks
    mintrack = 1    # only D9060/D9090 differ, they use 0
    # used for optional "interleave sides" optimization:
    allow_1571_8250_trick = False   # this flag says whether the var below can be changed
    sides_interleave    = 0 # track diff between sides, zero means "do not interleave sides"

    def __init__(self):
        if "blocks_total" not in self.__dir__():
            raise Exception("Only subclasses (1541, 1571, 1581, ...) can be instantiated!")

    # TODO: fix code so all this stuff can be moved to __init__()!
    def _populate(self, imagefile):
        """
        Called after correct subclass has been instantiated.
        """
        self.imagefile = imagefile
        # build lookup tables for "sectors per track" and "blocks before track":
        self._sectors_of_track = {}
        self._blocks_before_track = {}
        lba = 0
        sectors = None  # trigger exception if track_length_changes has no key for mintrack!
        for track in range(self.mintrack, self.maxtrack + 1):
            sectors = self.track_length_changes.get(track, sectors) # get new length or keep old one
            self._sectors_of_track[track] = sectors
            self._blocks_before_track[track] = lba
            lba += sectors
        # sanity check:
        if self.blocks_total != lba:
            sys.exit("BUG: total number of blocks is inconsistent!")
        # display error chunk
        if _debuglevel >= 2:
            self.errorchunk_display()

    def _check_track_num(self, ts: (int, int)) -> None:
        """_check_track_num((int, int))

        Throw exception if track number is invalid for this image type.
        """
        track, sector = ts
        if track < self.mintrack:
            raise IllegalTrackOrSectorException(track, sector)  # "Track number too low."
        if track > self.maxtrack:
            raise IllegalTrackOrSectorException(track, sector)  # "%d exceeds maximum track number (%d)." % (track, self.maxtrack)

    def _check_ts(self, ts: (int, int)) -> None:
        """_check_ts((int, int))

        Throw exception if track/sector address is invalid for this image type.
        """
        self._check_track_num(ts)
        track, sector = ts
        # now check sector number
        if sector < 0:
            raise IllegalTrackOrSectorException(track, sector)  # "Given sector number was negative."
        maxsector = self.sectors_of_track(track) - 1
        if sector > maxsector:
            raise IllegalTrackOrSectorException(track, sector)  # "%d exceeds maximum sector number of track %d (%d)." % (sector, track, maxsector)

    def sectors_of_track(self, track: int) -> int:
        """sectors_of_track(int) -> int

        Return number of sectors in given track. Sector numbers start
        at 0, so the maximum sector number is one less than this.
        """
        self._check_track_num((track, 0))
        return self._sectors_of_track[track]

    def ts_to_lba(self, ts: (int, int)) -> int:
        """ts_to_lba((int, int)) -> int

        Convert track and sector to logical block address
        (0-based block number).
        """
        self._check_ts(ts)
        track, sector = ts
        return self._blocks_before_track[track] + sector

    def block_read(self, ts: (int, int)) -> bytes or bytearray:
        """
        Read block indicated via track/sector tuple and return as bytes or
        bytearray.
        """
        _debug(2, "Reading t%d s%d." % (ts[0], ts[1]))
        self._check_ts(ts)
        return self.imagefile.block_read(self.ts_to_lba(ts))

    def block_write(self, ts, data: bytes) -> None:
        """
        Write data to block indicated via track/sector tuple.
        """
        if len(data) != 256:
            raise Exception("block should be 256 bytes")
        _debug(2, "Writing t%d s%d." % ts)
        self._check_ts(ts)
        self.imagefile.block_write(self.ts_to_lba(ts), data)

    # TODO: add a new class for easier handling of bam blocks:
    # when reading block from t/s, remember where data came from to ease
    # writing it back:
    #def dosblock_read(self, ts):
    #    data = self.block_read(ts)
    #    return DOSblock(ts, data, self)

    def writeback(self) -> None:
        """
        If data has been changed, write to file.
        """
        self.imagefile.writeback()

    def errorchunk_display(self) -> None:
        """
        Show contents of error chunk
        """
        if self.imagefile.has_error_chunk():
            print("Error chunk:")
            offset = 0
            for track in range(self.mintrack, self.maxtrack + 1):
                out = "%3d: " % track
                for sector in range(self.sectors_of_track(track)):
                    code = self.imagefile.errorchunk_read_entry(offset)
                    offset += 1
                    code = _errorcodes_map.get(code, 0) # 0 is invalid
                    char = _errorcodes_chars[code]
                    out += char
                print(out)
        else:
            print("Image does not have an error chunk.")

    def _check_totals(self, ts: (int, int), first_byte_offset: int, size: int, howmanytracks: int, firsttrack: int) -> None:
        """
        Helper function to compare free blocks totals to free blocks bitmaps.

        This fn is not used for second side of 1571, where data is split over two blocks.

        ts: track and sector where data is stored
        first_byte_offset: offset of data in block
        size: bytes per entry
        howmanytracks: number of entries to process
        firsttrack: number to add to zero-based loop counter to get meaningful debug/error messages
        """
        block = self.block_read(ts)
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

    def _bam_fill_free_blocks_dict(self, d: dict, ts: (int, int), offset: int, step: int, howmanytracks: int, firsttrack: int) -> None:
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
        block = self.block_read(ts)
        for entry in range(howmanytracks):
            track = entry + firsttrack
            freeblocks = block[offset + entry * step]
            d[track] = freeblocks
            d["all"] += freeblocks
            if track == self.directory_ts[0]:
                d["reserved"] = freeblocks
            else:
                d["shown"] += freeblocks

    # FIXME: define *exactly* what this fn is for!
    # it was initially written as an easy sanity check to see if the BAM was found.
    # atm it just prints errors to stdout, so it's just info for the user - but
    #   then it should be explicitly requested/disabled by cli arg.
    # it _could_ set a flag for "BAM is not reliable" which then would have to
    #   be checked by other (write) functions, but then this fn should be called
    #   automatically when opening the image file.
    # it _could_ switch the image to read-only mode, but that's bad because not
    #   all writes need a valid BAM.
    # ...and if this tool ever gets a "validate" fn, it needs to be fine-grained
    #   about what the user actually wants to be fixed.
    def check_bam_counters(self, alt_maxtrack=None) -> None:
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

    def bam_read_header_fields(self) -> (int, bytes, bytes):
        """
        Return "drive", disk name and five-byte "pseudo id".
        """
        block = self.block_read(self.header_ts)
        of = self.header_offset
        return 0, block[of:of+16], block[of+18:of+23]

    def bam_read_free_blocks(self, alt_maxtrack=None) -> dict:
        """
        Return dictionary of free blocks per track.
        There are three additional keys, "shown", "reserved" and "all", where
        "shown" + "reserved" == "all"
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
            self._bam_fill_free_blocks_dict(d, bamblock_ts, startoffset, size, maxtracks, starttrack)
            starttrack += maxtracks
            tracks_left -= maxtracks
            maxtracks = min(maxtracks, tracks_left)
        # sanity check:
        if tracks_left:
            sys.exit("BUG: inconsistent number of tracks when reading free blocks!")
        return d

    def _bam_offset_and_bit(self, sector: int) -> (int, int):
        """
        Convert sector number to byte offset and bit value for accessing BAM bitmap.
        """
        # bitmap bytes are little-endian, lsb first:
        return sector >> 3, 1 << (sector & 7)
        # ...except for CMD native partitions, which is why they have their own version of this fn.

    def _bam_bitmap_to_string(self, bitmap: bytes, length: int) -> str:
        """
        Convert first "length" bits of BAM bitmap to string of zeroes and ones.
        """
        s = ""
        for byte in bitmap:
            byte |= 256
            while byte != 1:
                s += "." if byte & 1 else "#"
                byte >>= 1
                length -= 1
                if length == 0:
                    return s
        raise Exception("BUG in _bam_bitmap_to_string()")

    def _bamblock_release_block(self, ts: (int, int), bam_block: bytes, first_track: int, bam_offset: int, size: int) -> None:
        """
        Helper function to release a single block in BAM. This is used by 1541,
        8050, 8250, 1581 and first side of 1571. Other formats and second side
        of 1571 differ too much and therefore have their own functions.

        ts: track and sector of block to release
        bam_block: block with "totals/bitmaps" structs
        first_track: lowest track number described by this bam block
        bam_offset: offset in block where structs begin
        size: number of bytes per track struct (one "totals" byte and several "bitmaps" bytes)
        """
        # process args
        track, sector = ts
        # calculate offsets
        totals_offset = bam_offset + (track - first_track) * size
        byte_offset, bit_value = self._bam_offset_and_bit(sector)
        bitmap_offset = totals_offset + 1 + byte_offset # bitmaps start after totals -> +1
        # check/release:
        if bam_block[bitmap_offset] & bit_value:
            raise Exception("Attempted to free a block (t%ds%d) that is already free." % (track, sector))
        else:
            bam_block[bitmap_offset] |= bit_value
            bam_block[totals_offset] += 1
        # TODO - compare totals to maximum for this track?
        #raise Exception("BAM is corrupt, totals do not match bitmap.")

    def _bamblock_allocate_block(self, wanted_ts: (int, int), bam_ts: (int, int), table_offset: int, step: int, entry: int, exact=True) -> (int, int) or None:
        """
        Helper function to allocate a single block on given track in BAM. This
        is used by 1541, 8050, 8250, 1581 and first side of 1571. Other formats
        and second side of 1571 differ too much and therefore have their own
        functions.
        If block is available, allocate it and return t/s.
        If block is not available, return None.

        wanted_ts: track and sector to allocate
        bam_ts: track and sector where totals/bitmaps are stored
        table_offset: where totals and bitmaps start in bam block
        step: size of each entry
        entry: which entry to process
        exact: if False, function may allocate and return a different sector from this track
        """
        track, wanted_sector = wanted_ts
        bam_block = self.block_read(bam_ts)
        # calculate offsets
        totals_offset = table_offset + entry * step
        cand_sector = wanted_sector # we start the search with the wanted sector...
        num_sectors = self.sectors_of_track(track)  # ...and this is where we wrap around
        while True:
            byte_offset, bit_value = self._bam_offset_and_bit(cand_sector)
            bitmap_byte_offset = totals_offset + 1 + byte_offset
            # bit clear: sector is not available, i.e. is allocated or does not even exist
            # bit set: sector is available, i.e. can be allocated
            available = bool(bam_block[bitmap_byte_offset] & bit_value)
            if available:
                # allocate block
                bam_block[bitmap_byte_offset] &= (255 - bit_value)
                if bam_block[totals_offset]:
                    bam_block[totals_offset] -= 1
                else:
                    raise Exception("BAM is corrupt, totals do not match bitmap.")
                # write bam block
                self.block_write(bam_ts, bam_block)
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

    def _bam_to_list(self, ts: (int, int), offset: int, table_step: int, track1st: int, limit: int) -> list:
        """
        convert BAM bitmaps to human-readable version by returning list of
        tuples with track number and string.
        """
        ret = []
        bamblock = self.block_read(ts)
        for track in range(track1st, limit):
            bitmap = bamblock[offset:offset + table_step]
            ret.append((track, self._bam_bitmap_to_string(bitmap, self.sectors_of_track(track))))
            offset += table_step
        return ret

    def _bam_print_list(self, li: list) -> None:
        """
        print list generated by function above.
        """
        for track, string in li:
            print("t%03d:" % track, string)

    def _bam_print_transposed(self, li: list) -> None:
        """
        print list, but transpose to better use console
        """
        maxtrack = 0
        maxlength = 0
        fit = {}
        for track, string in li:
            maxtrack = max(maxtrack, track)
            maxlength = max(maxlength, len(string))
        for track, string in li:
            fit[track] = string.ljust(maxlength, " ")
        # TODO: show track numbers!
        for se in range(maxlength):
            print("s%02d: " % se, end="")
            for tr in range(1, 1+maxtrack):
                print(fit[tr][se], end="")
            print()

    def bam_display(self) -> None:
        print("Sorry, display of BAM has not yet been implemented for this image type.")

    def _virtualfn(self) -> None:
        raise Exception("BUG: A virtual function was called!")

    def bam_allocate_block(self, track: int, sector: int, exact: bool) -> (int, int) or None:
        """
        Helper function to allocate a single block in BAM.
        If block is available, allocate it and return t/s.
        If block is not available, return None.
        exact: if False, function may allocate and return a different sector from this track.
        """
        self._virtualfn()

    def bam_release_blocks(self, set_of_ts: list) -> None:
        """
        Free all blocks given as t/s tuples.
        """
        self._virtualfn()

    def linkptrs_follow(self, ts: (int, int), display=False, include_blocks=True) -> ((int, int), bytes or bytearray):
        """
        Follow link pointers and return each block as (ts, data) tuple.
        If include_blocks is False, only return ts tuples.
        """
        last_track = -1 # for display
        all_used_ts = set() # for sanity check
        while True:
            if ts[0] == 0:  # TODO: this might fail for DOS 7 where 0xff marks the end?
                if display:
                    print()
                break
            if display:
                if ts[0] != last_track:
                    print()
                    print("t %02x s" % ts[0], end="")
                    last_track = ts[0]
                print(" %02x" % ts[1], end="")
            # sanity checks:
            if ts in all_used_ts:
                if display:
                    print()
                raise Exception("Block chain loops back to itself, please check disk image!")
                # TODO: make this into some "d64.loop_exception" so caller can
                # check for it, or add some return_None_in_case_of_loop arg!
            try:
                self._check_ts(ts)
            except Exception as e:
                if display:
                    print()
                print("WARNING, stopped reading t/s chain:", e)
                return
                # TODO: make this into some "d64.illegal_track_or_sector_exception"
                # otherwise caller might try to interpret bogus data!
            all_used_ts.add(ts)
            # read block, go on with link
            block = self.block_read(ts)
            # read link before delivering block, because caller might alter it!
            next_ts = block[0], block[1]
            yield (ts, block) if include_blocks else ts
            ts = next_ts    # go on

    def geos_get_border_block_ts(self) -> (int, int) or None:
        """
        Return t/s of GEOS border block or None if there isn't one.
        """
        # FIXME: check some "supported_by_GEOS" var so this only works for 1541
        # (and 40t/1571 friends), 1581 and CMD native.
        block = self.block_read(self.header_ts)
        if block[0xad:0xba] == b"GEOS format V":    # goes on with "1.0" or "1.1"...
            return block[0xab], block[0xac]
        return None

    def _directory_iter(self, start_ts: (int, int)):
        """
        Return directory entries as (30-byte sequence, wherefound) tuples.
        """
        for ts, block in self.linkptrs_follow(start_ts):
            if block[0] == 0:
                _debug(2, "Last! Link is t%d s%d." % (block[0], block[1]))
                if block[1] != 255:
                   print("WARNING: sector value of final dir block is %d instead of 255!" % block[1], file=sys.stderr)
            readidx = 2
            for i in range(8):
                yield block[readidx:readidx + 30], (ts, readidx)
                readidx += 32

    def directory_read_entries(self, include_invisible=False):
        """
        Return directory entries as tuples:
        (index, raw entry (30 bytes), location, tuple might grow in future...)
        Setting include_invisible to True yields even the empty entries.
        """
        entry_number = 0    # index, so caller can unambiguously reference each entry
        # process directory
        for bin30, location in self._directory_iter(self.directory_ts):
            if bin30[0] or include_invisible:
                yield entry_number, bin30, location # CAUTION, maybe more fields will get added in future!
            entry_number += 1
        # process GEOS "border block", if there is one:
        borderblock_ts = self.geos_get_border_block_ts()
        if borderblock_ts:
            _debug(1, "Contents of GEOS border block:")
            for bin30, location in self._directory_iter(borderblock_ts):
                if bin30[0] or include_invisible:
                    yield entry_number, bin30, location # CAUTION, maybe more fields will get added in future!
                entry_number += 1

    def direntry_get(self, which: int) -> bytes or bytearray:
        """
        Return a single dir entry, specified as index.
        """
        for entry in self.directory_read_entries(include_invisible=True):
            if entry[0] == which:
                return entry[1] # only return bin30 data
            if entry[0] > which:
                sys.exit("BUG: This cannot happen!")
        sys.exit("Error: Specified directory index does not exist.")

    def direntry_cook(self, index: int, bin30: bytes or bytearray) -> (bool, int, str, str, str):
        """
        Convert directory entry to what is shown to user
        """
        # read fields
        cbm_type = bin30[0]
        t_s = bin30[1:3]
        file_name = bin30[3:19]
        second_t_s = bin30[19:21]   # t/s of side sector (REL) or info block (GEOS)
        record_length = bin30[21:22]    # GEOS: 0=normal file, 1=VLIR file
        geos_type = bin30[22:23]
        year = bin30[23]    # year
        month_day = bin30[24:26]    # month, day
        time = bin30[26:28] # hour, minute
        block_count = bin30[28:30]
        # process
        in_use = cbm_type != 0
        if ((cbm_type & 15) == 3) and geos_type[0] != 0:
            file_type = geos_file_types.get(geos_type[0], b" ??? ")
        else:
            file_type = self.filetype(cbm_type)
        if year >= 80:
            year += 1900    # one source claims "years since 1900"
        else:
            year += 2000    # but "wheels" seems to use 2 for 2002
        block_count = int.from_bytes(block_count, "little")
        optional = " : " + t_s.hex(" ") + " : " + second_t_s.hex(" ") + " " + record_length.hex(" ") + " " + geos_type.hex(" ") + " "
        if all(month_day):
            optional += ": %04d-%02d-%02d %02d:%02d" % (year, month_day[0], month_day[1], time[0], time[1])
        else:
            optional += bin30[23:28].hex(" ")
        return in_use, block_count, file_name, file_type, optional

    def directory_write(self, new_dir: list) -> None:   # TODO: add flag for "accept oversized dirs and use other tracks"
        """
        Overwrite directory with data given as list of 30-byte entries.
        """
        if len(new_dir) > self.std_max_dir_entries:
            raise DiskFullException()   # "New dir has too many entries (%d > %d)." % (len(new_dir), self.std_max_dir_entries)
        ts = self.directory_ts
        all_used_ts = set() # for sanity check
        init_link_ptrs = False  # flag needed when growing directory (to init fresh block)
        while ts[0]:
            # sanity check
            if ts in all_used_ts:
                raise Exception("Directory loops back to itself, please check disk image!")
            all_used_ts.add(ts)
            # read a directory sector
            block = self.block_read(ts)
            if init_link_ptrs:
                block[0:2] = 0x00, 0xff # init fresh dir block
            elif block[0] == 0:
                _debug(2, "Last! Link is t%d s%d." % (block[0], block[1]))
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
            next_ts = block[0], block[1]
            # now for the interesting part, enlarging/shrinking directory:
            if len(new_dir) and next_ts[0] == 0:
                # we have more data but no block to put it
                next_ts = self.bam_get_block_on_track(ts, self.std_directory_interleave)    # get a new block and let current block point to it
                if next_ts == None:
                    raise Exception("Directory track is full!") # TODO: raise DiskFullException
                block[0:2] = next_ts
                init_link_ptrs = True   # make sure link pointers get overwritten with 00/ff from now on
            elif len(new_dir) == 0 and next_ts[0]:
                # we have more block(s) but no data for them
                self.linkptrs_release_all(next_ts)  # free all following blocks
                block[0:2] = 0x00, 0xff # mark this block as the last
            # write back
            self.block_write(ts, block)
            # get (potentially modified) link ptr
            ts = block[0], block[1]
        # "track" is zero, so there is no next block
        if ts[1] != 255:
            print("WARNING: sector value of final dir block is %d instead of 255!" % ts[1], file=sys.stderr)
        # sanity check
        if len(new_dir):
            sys.exit("BUG: still data left after writing dir!")

    def direntry_make_dummy(self, namebytes=b"") -> bytes:
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
        # and then there are 11 bytes for various purposes, all 0 here:
        #   2 bytes for number of blocks
        #   3 bytes for REL files (or info block in case of GEOS)
        #   1 byte for GEOS file type
        #   5 bytes for year/month/day/hour/minute (GEOS/CMD date stamp)
        return entry

    def _bam_get_new_start_block(self) -> (int, int):
        """
        Find and allocate a block for a new file.
        Return ts or raise DiskFullException.
        """
        maxtrack = self.maxtrack
        if self.sides_interleave:
            maxtrack -= self.sides_interleave
            # optimization for 8250: first check dir track on second side
            # (won't work on 1571 because track 53 is completely allocated):
            ts = self.bam_allocate_block(self.directory_ts[0] + self.sides_interleave, 0, exact=False)
            if ts is not None:
                return ts
        # start next to directory track, alternating between lower and higher tracks:
        lower = self.directory_ts[0] - 1
        higher = self.directory_ts[0] + 1
        while True:
            checks = 0
            # check below directory:
# FIXME - 9060/9090 have a track 0, but it is not used for allocations, so this is wrong:
            if lower >= self.mintrack:
                ts = self.bam_allocate_block(lower, 0, exact=False)
                if ts is not None:
                    return ts
                if self.sides_interleave:
                    # optimization for 1571/8250: check corresponding track on other side:
                    ts = self.bam_allocate_block(lower + self.sides_interleave, 0, exact=False)
                    if ts is not None:
                        return ts
                lower -= 1
                checks += 1
            # check above directory:
            if higher <= maxtrack:
                ts = self.bam_allocate_block(higher, 0, exact=False)
                if ts is not None:
                    return ts
                if self.sides_interleave and higher <= self.sides_interleave:
                    # optimization for 1571/8250: check corresponding track on other side:
                    ts = self.bam_allocate_block(higher + self.sides_interleave, 0, exact=False)
                    if ts is not None:
                        return ts
                higher += 1
                checks += 1
            # if tracks are exhausted in both directions, give up:
            if checks == 0:
                break
        raise DiskFullException()

    def bam_get_block_on_track(self, prev_ts: (int, int) or None, interleave=None, cand_track=None) -> (int, int) or None:
        """
        Find and allocate a new block after given block on specific track.
        Return ts or None if track is full.
        """
        prev_track, prev_sector = prev_ts
        if interleave is None:
            interleave = self.std_file_interleave
        if cand_track is None:
            cand_track = prev_track
        # to get a candidate sector, add interleave:
        cand_sector = prev_sector + interleave
        # and wrap around if out of range:
        num_sectors = self.sectors_of_track(cand_track)
        # some clown may try to use an interleave of "100", so use "while" instead of "if":
        while cand_sector >= num_sectors:
            cand_sector -= num_sectors
        # just do what the DOS does, allocate the next possible sector:
        return self.bam_allocate_block(cand_track, cand_sector, exact=False)

    def bam_get_additional_block(self, prev_ts: (int, int) or None, interleave=None) -> (int, int):
        """
        Find and allocate a new block after given t/s tuple.
        If given t/s is None, find a block near directory (as CBM DOS does it).
        Return ts or raise DiskFullException.
        """
# CBM DOS algorithm when searching for free blocks:
# step 1: scan "away from directory", i.e. go down if lower and go up if higher.
# step 2: when reaching track limit, reverse direction and scan on other side of directory.
# step 3: when reaching track limit, reverse again and scan until reaching original track.
# i.e.:
#            dirtrack
#t1             |             tmax  (when starting below dir)
#tmax           |              t1   (when starting above dir)
#|              |               |
#<=======........................   step 1, scan from start pos to near edge
#................===============>   step 2, scan from dir to far edge
#........<======.................   step 3, scan from dir to start pos
        # if this is the first block, find one near directory:
        if prev_ts is None:
            return self._bam_get_new_start_block()
        # if optional "interleave sides" optimization for 1571/8250 is in
        # effect, adjust values:
        maxtrack = self.maxtrack
        original_track = prev_ts[0]
        if self.sides_interleave:
            maxtrack -= self.sides_interleave
            if original_track > maxtrack:
                original_track -= self.sides_interleave # convert "track" to "cylinder"
        # scan current track:
        ts = self.bam_get_block_on_track(prev_ts, interleave, original_track)
        if ts is not None:
            return ts
        if self.sides_interleave:
            # try the second head:
            ts = self.bam_get_block_on_track(prev_ts, interleave, original_track + self.sides_interleave)
            if ts is not None:
                return ts
        # now scan other tracks:
        candidate = original_track
        if candidate < self.directory_ts[0]:
            change = -1
        else:
            change = 1
        while True:
            # get new candidate track
            candidate += change
            # back to where we started? -> give up
            if candidate == original_track:
                raise DiskFullException()
            # at limit? -> go to directory and reverse search direction
            if candidate < self.mintrack or candidate > maxtrack:
                candidate = self.directory_ts[0]
                change = -change
                continue
            # scan candidate track
            ts = self.bam_get_block_on_track(prev_ts, interleave, candidate)
            if ts is not None:
                return ts
            if self.sides_interleave:
                # try the second head:
                ts = self.bam_get_block_on_track(prev_ts, interleave, candidate + self.sides_interleave)
                if ts is not None:
                    return ts

    def linkptrs_release_all(self, start_ts: (int, int)) -> None:
        """
        Free all connected blocks, following the chain of link pointers.
        """
        all_ts = set()
        for ts, block in self.linkptrs_follow(start_ts):
            all_ts.add(ts)
        self.bam_release_blocks(all_ts)

    def filetype(self, filetype: int) -> bytes:
        """
        Convert CBM file type to text representation.
        """
        ret = b" " if filetype & 128 else b"*"
        ft = filetype & 15
        if ft in self.filetypes:
            ret += self.filetypes[ft]   # supported file types are decoded
        else:
            ret += b"0X%X" % ft # unsupported file types are shown as hex digit
        ret += b"<" if filetype & 64 else b" "
        return ret

    def enter(self, which: int) -> object:
        """
        Enter partition or subdirectory.
        This works for 1581-style "CBM" partitions, for CMD-style partitions
        and for subdirectories in CMD native partitions.
        The partition/subdirectory is specified via its directory index.
        The return value is the disk image to use for all further accesses.
        """
        sys.exit("Error: This image type does not support entering partitions/directories.")

# FIXME - a lot of this code could be re-used for _creating_ side sectors, so
# better split this up into "make correct side sectors" and "compare correct
# side sectors to actual side sectors":
    def check_side_sectors(self, std_chain: list, second_chain: list, record_length: int) -> None:
        """
        Check all contents of a REL file's side sectors.
        "std_chain" is the t/s list of file's normal blocks, without actual block data.
        "second_chain" is the t/s list of side sectors, along with content.
        """
        # is the first side sector a "super side sector"?
        # (used on 8250, 1581 and CMD native partitions):
        if second_chain[0][1][2] == 0xfe:
            print("First side sector is a super side sector.")
            super_side_sector = second_chain[0][1]  # [1] -> contents
            second_chain = second_chain[1:]
        else:
            super_side_sector = None
            if len(second_chain) > 6:
                print(f"Found {len(second_chain)} side sectors, max is 6!")
                return
        # convert all t/s values to bytes for easier checking:
        std_table = bytes()
        for ts in std_chain:
            std_table += bytes(ts)
        ss_table = bytes()
        for ts, block in second_chain:
            ss_table += bytes(ts)
        # check contents of super side sector:
        if super_side_sector is not None:
            sss_table = bytes()
            ss_copy = ss_table
            while ss_copy:
                sss_table += ss_copy[0:2]   # pointer to first of six
                ss_copy = ss_copy[12:]  # ignore the other five
            sss_table += bytes(126*2)   # append zeroes
            sss_table = sss_table[:126*2]   # cut to correct length
            if super_side_sector[3:255] != sss_table:
                print("Super side sector has wrong pointers to side sector groups.")
        # check contents of "normal" side sectors:
        for index, pair in enumerate(second_chain):
            ts, block = pair
            if block[0] == 0:
                if block[1] != 16 + len(std_table) - 1:
                    print(f"Final side sector ({index}) has wrong end offset ({block[1]}).")
            if block[2] != index % 6:
                print(f"Side sector {index} has wrong number ({block[2]}).")
            if block[3] != record_length:
                print(f"Side sector {index} has wrong record length ({block[3]}).")
            # is this a new group? -> get a new table of pointers to side sectors
            if index % 6 == 0:
                ss6_table = (ss_table[:12] + bytes(12))[:12]
                ss_table = ss_table[12:]
            # check side sector's pointers to side sectors:
            if block[4:16] != ss6_table:
                print(f"Side sector {index} has wrong pointers to side sectors.")
            # check side sector's pointers to data blocks:
            pointer_bytes = std_table[:240]
            std_table = std_table[240:]
            if block[16:16+len(pointer_bytes)] != pointer_bytes:
                print(f"Side sector {index} has wrong pointers to data blocks.")

    def check_file(self, index: int) -> None:
        """
        check blocks of file (only really useful for REL files)
        """
        bin30 = self.direntry_get(index)
        # read fields
        cbm_type = bin30[0]
        t_s = bin30[1], bin30[2]
        #file_name = bin30[3:19]
        second_t_s = bin30[19], bin30[20]   # t/s of side sector (REL) or info block (GEOS)
        record_length = bin30[21]   # GEOS: 0=normal file, 1=VLIR file
        geos_type = bin30[22:23]
        block_count = int.from_bytes(bin30[28:30], "little")
        # std file body:
        print("Checking block chain:", end="")
        std_chain = list(self.linkptrs_follow(t_s, display=True, include_blocks=False))
        actual_blocks = len(std_chain)
        if second_t_s[0]:
            # side sectors:
            print("Checking side sector chain:", end="")
            second_chain = list(self.linkptrs_follow(second_t_s, display=True))
            actual_blocks += len(second_chain)
        # check block count
        if actual_blocks != block_count:
            print(f"Block count in directory is wrong ({block_count}), should be {actual_blocks}.")
        # for REL files, check contents of side sectors:
        if (cbm_type & 15) == 4:
            self.check_side_sectors(std_chain, second_chain, record_length)
        #

    def delete_file(self, index: int) -> None:
        """
        release all blocks of file and remove its directory entry
        """
        bin30 = self.direntry_get(index)
        # read fields
        cbm_type = bin30[0]
        t_s = bin30[1], bin30[2]
        #file_name = bin30[3:19]
        second_t_s = bin30[19], bin30[20]   # t/s of side sector (REL) or info block (GEOS)
        #record_length = bin30[21]   # GEOS: 0=normal file, 1=VLIR file
        #geos_type = bin30[22:23]
        #block_count = int.from_bytes(bin30[28:30], "little")
        if cbm_type == 0:
            print("Directory entry is unused.")
            return
        # std file body:
        print("Reading block chain:", end="")
        chain = list(self.linkptrs_follow(t_s, display=True, include_blocks=False))
        if second_t_s[0]:
            # side sectors:
            print("Reading side sector chain:", end="")
            chain += list(self.linkptrs_follow(second_t_s, display=True, include_blocks=False))
        print("Releasing blocks.")
        self.bam_release_blocks(chain)
        # remove dir entry
        for entry in self.directory_read_entries():
            if entry[0] == index:
                # entry[1] is actual bin30 entry
                ts, offset = entry[2]
                dirblock = self.block_read(ts)
                dirblock[offset] = 0    # mark entry as free
                self.block_write(ts, dirblock)
                return
        raise Exception("This cannot happen, dir entry was found before but not now?!")

    # TODO: maybe add optional args for
    #   - use alternative sector interleave
    #   - interleave sides on 1571 (and 8250?)
    #   - only make a super side sector if needed
    def allocate_file(self, filesize: int, record_length: int=None, interleave: int=None) -> (list, (int, int)):
        """
        allocate enough blocks so a file of the given size can be stored.
        if record_length is given, create a chain of side sectors pointing to
        the data blocks.
        return a tuple: first item is list of track/sector tuples so caller can
        do the actual block writes, second item is t/s of first side sector.
        """
        # first check how many blocks this file would need:
        data_block_count = _ceil_div(filesize, 254)
        total_block_count = data_block_count
        if record_length is not None:
            if filesize % record_length != 0:
                raise Exception(f"File size ({filesize}) is not a multiple of record length ({record_length}).")
            side_sector_count = _ceil_div(data_block_count, 120)
            total_block_count += side_sector_count
            if self.use_super_side_sector:
                total_block_count += 1
        # does it fit?
        free_blocks = self.bam_read_free_blocks()
        if total_block_count > free_blocks["shown"]:
            raise DiskFullException()   # "Not enough free blocks to add file."
        # allocate blocks
        data_blocks = []
        side_sectors = []
        super_side_sector = None
        ts = None
        for ii in range(data_block_count):
            ts = self.bam_get_additional_block(ts, interleave)
            data_blocks.append(ts)
            # did we just start a fresh batch of 120 data blocks?
            if ii % 120 == 0:
                if record_length is not None:
                    side_sectors.append(self.bam_get_additional_block(ts, interleave))
                    # do we want a super side sector?
                    if self.use_super_side_sector and super_side_sector is None:
                        super_side_sector = self.bam_get_additional_block(ts, interleave)
        # TODO:
        # write data to sss and ss.
        # return list of all data ts and ts of first ss.

    def add_file(self, name: bytes, body: bytes, record_length: int=None, interleave: int=None) -> None:
        """
        add file to image
        """
        tup = self.allocate_file(len(body), record_length, interleave)
        print("not finished")

    def set_interleave_sides(self, new_state: bool) -> None:
        """
        set new state for "interleave sides" optimization, which is only useful
        on 1571 and 8250.
        """
        if self.allow_1571_8250_trick and new_state:
            self.sides_interleave = self.maxtrack / 2   # interleave sides
        else:
            self.sides_interleave = 0   # default

    # TODO: add function for "force 1571 interleave (6) on 1541 image"

    # TODO: add function for "force 1541 interleave (10) on 1571 image"

################################################################################
# CBM DOS 1.0 disk format, as used in non-upgraded CBM 2040 and CBM 3040 units.
# tracks 18..24 have 20 sectors, which is a very tight fit.
# 690 blocks total, 670 free, 152 dir entries.
# no REL files.
# file extension is .d67 (at least in VICE)
# t18s0, offset 2: 0x01 (-> format 1)

class _dos1(d64):
    name = "2040/3040 (DOS 1.0)"
    blocks_total = 690  # 670 free
    maxtrack = 35
    track_length_changes = {1: 21, 18: 20, 25: 18, 31: 17}
    header_ts = (18, 0) # where to find diskname (dos 1 had no pseudo id?)
    header_offset = 144 # where to find diskname (dos 1 had no pseudo id?)
    bam_blocks = [(18, 0)]
    bam_start_size_maxperblock = (4, 4, 35)
    directory_ts = (18, 1)
    std_max_dir_entries = 152   # for writing directory
    #std_directory_interleave =
    #std_file_interleave =
    filetypes = { 0:b"DEL", 1:b"SEQ", 2:b"PRG", 3:b"USR"}   # decoded file types
    use_super_side_sector = False

    def bam_display(self):
        """
        Display block availability map in human-readable format.
        """
        li = self._bam_to_list((18, 0), 5, 4, 1, 36)
        self._bam_print_transposed(li)

################################################################################
# CBM DOS 2.1 was used in CBM 4040 units and in upgraded 2040/3040 units.
# tracks 18..24 have 19 sectors, which is more robust.
# 683 blocks total, 664 free, 144 dir entries.
# new: REL files and @SAVE
# t18s0, offset 2: 0x41 0x00 (-> format "A"), directory says "2A"

class _dos2p1(d64):
    filetypes = { 0:b"DEL", 1:b"SEQ", 2:b"PRG", 3:b"USR", 4:b"REL" }    # decoded file types
    use_super_side_sector = False
    # this class does not contain anything else because the disk format is the
    # same as the one that was later used by 1541 and friends, so see below.
    # the only reason for this class is so both "1541" and "8050" can inherit
    # the REL file definition.
    # FIXME: DOS capabilities differ too much to be represented by inheritance,
    # so just add "supports_REL/supports_CBM/supports_DIR" vars!

################################################################################
# CBM DOS 2.6 was a merger of 2.1 and 2.5 (see below): basically it was a
# downgrade of CBM DOS 2.5 to a single drive, but using the same disk format as
# CBM DOS 2.1.
# so the same disk format was used in CBM 4040 units, upgraded 2040/3040 units,
# 2031, 4031, 1540, 1541, 1551 and 1570 units.
# file extension is mostly .d64, sometimes .d41
# t18s0, offset 2: 0x41 0x00 (-> format "A"), directory says "2A"

class _1541(_dos2p1):
    name = "1541"
    blocks_total = 683  # 664 free
    maxtrack = 35
    track_length_changes = {1: 21, 18: 19, 25: 18, 31: 17}
    header_ts = (18, 0) # where to find diskname and five-byte "pseudo id"
    header_offset = 144 # where to find diskname and five-byte "pseudo id"
    bam_blocks = [(18, 0)]
    bam_start_size_maxperblock = (4, 4, 35)
    directory_ts = (18, 1)
    std_max_dir_entries = 144   # for writing directory
    std_directory_interleave = 3
    std_file_interleave = 10

    def bam_release_blocks(self, set_of_ts: list) -> None:
        bamblock = self.block_read((18, 0))
        dirty = False
        for track, sector in set_of_ts:
            self._bamblock_release_block((track, sector), bamblock, 1, 4, 4)
            dirty = True
        if dirty:
            self.block_write((18, 0), bamblock)

    def bam_allocate_block(self, track: int, sector: int, exact: bool) -> (int, int) or None:
        return self._bamblock_allocate_block((track, sector), (18, 0), 4, 4, track - 1, exact=exact)

    def bam_display(self):
        """
        Display block availability map in human-readable format.
        """
        li = self._bam_to_list((18, 0), 5, 4, 1, 36)
        self._bam_print_transposed(li)

################################################################################
# info taken from: http://mhv.bplaced.net/diskettenformate.html
# there was a "CBM 2030" disk drive unit that used disk format "2B".
# it had the same number of tracks, but fewer sectors per track than "2A".
# capacity was 130k instead of 170k.
# unknown, but likely: t18s0, offset 2: 0x42 0x00 (-> format "B"), directory says "2B"

class _2030(_1541):
    name = "2030"
    #track_length_changes = {1: 19, ..., ...: 15}
    #blocks_total = ...so maybe about 610? with about 590 free? but "130k" would be 520!
    #std_max_dir_entries =

################################################################################
# CBM DOS 2.5 was used in CBM 8050 units.
# new: disk changes were automatically detected
# REL files were still limited to 720 data blocks (180 kB), but these disks may
# have been used in 8250 units, therefore REL files on a dos 2.5 disk may use
# the dos 2.7 format, i.e. with a super side sector - beware!
# file extension is .d80
# header block links to first bam block, last bam block links to first dir block.
# t39s0, offset 2: 0x43 0x00 (-> format "C"), directory says "2C"

class _8050(_dos2p1):
    name = "8050"
    blocks_total = 2083 # 2052 free
    maxtrack = 77
    track_length_changes = {1: 29, 40: 27, 54: 25, 65: 23}
    header_ts = (39, 0) # where to find diskname and five-byte "pseudo id"
    header_offset = 6   # where to find diskname and five-byte "pseudo id"
    bam_blocks = [(38, 0), (38, 3)] # one source says sector 1 instead of 3!
    bam_start_size_maxperblock = (6, 5, 50)
    directory_ts = (39, 1)
    std_max_dir_entries = 224   # for writing directory
    #std_directory_interleave = 3   ?
    #std_file_interleave = 10       ?

    def bam_release_blocks(self, set_of_ts: list) -> None:
        bamblock380 = self.block_read((38, 0))
        bamblock383 = self.block_read((38, 3))
        dirty380 = False
        dirty383 = False
        for track, sector in set_of_ts:
            if track < 51:
                self._bamblock_release_block((track, sector), bamblock380, 1, 6, 5)
                dirty380 = True
            else:
                self._bamblock_release_block((track, sector), bamblock383, 51, 6, 5)
                dirty383 = True
        if dirty380:
            self.block_write((38, 0), bamblock380)
        if dirty383:
            self.block_write((38, 3), bamblock383)

    def bam_allocate_block(self, track: int, sector: int, exact: bool) -> (int, int) or None:
        if track < 51:
            ts = self._bamblock_allocate_block((track, sector), (38, 0), 6, 5, track - 1, exact=exact)
        else:
            ts = self._bamblock_allocate_block((track, sector), (38, 3), 6, 5, track - 51, exact=exact)
        return ts

    def bam_display(self):
        """
        Display block availability map in human-readable format.
        """
        li  = self._bam_to_list((38, 0), 7, 5,  1, 51)
        li += self._bam_to_list((38, 3), 7, 5, 51, 78)
        self._bam_print_list(li)

################################################################################
# CBM DOS 2.7 was used in CBM 8250 units and SFD-1001 units.
# new: support for double-sided discs
# new: super side sector, so REL files can use the whole disk. but as this is
# incompatible with CBM DOS 2.5, support for sss can be switched on/off via M-W
# command. so image files of this type may hold REL file with or without sss!
# file extension is .d82
# t39s0, offset 2: 0x43 0x00 (-> format "C"), directory says "2C"
# the "double-sidedness" may have been deduced from the bam blocks (number and content).

class _8250(_8050):
    name = "8250"
    blocks_total = 4166 # 4133 free
    maxtrack = 154
    track_length_changes = {1: 29, 40: 27, 54: 25, 65: 23, 77+1: 29, 77+40: 27, 77+54: 25, 77+65: 23}
    bam_blocks = [(38, 0), (38, 3), (38, 6), (38, 9)]
    allow_1571_8250_trick = True
    use_super_side_sector = True

    def bam_release_blocks(self, set_of_ts: list) -> None:
        bamblock380 = self.block_read((38, 0))
        bamblock383 = self.block_read((38, 3))
        bamblock386 = self.block_read((38, 6))
        bamblock389 = self.block_read((38, 9))
        dirty380 = False
        dirty383 = False
        dirty386 = False
        dirty389 = False
        for track, sector in set_of_ts:
            if track <= 50:
                self._bamblock_release_block((track, sector), bamblock380, 1, 6, 5)
                dirty380 = True
            elif track <= 100:
                self._bamblock_release_block((track, sector), bamblock383, 51, 6, 5)
                dirty383 = True
            elif track <= 150:
                self._bamblock_release_block((track, sector), bamblock386, 101, 6, 5)
                dirty386 = True
            else:
                self._bamblock_release_block((track, sector), bamblock389, 151, 6, 5)
                dirty389 = True
        if dirty380:
            self.block_write((38, 0), bamblock380)
        if dirty383:
            self.block_write((38, 3), bamblock383)
        if dirty386:
            self.block_write((38, 6), bamblock386)
        if dirty389:
            self.block_write((38, 9), bamblock389)

    def bam_allocate_block(self, track: int, sector: int, exact: bool) -> (int, int) or None:
        if track <= 50:
            ts = self._bamblock_allocate_block((track, sector), (38, 0), 6, 5, track - 1, exact=exact)
        elif track <= 100:
            ts = self._bamblock_allocate_block((track, sector), (38, 3), 6, 5, track - 51, exact=exact)
        elif track <= 150:
            ts = self._bamblock_allocate_block((track, sector), (38, 6), 6, 5, track - 101, exact=exact)
        else:
            ts = self._bamblock_allocate_block((track, sector), (38, 9), 6, 5, track - 151, exact=exact)
        return ts

    def bam_display(self):
        """
        Display block availability map in human-readable format.
        """
        li  = self._bam_to_list((38, 0), 7, 5,   1,  51)
        li += self._bam_to_list((38, 3), 7, 5,  51, 101)
        li += self._bam_to_list((38, 6), 7, 5, 101, 151)
        li += self._bam_to_list((38, 9), 7, 5, 151, 155)
        self._bam_print_list(li)

################################################################################
# disk format of 1541-with-40-tracks-support
# file extension is mostly .d64, sometimes .d41
#
# there are three or four different ways to store the additional 5*4 BAM bytes:
#   SpeedDOS uses offset 192, in the unused part of t18s0,
#   DolphinDOS uses offset 236? in the unused part of t18s0,
#   another system puts the data right after the old bam and moves the disk name!
# t18s0, offset 2: 0x41 0x00 (-> format "A"), directory says "2A"
# TODO: is there a way to recognize this format, short of trying to access tracks 36..40?

class _40track(_1541):
    name = "40-track 1541"
    blocks_total = 768  # 749 free
    maxtrack = 40

    def check_bam_counters(self):
        super().check_bam_counters(alt_maxtrack=35) # let 1541 class check the first side...
        # ...and now check extra tracks:
        print("FIXME: Checking BAM counters of tracks 36..40 is not implemented.")  # TODO
        # TODO: we could check all alternatives and if only one of those is
        # valid, consider that a successful format detection.
        # The problem is: "all zeroes" is a valid pattern in any case, it would
        # mean "all blocks are allocated".

    def bam_read_free_blocks(self):
        d = super().bam_read_free_blocks(alt_maxtrack=35)   # let 1541 class do the first side...
        # ...and now add extra tracks:
        print("FIXME: Reading BAM counters of tracks 36..40 is not implemented.")  # TODO
        #self._bam_fill_free_blocks_dict(d, (18, 0), FIXME, 4, 5, 36)
        # TODO - find a way to include "would show XYZ in a 1541 drive" info!
        return d

    #def bam_release_blocks(self, set_of_ts: list) -> None:
    # TODO: add this fn or a dummy because the 1541 method would overwrite the
    # disk name when track grows beyond 35!

    def bam_allocate_block(self, track: int, sector: int, exact: bool) -> (int, int) or None:
        if track <= 35:
            return super().bam_allocate_block(track, sector, exact) # just call 1541 method
        raise Exception("Allocation of tracks 36..40 is not yet supported!")    # TODO

################################################################################
# CBM DOS 3.0 in 1571 units and CBM DOS 3.1 in C128 DCRs extended the 1541 format
# to double-sided discs (-> 1328 blocks free), but:
# - on disc the format is still called "2A".
# - there are no super side sectors, so REL files are limited to 720 data blocks.
# - the BAM for the second side is split:
#       counters are at t18s0 (offset 221), bitmaps are at t53s0 (offset 0).
# - t53 (t18 on second side) is allocated completely, though only the first 105
#       bytes of t53s0 are used.
# (this "DOS 3" has nothing to do with the "DOS 3" used in D9060/D9090 hard disks)
# file extension is .d71 or .d64
# t18s0, offset 2: 0x41 0x80 (-> format "A" and flag bit), directory says "2A"
# the "0x80" above is the flag that signals "this disk is double-sided".

class _1571(_1541):
    name = "1571"
    blocks_total = 1366 # 1328 free
    maxtrack = 70
    track_length_changes = {1: 21, 18: 19, 25: 18, 31: 17, 35+1: 21, 35+18: 19, 35+25: 18, 35+31: 17}
    std_file_interleave = 6
    allow_1571_8250_trick = True

    def check_bam_counters(self):
        super().check_bam_counters(alt_maxtrack=35) # let 1541 class check the first side...
        # ...and now check second side:
        #_debug(1, "Checking 1571 BAM (second side)")
        # BAM for second side is split into two parts (totals at 18,0 and bitmaps at 53,0),
        # so we cannot use the std function:
        totals_block = self.block_read((18, 0))
        bits_block = self.block_read((53, 0))
        for entry in range(35):
            track = entry + 36
            total = totals_block[221 + entry]
            sum = 0
            for i in range(3):
                sum += _popcount(bits_block[entry * 3 + i])
            if total == sum:
                _debug(9, "track %d: %d free blocks" % (track, sum))
            else:
                print("BAM error for track %d: counter says %d, bitfield says %d!" % (track, total, sum), file=sys.stderr)
        # if t53 is free'd, VALIDATE on realhw throws "71, dir error" and does
        # not write fixed BAM to disk, so let's add a check for that:
        totals_t53 = totals_block[221 + 17] # offset 17 means "cylinder 18" i.e. track 53
        if totals_t53:
            print("BAM shows %d free blocks on t53, drive does not like this!" % totals_t53)

    def bam_read_free_blocks(self) -> dict:
        d = super().bam_read_free_blocks(alt_maxtrack=35)   # let 1541 class do the first side...
        # ...and now add second side:
        self._bam_fill_free_blocks_dict(d, (18, 0), 221, 1, 35, 36)
        # tested on real hw: if t53 is freed, its blocks are included in "BLOCKS FREE"
        # (but then VALIDATE will fail with "71, dir error, 53, 55", so don't do that).
        # TODO - find a way to include "would show XYZ in a 1541 drive" info!
        return d

    def release_1571_side2_block(self, ts: (int, int), totals_block: bytes, bitmaps_block: bytes) -> None:
        """
        Helper function to release a single block in BAM. This function is only
        used for blocks on second side of 1571 because their BAM is split over
        two blocks.

        ts: track and sector of block to release
        totals_block: block with totals (t18s0)
        bitmaps_block: block with bitmaps (t53s0)
        """
        # process args
        track, sector = ts
        # calculate offsets
        entry_index = track - 36    # this fn handles logical tracks 36..70
        totals_offset = 221 + entry_index   # table of totals is located at end of t18s0
        # bitmaps (three bytes each) are at start of t53s0:
        byte_offset, bit_value = self._bam_offset_and_bit(sector)
        bitmap_offset = entry_index * 3 + byte_offset
        # check/release:
        if bitmaps_block[bitmap_offset] & bit_value:
            raise Exception("Attempted to free a block (t%ds%d) that is already free." % (track, sector))
        else:
            bitmaps_block[bitmap_offset] |= bit_value
            totals_block[totals_offset] += 1
        # TODO - compare totals to maximum for this track?
        #raise Exception("BAM is corrupt, totals do not match bitmap.")

    def allocate_1571_side2_block(self, track: int, wanted_sector: int, exact: bool) -> (int, int) or None:
        """
        Helper function to allocate a single block in BAM. This function is only
        used for blocks on second side of 1571 because their BAM is split over
        two blocks.
        If block is available, allocate it and return t/s.
        If block is not available, return None.

        track, wanted_sector: track and sector of block to allocate
        totals_block: block with totals (t18s0)
        exact: if False, function may allocate and return a different sector from this track
        """
        # read bam blocks
        totals_block = self.block_read((18, 0))
        bitmaps_block = self.block_read((53, 0))
        # calculate offsets
        entry = track - 36  # zero-based offset, so subtract 36 instead of 35
        totals_offset = 221 + entry
        bitmaps_offset = 3 * entry
        cand_sector = wanted_sector # we start the search with the wanted sector...
        num_sectors = self.sectors_of_track(track)  # ...and this is where we wrap around
        while True:
            byte_offset, bit_value = self._bam_offset_and_bit(cand_sector)
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
                # write bam blocks
                self.block_write((18, 0), totals_block)
                self.block_write((53, 0), bitmaps_block)
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

    def bam_release_blocks(self, set_of_ts: list) -> None:
        bamblock180 = self.block_read((18, 0))
        bamblock530 = self.block_read((53, 0))
        dirty180 = False
        dirty530 = False
        for track, sector in set_of_ts:
            if track <= 35:
                self._bamblock_release_block((track, sector), bamblock180, 1, 4, 4)
                dirty180 = True
            else:
                self.release_1571_side2_block((track, sector), bamblock180, bamblock530)
                dirty180 = True
                dirty530 = True
        if dirty180:
            self.block_write((18, 0), bamblock180)
        if dirty530:
            self.block_write((53, 0), bamblock530)

    def bam_allocate_block(self, track: int, sector: int, exact: bool) -> (int, int) or None:
        if track <= 35:
            ts = self._bamblock_allocate_block((track, sector), (18, 0), 4, 4, track - 1, exact=exact)    # all at t18s0
        else:
            ts = self.allocate_1571_side2_block(track, sector, exact)   # totals at t18s0, bitmaps at t53s0
        return ts

    def bam_display(self):
        """
        Display block availability map in human-readable format.
        """
        li  = self._bam_to_list((18, 0), 5, 4,  1, 36)
        li += self._bam_to_list((53, 0), 0, 3, 36, 71)
        self._bam_print_transposed(li)

################################################################################
# CBM DOS 3.0 was used in D9060 and D9090 hard disk units
# (this "DOS 3" has nothing to do with the "DOS 3" used in 1570/1571 units)
# shit, these things actually have a track zero!
# first block (t0s0) seems to hold
#   pointers to bad sector list, dir, header, bam
#   and then ID:
# 00000000  00 01 00 ff 4c 0a 4c 14  01 00 39 30 00 00 00 00  |....L.L...90....|
# the directory says "3A", and there does not seem to be a format indicator at
# the start of the header block (...also not really needed with hard disks)

# D9090: 20 bam blocks, one every 8 tracks, starting at t1s0
# each contains 48 five-byte entries (counter byte and four bitmap bytes)
# 48*32 -> 12*128 -> each bam block manages twelve tracks * 128 blocks.
# 48*32 -> 8*192 -> each bam block manages eight tracks * 192 blocks.

# 64er 1988-11 says:
# "Nach Umrüstung [...] betragen die Speicherkapazitäten zwischen 2 mal 5 und 2 mal 16 MByte"
# picture shows:
# 2 logical drives (0 and 1), each with "38884 BLOCKS FREE".
# ...those could be: 8 heads, 256 blocks per cylinder, 153 tracks, 39168 blocks total

class _d90(_dos2p1):
    mintrack = 0
    maxtrack = 152
    use_super_side_sector = True    # this is just a guess...

    def _populate(self, imagefile):
        # first call parent class
        super()._populate(imagefile)
        # then get parameters from t0s0:
        t0s0 = self.block_read((0, 0))
        self.header_ts = (t0s0[6], t0s0[7]) # where to find diskname and five-byte "pseudo id"
        self.header_offset = 6  # where to find diskname and five-byte "pseudo id"
        self.directory_ts = (t0s0[4], t0s0[5])
        self.first_bam_ts = (t0s0[8], t0s0[9])

    def check_bam_counters(self) -> None:
        _debug(1, "Checking " + self.name + " BAM for internal consistency")
        bam_ts = self.first_bam_ts
        entry = 48  # so we know when to fetch the next bam block
        for cyl in range(self.mintrack, self.maxtrack+1):
            if entry == 48:
                bam_block = self.block_read(bam_ts)
                bam_ts = (bam_block[0], bam_block[1])
                entry = 0
            for head in range(self.heads):
                start = 16 + entry * 5
                total = bam_block[start]
                sum = 0
                for i in range(4):
                    sum += _popcount(bam_block[start + 1 + i])
                if total == sum:
                    _debug(9, "cyl %d, head %d: %d free blocks" % (cyl, head, sum))
                else:
                    print("BAM error for cyl %d, head %d: counter says %d, bitfield says %d!" % (cyl, head, total, sum), file=sys.stderr)
                entry += 1

    def bam_read_free_blocks(self) -> dict:
        _debug(3, "Reading free blocks counters")
        d = dict()
        d["all"] = 0
        d["reserved"] = 0
        d["shown"] = 0
        bam_ts = self.first_bam_ts
        entry = 48  # so we know when to fetch the next bam block
        for cyl in range(self.mintrack, self.maxtrack+1):
            if entry == 48:
                bam_block = self.block_read(bam_ts)
                bam_ts = (bam_block[0], bam_block[1])
                entry = 0
            freeblocks = 0
            for head in range(self.heads):
                start = 16 + entry * 5
                freeblocks += bam_block[start]
                entry += 1
            d[cyl] = freeblocks
            d["all"] += freeblocks
            if cyl:
                d["shown"] += freeblocks
            else:
                d["reserved"] += freeblocks
        return d

    def bam_release_blocks(self, set_of_ts: list) -> None:
        raise Exception("releasing blocks not implemented for D9060/D9090!")

    def bam_allocate_block(self, track: int, sector: int, exact: bool) -> (int, int) or None:
        raise Exception("allocating blocks not implemented for D9060/D9090!")

class _d9060(_d90):
    name = "D9060"
    blocks_total = 19584    # 19441 free -> track 0 is not counted though it is mostly free in BAM!
    # ...so 143 allocated, that's 128 (t0) + 2 (header + dir) + 13 (bam)
    heads = 4
    track_length_changes = {0: 128} # all tracks have 32 sectors (times 4 heads -> 128)

class _d9090(_d90):
    name = "D9090"
    blocks_total = 29376    # 29162 free -> track 0 is not counted though it is mostly free in BAM!
    # ...so 214 allocated, that's 192 (t0) + 2 (header + dir) + 20 (bam)
    heads = 6
    track_length_changes = {0: 192} # all tracks have 32 sectors (times 6 heads -> 192)

################################################################################
# CBM DOS 7.0 was used in 8061 and 8062 units (dual 8-inch floppy disk drives)
# docs say:
# tracks start at 0(!), sectors start at 1(!)
# t0s1 holds "volume label" with pointers to bad sector list, bam, directory.
# first byte holds sector number(!), "actually head *and* sector number",
# 0xff means "end of list",
# second byte holds track number or points to last byte used in block.
# ...so it is a lot like the hard disk format, except for the track/sector confusion.
# single-sided format uses 3 bam blocks, double-sided uses 6.

################################################################################
# CBM DOS 10.0 was used in 1581 units:
# - it knows about super side sectors, so REL files can use the whole disk.
# new filetype: "CBM" for a sequence of consecutive blocks
# file extension is .d81 or .d64
# TODO: display the "1581 autoboot" flag from the header block as debug info with the directory?
#   (that is bit 6 of byte 7 in t40s1 and t40s2)
# t40s0, offset 2: 0x44 0x00 (-> format "D"), directory says "3D"

class _1581(d64):
    name = "1581"
    blocks_total = 3200 # 3160 free
    maxtrack = 80
    track_length_changes = {1: 40}  # all tracks have 40 sectors
    header_ts = (40, 0) # where to find diskname and five-byte "pseudo id"
    header_offset = 4   # where to find diskname and five-byte "pseudo id"
    bam_blocks = [(40, 1), (40, 2)]
    bam_start_size_maxperblock = (16, 6, 40)
    directory_ts = (40, 3)
    std_max_dir_entries = 296   # for writing directory
    std_directory_interleave = 1    # 1581 uses interleave 1 because of track cache
    std_file_interleave = 1 # 1581 uses interleave 1 because of track cache
    filetypes = { 0:b"DEL", 1:b"SEQ", 2:b"PRG", 3:b"USR", 4:b"REL", 5:b"CBM" }  # decoded file types
    use_super_side_sector = True

    def bam_release_blocks(self, set_of_ts: list) -> None:
        bamblock1 = self.block_read((40, 1))
        bamblock2 = self.block_read((40, 2))
        dirty1 = False
        dirty2 = False
        for track, sector in set_of_ts:
            if track <= 40:
                self._bamblock_release_block((track, sector), bamblock1, 1, 16, 6)
                dirty1 = True
            else:
                self._bamblock_release_block((track, sector), bamblock2, 41, 16, 6)
                dirty2 = True
        if dirty1:
            self.block_write((40, 1), bamblock1)
        if dirty2:
            self.block_write((40, 2), bamblock2)

    def bam_allocate_block(self, track: int, sector: int, exact: bool) -> (int, int) or None:
        if track <= 40:
            ts = self._bamblock_allocate_block((track, sector), (40, 1), 16, 6, track - 1, exact=exact)
        else:
            ts = self._bamblock_allocate_block((track, sector), (40, 2), 16, 6, track - 41, exact=exact)
        return ts

    def enter(self, which: int) -> object:
        # enter 1581-style "CBM" partition
        bin30 = self.direntry_get(which)
        filetype = bin30[0] & 7
        if filetype != 5:
            sys.exit('Error: Chosen directory entry has type %d instead of 5 ("CBM")!' % filetype)
        start_track = bin30[1]
        start_sector = bin30[2]
        if start_sector != 0:
            sys.exit('Error: Chosen CBM partition starts at sector %d, not 0!' % start_sector)
        block_count = int.from_bytes(bin30[28:30], "little")
        if block_count < 120:
            sys.exit('Error: Chosen CBM partition has a size of %d blocks (need at least 120 blocks)!' % block_count)
        if block_count % 40:
            sys.exit('Error: Chosen CBM partition has a size of %d blocks (not a multiple of 40)!' % block_count)
        # ok, let's do it then
        self.header_ts = (start_track, start_sector)
        self.header_offset = 4
        self.bam_blocks = [(start_track, 1), (start_track, 2)]
        self.directory_ts = (start_track, 3)
        return self # we keep using this "1581" object

    def bam_display(self):
        """
        Display block availability map in human-readable format.
        """
        li  = self._bam_to_list((40, 1), 17, 6,  1, 41)
        li += self._bam_to_list((40, 2), 17, 6, 41, 81)
        self._bam_print_transposed(li)

################################################################################
# disk format of CMD native partitions
# new filetype: "DIR" for sub-directories (but no support for 1581-style "CBM" partitions)
# file extension is .dnp
# TODO: display the "native autoboot" flag (byte 7 of t1s2) as debug info with the directory?
# TODO: read and check "track number of last available track in partition" (byte 8 of t1s2)
# t1s1, offset 2: 0x48 0x00 (-> format "H"), directory says "1H"

# the first 64 blocks are allocated/reserved for bootblock/bam/rootdir, so:
# a 1-track DNP would use 256 blocks with 192 free
# a 12-track DNP (DD FD) uses 3072 blocks with 3008 free
# a 25-track DNP (HD FD) uses 6400 blocks with 6336 free
# a 50-track DNP (ED FD) uses 12800 blocks with 12736 free
# a 255-track DNP (maximum) uses 65280 blocks with 65216 free

class _cmdnative(d64):
    name = "CMD native"
    blocks_total = -1   # dynamic, multiple of 256
    # first 1+1+32+1 blocks are always allocated (bootblock, header, bam blocks, first dir block)
    # the first 64 blocks are not included in "blocks free" number (reserved for root dir)
#    maxtrack = -1   # dynamic!
    track_length_changes = {1: 256}  # all tracks have 256 sectors
    header_ts = (1, 1)  # where to find diskname and five-byte "pseudo id"
    header_offset = 4   # where to find diskname and five-byte "pseudo id"
#    bam_blocks = [(1, 2), (1, 3)]   # ...up to and including (1, 33)!
    bam_start_size_maxperblock = (0, 32, 8)
    directory_ts = (1, 34)
#    std_max_dir_entries = -1    # for writing directory (dynamic!)
    std_directory_interleave = 1
    std_file_interleave = 1
    filetypes = { 0:b"DEL", 1:b"SEQ", 2:b"PRG", 3:b"USR", 4:b"REL", 6:b"DIR" }  # decoded file types
    use_super_side_sector = True

    def __init__(self, filesize):
        self.maxtrack = filesize // 65536   # each track is 64 KiB
        if filesize % 65536:
            print("CAUTION, file size is not a multiple of 64 KiB, so there is unused space at the end!")
        self.blocks_total = filesize // 256
        self.bam_counters = dict()

    def check_bam_counters(self) -> None:
        # CMD native format does not have "free blocks on track" counters,
        # so they cannot disagree with bitmaps, so check always succeeds:
        pass

    def _fake_bam_counters(self) -> None:
        # BAM bitmaps are in blocks (1,2) up to and including (1,33),
        # but there are no counters (because 0..256 does not fit in a byte),
        # so we count the bits ourself:
        bits_block = None
        bytes_done = 0  # counter so we know when we reached end of "reserved" area
        for track in range(1, self.maxtrack + 1):
            entry = track & 7
            if (entry == 0) or (bits_block == None):
                bits_block = self.block_read((1, (track >> 3) + 2))
            sum = 0
            for i in range(32):
                sum += _popcount(bits_block[entry * 32 + i])
                bytes_done += 1
                if bytes_done == 8:
                    reserved_area_sum = sum # first 64 blocks
            _debug(9, "track %d: %d free blocks" % (track, sum))
            self.bam_counters[track] = sum
        # a correct directory display does not include the first 64 blocks in
        # its "free" result, so we store that number as an extra key and later
        # it can be subtracted when needed:
        self.bam_counters["reserved_area"] = reserved_area_sum

    def bam_read_free_blocks(self) -> dict:
        """
        Return dictionary of free blocks per track.
        There are three additional keys, "shown", "reserved" and "all", where
        "shown" + "reserved" == "all"
        """
        _debug(3, "Reading free blocks dict")
        if not self.bam_counters:
            self._fake_bam_counters()
        d = dict()
        sum = 0
        for track in range(1, self.maxtrack + 1):
            count = self.bam_counters[track]
            d[track] = count
            sum += count
        reserved_area_sum = self.bam_counters["reserved_area"]
        d["all"] = sum
        d["reserved"] = reserved_area_sum
        d["shown"] = sum - reserved_area_sum
        return d

    def _bam_offset_and_bit(self, sector: int) -> (int, int):
        """
        Convert sector number to byte offset and bit value for accessing BAM bitmap.
        """
        # bitmap bytes are little-endian (just like in all other formats), but
        # most significant bit goes first (which differs from all other formats):
        return sector >> 3, 128 >> (sector & 7)

    #def bam_release_blocks(self, set_of_ts: list) -> None:
        # TODO

    #def bam_allocate_block(self, track: int, sector: int, exact: bool) -> (int, int) or None:
        # TODO

    def enter(self, which: int) -> object:
        # enter CMD-style sub-directory
        bin30 = self.direntry_get(which)
        filetype = bin30[0] & 7
        if filetype != 6:
            sys.exit('Error: Chosen directory entry has type %d instead of 5 ("DIR")!' % filetype)
        start_track = bin30[1]
        start_sector = bin30[2]
        # ok, let's do it then
        self.header_ts = (start_track, start_sector)
        self.header_offset = 4
        headerblock = self.block_read(self.header_ts)
        self.directory_ts = (headerblock[0], headerblock[1])
        return self # we keep using this "CMD native" object

    def bam_display(self) -> None:
        """
        Display block availability map in human-readable format.
        """
        print("CAUTION, this BAM is shown in hex format and values are big-endian!")
        print("(Example: 0x1fffffff... means sectors 0, 1 and 2 are allocated)")
        offset = 32 # skip "track 0" entry of bam because it holds other data
        track = 1
        for bamsector in range(2, 34): # t1s2..33 are bam blocks
            bamblock = self.block_read((1, bamsector))
            while offset < 256:
                s = "t%03d: 0x" % track
                for byte in range(32):
                    s += "%02x" % bamblock[offset + byte]
                print(s)
                if track >= self.maxtrack:
                    return
                offset += 32
                track += 1
            offset = 0
        raise Exception("BUG - too many tracks in DNP")

################################################################################
# disk formats of CMD HD/FD images (CMD called them "partitionable formats"):
# this format started with the CMD HD and was also used in CMD FDs and RAMLink.
# up to 254 partitions possible on HD, 31 on RAMLink and FDs.
# partitions are stored sequentially without any gaps, beginning with the
# system partition. the system partition holds the DOS in case of HD and is
# empty in case of FD.
# 0 specifies the currently selected partition,
# 1..254 specify partitions 1..254 (one of which is the "currently selected
#   partition" and one of which is selected per default after power-on)
# 255 specifies the system partition (located first, so stored at offset 0)
# FD formats:
#   file extension .d1m, size is 829440
#   file extension .d2m, size is 1658880
#   file extension .d4m, size is 3317760
# all use 81 tracks, where track 81 holds the partition table
# supported partition types are:
#   0: not created
#   1: native   (takes N * 256 blocks)
#   2: 1541emu  (takes 683+1 blocks, see below why)
#   3: 1571emu  (takes 1366 blocks on HD, 1366+2 on FD, see below why)
#   4: 1581emu  (takes 3200 blocks)
#       -> each partition must start on a hardware sector.
#       HD uses 512-byte sectors, FD uses 1024-byte sectors, therefore sizes of
#       1541- and 1571-emulation partitions are rounded up.
# only supported by HD:
#   5: 1581 CP/M emulation
#   6: print buffer
#   7: foreign mode
# TODO: do not display any "blocks free" number (or maybe do, but use sensible data)
# TODO: mark the default partition in the directory

class _cmdpartitionable(d64):
    """
    virtual base class for the three CMD FD series formats.
    """
    maxtrack = 81
    directory_ts = (1, 0)

    def _populate(self, imagefile):
        # first call parent class
        super()._populate(imagefile)
        # then set partition values so we can access partition directory
        self.imagefile.partition_set(self.partition)

    # partition table does not really have header and id fields, so fake them:
    def bam_read_header_fields(self) -> (int, bytes, bytes):
#       return 255, b"CMD HD          ", b"HD 1H"   # HD
        return 255, b"CMD FD          ", b"FD 1H"   # FD
#       return 255, b"CMD RAMLINK     ", b"RL 1H"   # RAMLink   FIXME: check values!

    def bam_read_free_blocks(self) -> dict:
        # just use all zero for now:
        d = dict()
        for track in range(1, self.maxtrack + 1):
            d[track] = 0
        d["all"] = 0
        d["reserved"] = 0
        d["shown"] = 0
        return d

    def check_bam_counters(self):
        # partitions use consecutive space, so there is not really a BAM, so
        # check always succeeds:
        pass

    partition_types = { 0:b"*DEL", 1:b" NAT", 2:b"  41", 3:b"  71", 4:b"  81", 5:b" 81C", 255:b" SYS" }
    def filetype(self, filetype: int) -> bytes:
        """
        Convert CBM file type to text representation.
        """
        if filetype in self.partition_types:
            t = self.partition_types[filetype] + b" "
        else:
            t = b" 0X%2X" % filetype    # unsupported types are shown as hex
        return t

    def direntry_cook(self, index: int, bin30: bytes or bytearray) -> (bool, int, str, str, str):
        """
        Convert entry from partition directory to what is shown to user
        """
        in_use = bin30[0] != 0
        partition_type = self.filetype(bin30[0])
        unused1 = bin30[1:3]    # normally t/s
        partition_name = bin30[3:19]
        start_index = int.from_bytes(bin30[19:22], "big") * 512
        unused2 = bin30[22:27]
        size_bytes = int.from_bytes(bin30[27:30], "big") * 512
        optional = " : " + unused1.hex(" ")
        optional += " : %08x : " % start_index
        optional += unused2.hex(" ")
        optional += " : %08x " % size_bytes
        return in_use, index, partition_name, partition_type, optional

    def enter(self, which: int) -> object:
        # enter one of the partitions
        supported_types = {
            1: _cmdnative,
            2: _1541,
            3: _1571,
            4: _1581
        }
        bin30 = self.direntry_get(which)
        partition_type = bin30[0]
        if partition_type not in supported_types:
            sys.exit('Error: Chosen directory entry has type %d instead of 1/2/3/4!' % partition_type)
        chosen_type = supported_types[partition_type]
        # ok, let's do it then
        start_block = int.from_bytes(bin30[19:22], "big") * 2   # convert 512-byte sectors to 256-byte blocks
        block_count = int.from_bytes(bin30[27:30], "big") * 2   # convert 512-byte sectors to 256-byte blocks
        # create new handler object:
        if chosen_type == _cmdnative:
            new_obj = _cmdnative(block_count * 256)
        else:
            new_obj = chosen_type()
        # set partition values so new handler works correctly:
        self.imagefile.partition_set((start_block, block_count))
        # pass "partition" to handler:
        new_obj._populate(self.imagefile)
        return new_obj  # tell caller to use the new handler

    def geos_get_border_block_ts(self) -> (int, int) or None:
        return None

class _cmdfd1m(_cmdpartitionable):
    name = "CMD DD partitioned"
    blocks_total = 3240 # 3200 for partitions plus 40 on track 81 (holds partition table)
    # a DNP on a DD disk only uses 3072 blocks, because 3200 is not a multiple of 256.
    # so the last 128 blocks are unused in this case.
    track_length_changes = {1: 40}  # each track has 40 blocks (2 heads * 5 sectors * 1024 bytes)
    partition = (80 * 40 + 8, 4)

class _cmdfd2m(_cmdpartitionable):
    name = "CMD FD-2000 partitioned"
    blocks_total = 6480 # 6400 for partitions plus 80 on track 81 (holds partition table)
    track_length_changes = {1: 80}  # each track has 80 blocks (2 heads * 10 sectors * 1024 bytes)
    partition = (80 * 80 + 8, 4)

class _cmdfd4m(_cmdpartitionable):
    name = "CMD FD-4000 partitioned"
    blocks_total = 12960    # 12800 for partitions plus 160 on extra track (hold partition table)
    track_length_changes = {1: 160} # each track has 160 blocks (2 heads * 20 sectors * 1024 bytes)
    partition = (80 * 160 + 8, 4)

################################################################################
# wrapper stuff

# collect supported sizes so files can be identified.
# cmd native has a block count of -1 because its size is dynamic. so it cannot
# be identified by file size, but we still need it in this list so it is shown
# by list_formats().
_supported = (_dos1, _1541, _40track, _1571, _1581, _8050, _8250, _d9060, _d9090, _cmdnative, _cmdfd1m, _cmdfd2m, _cmdfd4m)
_type_of_size = dict()
for imgtype in _supported:
    _type_of_size[imgtype.blocks_total * 256] = (imgtype, False)    # no error info
    _type_of_size[imgtype.blocks_total * 257] = (imgtype, True) # with error info
_largest_size = 255*256*256 # for DNP, originally was max(_type_of_size.keys())
del imgtype # we do not want this to pop up in the online help...

def list_formats():
    for type in _supported:
        print(type.name)

def DiskImage(filename: str, img_mode = ImgMode.READONLY) -> d64:
    """
    Identify disk image file and return as d64 object.
    If img_mode is 'ImgMode.WRITEBACK', you need to call writeback() to flush
    changes to file.
    """
    # open file and read data
    mode = "r+b" if img_mode in (ImgMode.WRITEBACK, ImgMode.WRITETHROUGH) else "rb"
    fh = open(filename, mode)
    body = fh.read(_largest_size + 1)   # 1 more than largest supported type
    # FIXME - old version intercepted OSError! fix callers so they do that!
    body = bytearray(body)  # allow changes
    # determine type of disk image file
    filesize = len(body)
    if filesize in _type_of_size:
        img_type, error_info = _type_of_size[filesize]
        obj = img_type()
        hint = ""
    else:
        # FIXME: a d81 with 81 tracks (as used by some geos install disk) is the
        # same size as a d1m, so add code to tell them apart!
        # valid sizes of CMD native partitions are 1..255 * 64 KiB:
        if (0 < filesize <= 0xff0000) and (filesize & 0xffff == 0):
            obj = _cmdnative(filesize)
            error_info = False
            hint = f" ({filesize >> 16} tracks)"
        else:
            raise Exception("Could not process " + filename + ": Image type not recognised")
    _debug(1, filename, "is a", obj.name, "disk image" + hint + (" with error info." if error_info else "."))
    # FIXME: all the stuff above could be moved to imagefile.__init__(), right?
    img = imagefile(fh, body, img_mode, has_error_chunk=error_info) # create image file object...
    obj._populate(img)  # ...and pass to dos/drive object
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
_ANSI_REVERSE = "\033[7m"
_ANSI_RVS_OFF = "\033[27m"

def from_petscii(bindata: bytes, second_charset: bool) -> str:
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
            ret += _ANSI_REVERSE
        ret += charset[b]
        if revs:
            ret += _ANSI_RVS_OFF
    return ret

def _quote(name16: bytes) -> bytes:
    """
    Helper function to add opening and closing quotes at correct positions.
    """
    name16 += b"\xa0"   # appending a shift-space and then replacing the first...
    return b'"' + name16.replace(b"\xa0", b'"', 1)   # ...should do the trick.

def show_directory(img, show_all=False, second_charset=False, long=False) -> None:
    """
    Show directory of image file.
    """
    drive, name, id5 = img.bam_read_header_fields()
    name = from_petscii(name, second_charset)
    id5 = from_petscii(id5, second_charset)
    print('    ', drive, _ANSI_REVERSE + '"' + name + '"', id5 + _ANSI_RVS_OFF)
    nonempty = 0
    empty = 0
    blocks_total = 0    # for debug output
    for entry in img.directory_read_entries(include_invisible=True):
        index = entry[0]
        bin30 = entry[1]
        in_use, line_number, name, type, optional = img.direntry_cook(index, bin30)
        if in_use:
            nonempty += 1
        else:
            empty += 1
            if show_all == False:
                continue
        line = ('%3d: ' % index) + str(line_number).ljust(4) + " "
        line += from_petscii(_quote(name), second_charset)
        line += from_petscii(type, second_charset)
        if long:
            line += optional
        print(line)
        if in_use:
            blocks_total += line_number
    freeblocks = img.bam_read_free_blocks()
    #print(freeblocks)
    shown_free = freeblocks["shown"]
    print("     %d blocks free (+%d reserved)" % (shown_free, freeblocks["reserved"]))
    print("(%d directory entries, +%d empty)" % (nonempty, empty))
    #_debug(1, "%d + %d = %d" % (blocks_total, shown_free, blocks_total + shown_free))

def extract_block_sequence(img: d64, ts: (int, int), outname: str, blockcount: int) -> None:
    """
    Extract a number of consecutive blocks to file (CBM "partitions")
    """
    print("TODO: extract %d blocks beginning at t%ds%d to <%s>." % (blockcount, ts[0], ts[1], outname))

def extract_block_chain(img: d64, ts: (int, int), outname: str) -> None:
    """
    Extract block chain to file
    """
    print("extracting", outname)
    body = bytes()
    for ts, datablock in img.linkptrs_follow(ts):
        if datablock[0]:
            body += datablock[2:]
        else:
            body += datablock[2:datablock[1]+1]
    try:
        file = open(outname, "wb")
        file.write(body)
    except OSError as e:
        print("Error: Could not open " + outname + ": " + e.strerror, file=sys.stderr)
        return
    except Exception as e:
        print(e, "\n", file=sys.stderr)
        return
    #print("length:", len(body))

def extract_all(img: d64, full=False) -> None:
    """
    Extract all files to current directory.
    """
    for entry in img.directory_read_entries(include_invisible=full):
        index = entry[0]
        bin30 = entry[1]
        filetype = bin30[0]
        ts = bin30[1], bin30[2]
        cbmname = from_petscii(bin30[3:19], second_charset=True)
        cbmname = cbmname.replace("=", "=3D")
        cbmname = cbmname.replace("/", "=2F")
        outname = ('%03d-' % index) + cbmname
        outname += "."
        outname += from_petscii(img.filetype(filetype), second_charset=True)
        alt_ts = bin30[19], bin30[20]   # side sector chain or GEOS info block
        # FIXME: make this into some "extract entry" method so there is no need
        # to check for REL/CBM (which would fail in CMD partition tables anyway)
        # REL file?
        if (filetype & 15) == 4:
            # for REL files, add record size to name
            outname += "%03d" % bin30[21]
            extract_block_chain(img, ts, outname)
            # and for now, also extract the side sector chain:
            extract_block_chain(img, alt_ts, outname + ".ssc")
            continue
        # "CBM"-style partition?
        if (filetype & 15) == 5:
            # CBM files are partitions, i.e. without link pointers
            blockcount = int.from_bytes(bin30[28:30], "little")
            extract_block_sequence(img, ts, outname, blockcount)
            continue
        # extract GEOS info block
        if alt_ts[0]:
            extract_block_chain(img, alt_ts, outname + ".infoblock")
        # GEOS VLIR file?
        if bin30[21] == 1 and bin30[22] != 0:
            vlir = list(img.linkptrs_follow(ts))
            if len(vlir) == 1:
                # extract all 127 possible records with name postfix:
                dummy, vlir_block = vlir[0]
                for fork in range(1, 128):
                    ts = vlir_block[2 * fork], vlir_block[2 * fork + 1]
                    if ts[0]:
                        extract_block_chain(img, ts, outname + '.%03d' % fork)
                    else:
                        pass    # FIXME - 00 ff means ignore, 00 00 means end?
                continue
            print("looks like a GEOS VLIR file, but isn't! (too long)")
        # everything else:
        extract_block_chain(img, ts, outname)

def _process_file(file: str, show_all=False, second_charset=False, long=False, extract=False, path=None) -> None:
    if path == None:
        path = []
    try:
        img = DiskImage(file)
    except OSError as e:
        print("Error: Could not open " + file + ": " + e.strerror, file=sys.stderr)
        return
    except Exception as e:
        print(e, "\n", file=sys.stderr)
        return
    while path:
        component = path.pop(0)
        img = img.enter(component)
    img.check_bam_counters()    # TODO: add cli arg(s) to enable/disable this
    if extract:
        extract_all(img, full=show_all)
    else:
        show_directory(img, show_all=show_all, second_charset=second_charset, long=long)

def process_path(path: str) -> list:
    if not path:
        return []
    ret = []
    try:
        parts = path.split("/")
        for p in parts:
            num = int(p)
            if num < 0:
                sys.exit("Error: Path must be given as slash-separated, non-negative integers.")
            ret.append(num)
    except Exception:
        sys.exit("Error: Path must be given as slash-separated integers.")
    return ret

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
    parser.add_argument("-e", "--enter", metavar="PATH", type=str, help="enter partition/directory/...")
#    parser.add_argument("-f", "--full", action="store_true", help="same as -a -l")
    parser.add_argument("-a", "--all", action="store_true", help="show deleted dir entries as well")
    parser.add_argument("-l", "--long", action="store_true", help="show long directory entries")
    parser.add_argument(      "--list", action="store_true", help="list supported formats")
    parser.add_argument("-x", "--extract", action="store_true", help="extract all files to current directory")
    parser.add_argument("files", metavar="IMAGEFILE.D64", nargs='+', help="Disk image file.")
    args = parser.parse_args()
    _debuglevel += args.debug
    if args.list:
        list_formats()
    path = process_path(args.enter)
    for file in args.files:
        _process_file(file, show_all=args.all, second_charset=not args.charset, long=args.long, extract=args.extract, path=path[:])

if __name__ == "__main__":
    _main()
