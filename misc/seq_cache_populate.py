#!/usr/bin/env python3

# The MIT License

# Copyright (c) 2025 Genome Research Ltd.
# Author: Ruben Vorderman <r.h.p.vorderman@lumc.nl>

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# This file is a faster reimplementation of seq_cache_populate.pl
# <https://github.com/samtools/samtools/tree/develop/misc/seq_cache_populate.pl>
"""
Import references into a cram reference cache from fasta files.

When run with a list of fasta files, this program reads the files and stores
the sequences within it in the reference cache under directory <dir>.  The
sequences in the cache are stored with names based on the MD5 checksum
of the sequence.

By default, sequences are stored in a hierarchy two directories deep, to
keep the number of items in a single directory to a reasonable number.  This
depth can be changed using the -subdirs option.

If the -find option is used, the program will scan the given directory tree.
Any files that appear to be fasta (by looking for a first line starting with
'>' followed by something that looks like DNA sequence) will be read and
added to the reference cache.  The traversal will ignore symbolic links.

Samtools/htslib can be made to use the cache by appropriate setting of the
REF_PATH environment variable.  For example, if seq_cache_populate was run
using options '-root /tmp/ref_cache -subdirs 2', setting REF_PATH to
'/tmp/ref_cache/%2s/%2s/%s' should allow samtools to find the references that
it stored.

Note that if no REF_PATH is specified, htslib will default to downloading from
the EBI reference server and caching locally (see the samtools(1) man page for
details), defaulting to $HOME/.cache/hts-ref/%2s/%2s/%s.  This is functionally
equivalent to running this tool with '-root $HOME/.cache/hts-ref -subdirs 2'.
"""
import argparse
import hashlib
import logging
import mmap
import os
import shutil
import sys
import tempfile
from typing import Iterable, Iterator, List, Tuple

DEFAULT_NUMBER_OF_SUBDIRS = 2
MAX_MEMORY_SIZE = 256 * 1024 * 1024


class SpooledFile:
    """A spooled file that can reuse an existing memory map."""
    def __init__(self, spool: mmap.mmap, dir: str):
        self._file = spool
        self._file.seek(0)
        self._filename = None
        self._fd = None
        self._mmap = spool
        self._dir = dir

    def write(self, b, /) -> int:
        try:
            return self._file.write(b)
        except ValueError as e:
            # If the mmap is full a ValueError data out of range will be raised.
            if self._filename:
                raise e
        self._fd, self._filename = tempfile.mkstemp(dir=self._dir)
        self._file = open(self._fd, "wb+")
        self._file.write(memoryview(self._mmap)[:self._mmap.tell()])
        return self._file.write(b)

    def save_to(self, filename):
        if self._filename:
            self._file.close()
            shutil.move(self._filename, filename)
            self.close()
        else:
            with open(filename, "wb") as f:
                f.write(memoryview(self._mmap)[:self._mmap.tell()])

    def getvalue(self) -> bytes:
        current_pos = self._file.tell()
        self._file.seek(0)
        return self._file.read(current_pos)

    def closed(self):
        return self._file is None

    def close(self):
        if self._fd:
            self._file.close()
        self._file = None
        self._filename = None
        self._mmap = None
        self._fd = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


def fasta_to_blocks(fasta_file, block_size: int = 128 * 1024):
    """
    Yields Fasta in binary blocks. To separate contigs, the names are output
    as a block starting with ">".
    """
    if not file_is_fasta(fasta_file):
        raise ValueError("File is not a valid FASTA.")
    with open(fasta_file, "rb") as fasta:
        while True:
            block = fasta.read(block_size)
            if not block:
                return
            while block:
                # Run a loop to detect if multiple contigs are in a block.
                name_index = block.find(b">")
                if name_index == -1:
                    # No contigs left
                    yield block
                    break
                if name_index > 0:
                    yield block[:name_index]
                name_end = block.find(b"\n", name_index)
                if name_end == -1:
                    # Name is only partially present in block. Use readline
                    # as a shorthand to move the file exactly to the start
                    # of the sequence.
                    name = block[name_index:] + fasta.readline()
                    yield name
                    break
                name = block[name_index:name_end]
                yield name
                block = block[name_end:]


def name_block_to_name(name_block: bytes):
    full_name = name_block.decode("utf-8")
    name, *comments = full_name.split()
    return name.lstrip(">")


def fasta_to_contig_files_and_md5sums(
        fasta_file,
        memory_spool: mmap.mmap,
        contig_dir=tempfile.gettempdir()
) -> Iterator[Tuple[SpooledFile, str, str]]:
    """Create a spooled in memory file for each contig and calculate each md5sum"""
    if contig_dir and not os.path.exists(contig_dir):
        os.makedirs(contig_dir, mode=0o775)
    block_iter = fasta_to_blocks(fasta_file)
    try:
        name = name_block_to_name(next(block_iter))
    except StopIteration:
        return
    eof = False
    while not eof:
        block = b""
        hasher = hashlib.md5()
        tmp = SpooledFile(memory_spool, contig_dir)
        try:
            for block in block_iter:
                if block.startswith(b">"):
                    # Discard name block and continue with sequence.
                    break
                # Remove newlines as the created contig files should be a
                # raw sequence.
                block = block.replace(b"\n", b"").upper()
                tmp.write(block)
                hasher.update(block)
            else:  # No break, so no new contig was found
                eof = True
            yield tmp, hasher.hexdigest(), name
        finally:
            tmp.close()
        name = name_block_to_name(block)


def md5sum_to_dir_and_file(
        md5sum: str,
        number_of_subdirs: int = DEFAULT_NUMBER_OF_SUBDIRS
) -> Tuple[str, str]:
    if number_of_subdirs < 1 or number_of_subdirs > 15:
        raise ValueError("Number of subdirs must be between 1 and 15")
    parts: List[str] = []
    for _ in range(number_of_subdirs):
        parts.append(md5sum[:2])
        md5sum = md5sum[2:]
    return os.path.join(parts[0], *parts[1:]), md5sum


def refcache_create(
    fasta_files: Iterable[str],
    refcache_dir: str,
    number_of_subdirs: int = DEFAULT_NUMBER_OF_SUBDIRS
):
    with mmap.mmap(
            fileno=-1,
            length=MAX_MEMORY_SIZE,
            access=mmap.ACCESS_WRITE | mmap.ACCESS_READ
    ) as memory_spool:
        for fasta_file in fasta_files:
            logging.info(f"Reading {fasta_file} ...")
            for contig_file, md5sum, name in fasta_to_contig_files_and_md5sums(
                    fasta_file, memory_spool, refcache_dir):
                dir_name, file_name = md5sum_to_dir_and_file(
                    md5sum, number_of_subdirs)
                store_dir = os.path.join(refcache_dir, dir_name)
                store_file = os.path.join(store_dir, file_name)
                if os.path.exists(store_file):
                    logging.info(f"Already exists: {md5sum} {name}")
                    continue
                os.makedirs(store_dir, mode=0o775, exist_ok=True)
                contig_file.save_to(store_file)
                os.chmod(store_file, 0o664)
                logging.info(f"{store_file} {name}")


def file_is_fasta(filepath) -> bool:
    with open(filepath, "rb") as f:
        try:
            data = f.read(4096).decode("utf-8")
        except UnicodeDecodeError:
            return False
        if len(data) == 0:
            return False
        if data[0] != ">":
            return False
        first_newline_pos = data.find("\n")
        if first_newline_pos == -1:
            return False
        second_newline_pos = data.find("\n", first_newline_pos + 1)
        if second_newline_pos == -1:
            second_newline_pos = len(data)
        sequence_line = data[first_newline_pos + 1 : second_newline_pos]
        return set(sequence_line).issubset(set("ACGTMRWSYKVHDBNacgtmrwsykvhdbn"))


def find_fasta_files(find_dir) -> Iterator[str]:
    for entry in os.scandir(find_dir):  # type: os.DirEntry[str]
        if entry.is_symlink():
            continue
        if entry.is_file():
            if file_is_fasta(entry.path):
                yield entry.path
        if entry.is_dir():
            yield from find_fasta_files(entry.path)


def main():
    parser = argparse.ArgumentParser(
        "Import references into a cram reference cache from fasta files."
    )
    parser.usage = (
        f"{sys.argv[0]} -root <dir> [-subdirs <n>] input1.fasta ...\n       "
        f"{sys.argv[0]} -root <dir> [-subdirs <n>] -find <dir>\n")

    input_files = parser.add_argument_group("Input files (mutually exclusive)")
    input_files.add_argument(
        "fasta",
        nargs="*",
        help="Fasta files to include in the reference cache."
    )
    input_files.add_argument(
        "-find",
        metavar="<dir>",
        help="Directory to recursively search for fasta files."
    )
    parser.add_argument(
        "-root",
        required=True,
        metavar="<dir>",
        help="Directory to create the refcache dir.",
    )
    parser.add_argument(
        "-subdirs",
        metavar="<n>",
        type=int,
        default=DEFAULT_NUMBER_OF_SUBDIRS,
        help="The number of subdirectories in the cache."
    )
    args = parser.parse_args()
    if not (args.find or args.fasta):
        raise ValueError("At least one fasta file or a directory with '-find' "
                         "needs to be supplied.")
    if args.find and args.fasta:
        raise ValueError(
            "-find and positional fasta files are mutually exclusive.")
    if args.find:
        fasta_files = find_fasta_files(args.find)
    else:
        fasta_files = args.fasta
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "{message}",
        style="{")
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    refcache_create(
        fasta_files=fasta_files,
        refcache_dir=args.root,
        number_of_subdirs=args.subdirs,
    )
    logging.info(
        f"\nUse environment REF_CACHE=$root_dir{'/%2s' * args.subdirs}/%s for "
        f"accessing these files."
    )
    logging.info(
        "See also "
        "https://www.htslib.org/workflow/cram#the-ref_path-and-ref_cache for\n"
        "further information."
    )


if __name__ == "__main__":
    main()
