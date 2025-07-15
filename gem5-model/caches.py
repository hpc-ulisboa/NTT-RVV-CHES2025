# Copyright (c) 2015 Jason Power
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer;
# redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution;
# neither the name of the copyright holders nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

""" Caches with options for a simple gem5 configuration script

This file contains L1 I/D and L2 caches to be used in the simple
gem5 configuration script. It uses the SimpleOpts wrapper to set up command
line options from each individual class.
"""

import m5
from m5.objects import Cache, StridePrefetcher, TaggedPrefetcher, WriteAllocator

from argparse import Namespace

# Some specific options for caches
# For all options see src/mem/cache/BaseCache.py

class L1Cache(Cache):
    """Simple L1 Cache with default values"""
    assoc = 4
    tag_latency = 2
    data_latency = 2
    response_latency = 2
    prefetcher = StridePrefetcher(degree=1, latency=1, prefetch_on_access=True)
    write_buffers = 32
    #write_allocator = WriteAllocator()


class L1ICache(L1Cache):
    # https://www.cpubenchmark.net/cpu.php?cpu=Intel+Xeon+Gold+6312U+%40+2.40GHz&id=4606
    size = "32KiB"
    tag_latency = 1
    data_latency = 1
    response_latency = 1
    mshrs = 8
    tgts_per_mshr = 20

class L1DCache(L1Cache):
    size = "64KiB"
    mshrs = 32
    tgts_per_mshr = 32

class L2Cache(Cache):
    size = "1MiB"
    assoc = 16
    tag_latency = 5
    data_latency = 5
    response_latency = 5
    mshrs = 32
    tgts_per_mshr = 12
    prefetcher = TaggedPrefetcher(degree=8, latency=1, prefetch_on_access=True)

class L3Cache(Cache):
    size = "32MiB"
    assoc = 16
    tag_latency = 20
    data_latency = 20
    response_latency = 20
    mshrs = 32
    tgts_per_mshr = 12
    write_buffers = 16
