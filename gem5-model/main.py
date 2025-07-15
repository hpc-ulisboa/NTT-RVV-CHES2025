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

""" This file creates a single CPU and a two-level cache system.
This script takes a single parameter which specifies a binary to execute.
If none is provided it executes 'hello' by default (mostly used for testing)

See Part 1, Chapter 3: Adding cache to the configuration script in the
learning_gem5 book for more information about this script.
This file exports options for the L1 I/D and L2 cache sizes.

IMPORTANT: If you modify this file, it's likely that the Learning gem5 book
           also needs to be updated. For now, email Jason <power.jg@gmail.com>

"""

import os
import argparse
# import the m5 (gem5) library created when gem5 is built
import m5
# import all of the SimObjects
from m5.objects import *

from caches import L1DCache, L1ICache, L2Cache, L3Cache
from func_units import CustomFUPool

parser = argparse.ArgumentParser(
    prog="gem5 model",
    description="Basic gem5 model with input parameters"
)

parser.add_argument("binary",   help="Binary to simulate")
parser.add_argument("args",     help="Remaining arguments", nargs='*') # 0 to N arguments
parser.add_argument("--simd-ldst")
parser.add_argument("--simd-fus")
parser.add_argument("--vlen", default=512)

args = parser.parse_args()

# create the system
system = System()

# Set the clock frequency
system.clk_domain = SrcClockDomain()
system.clk_domain.clock = "1GHz"
system.clk_domain.voltage_domain = VoltageDomain()

# Memory settings
system.mem_mode = "timing"  # Use timing accesses
system.mem_ranges = [AddrRange("512MiB")]  # Create an address range
# Create a memory bus
system.membus = SystemXBar()

# Create out-of-order RISC-V CPU
system.cpu = RiscvO3CPU()
system.cpu.isa = RiscvISA(vlen=args.vlen)
# Zen 5-like configuration
system.cpu.numPhysVecRegs = 384
system.cpu.numROBEntries = 630
system.cpu.issueWidth = 12 # should be 16, but is hard-capped at 12
system.cpu.renameWidth = 12 # should be 14
system.cpu.dispatchWidth = 8
system.cpu.squashWidth = 32
system.cpu.mmu.dtb.size = 96
system.cpu.LQEntries = 630
system.cpu.SQEntries = 100
system.cpu.numIQEntries = 630
system.cache_line_size = 128

# Functional units
system.cpu.fuPool = CustomFUPool(args)

# Create an L1 instruction and data cache
system.cpu.icache = L1ICache()
system.cpu.dcache = L1DCache()
# Connect the instruction and data caches to the CPU
system.cpu.icache.cpu_side = system.cpu.icache_port
system.cpu.dcache.cpu_side = system.cpu.dcache_port

# Create an L2 Cache and the respective bus
system.l2cache = L2Cache()
system.l2cache.xbar = L2XBar()

# Connect the L1 to the L2
system.cpu.icache.mem_side = system.l2cache.xbar.cpu_side_ports
system.cpu.dcache.mem_side = system.l2cache.xbar.cpu_side_ports
system.l2cache.cpu_side = system.l2cache.xbar.mem_side_ports

# Create an L3 Cache and the respective bus
system.l3cache = L3Cache()
system.l3cache.xbar = L2XBar()

# Connect the L2 to the L3
system.l2cache.mem_side = system.l3cache.xbar.cpu_side_ports
system.l3cache.cpu_side = system.l3cache.xbar.mem_side_ports

# Connect the L3 to the memory bus
system.l3cache.mem_side = system.membus.cpu_side_ports
# Connect the system up to the membus
system.system_port = system.membus.cpu_side_ports

# create the interrupt controller for the CPU
system.cpu.createInterruptController()

# Create a DDR3 memory controller
system.mem_ctrl = MemCtrl()
#system.mem_ctrl.dram = DDR5_4400_4x8()
system.mem_ctrl.dram = DDR3_1600_8x8()
system.mem_ctrl.dram.range = system.mem_ranges[0]
system.mem_ctrl.port = system.membus.mem_side_ports

system.workload = SEWorkload.init_compatible(args.binary)

# Create a process
process = Process()
process.cmd = [args.binary, *args.args]
process.output = "stdout.txt"
system.cpu.workload = process
system.cpu.createThreads()

# set up the root SimObject and start the simulation
root = Root(full_system=False, system=system)
# instantiate all of the objects we've created above
m5.instantiate()


print(f"Beginning simulation!")
while True:
    exit_cause = m5.simulate().getCause()

    if exit_cause == "checkpoint":
        continue
    else:
        break