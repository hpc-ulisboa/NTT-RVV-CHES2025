from m5.SimObject import SimObject
from m5.params import *
from m5.defines import buildEnv
from m5.objects.FuncUnit import *
from m5.objects.FuncUnitConfig import *
from m5.objects.FUPool import *

from m5.objects.FuncUnit import *

from argparse import Namespace

PIPELINED_FPU = True

# copied from configs/common/cores/arm/O3_ARM_v7a.py
class CustomSIMD(FUDesc):
    opList = [
        OpDesc(opClass="SimdAdd", opLat=2),
        OpDesc(opClass="SimdAddAcc", opLat=2),
        OpDesc(opClass="SimdAlu", opLat=4),
        OpDesc(opClass="SimdCmp", opLat=4),
        OpDesc(opClass="SimdCvt", opLat=3),
        OpDesc(opClass="SimdMisc", opLat=3),
        OpDesc(opClass="SimdMult", opLat=5),
        OpDesc(opClass="SimdMultAcc", opLat=5),
        OpDesc(opClass="SimdMatMultAcc", opLat=5),
        OpDesc(opClass="SimdShift", opLat=3),
        OpDesc(opClass="SimdShiftAcc", opLat=3),
        OpDesc(opClass="SimdSqrt", opLat=9),
        OpDesc(opClass="SimdFloatAdd", opLat=2),
        OpDesc(opClass="SimdFloatAlu", opLat=2),
        OpDesc(opClass="SimdFloatCmp", opLat=2),
        OpDesc(opClass="SimdFloatCvt", opLat=3),
        OpDesc(opClass="SimdFloatDiv", opLat=3),
        OpDesc(opClass="SimdFloatMisc", opLat=3),
        OpDesc(opClass="SimdFloatMult", opLat=3),
        OpDesc(opClass="SimdFloatMultAcc", opLat=5),
        OpDesc(opClass="SimdFloatMatMultAcc", opLat=5),
        OpDesc(opClass="SimdFloatSqrt", opLat=9),

        OpDesc(opClass="SimdReduceAdd"),
        OpDesc(opClass="SimdReduceAlu"),
        OpDesc(opClass="SimdReduceCmp"),
        OpDesc(opClass="SimdFloatReduceAdd"),
        OpDesc(opClass="SimdFloatReduceCmp"),
        OpDesc(opClass="SimdExt"),
        OpDesc(opClass="SimdFloatExt"),
        OpDesc(opClass="SimdConfig"),
    ]
    count = 2

    def __init__(self, opts=None):
        super().__init__()
        if not opts or not opts.simd_fus:
            return
        self.count = opts.simd_fus

class CustomSIMDLdSt(FUDesc):
    opList = [
        OpDesc(opClass="SimdUnitStrideLoad"),
        OpDesc(opClass="SimdUnitStrideStore"),
        OpDesc(opClass="SimdUnitStrideMaskLoad"),
        OpDesc(opClass="SimdUnitStrideMaskStore"),
        OpDesc(opClass="SimdStridedLoad"),
        OpDesc(opClass="SimdStridedStore"),
        OpDesc(opClass="SimdIndexedLoad"),
        OpDesc(opClass="SimdIndexedStore"),
        OpDesc(opClass="SimdUnitStrideFaultOnlyFirstLoad"),
        OpDesc(opClass="SimdWholeRegisterLoad"),
        OpDesc(opClass="SimdWholeRegisterStore"),
    ]
    count = 8

    def __init__(self, opts=None):
        super().__init__()
        if not opts or not opts.simd_ldst:
            return
        self.count = opts.simd_ldst


class CustomFPU(FUDesc):
    opList = [
        OpDesc(opClass="FloatAdd", opLat=2),
        OpDesc(opClass="FloatCmp", opLat=2),
        OpDesc(opClass="FloatCvt", opLat=2),
        OpDesc(opClass="FloatMult", opLat=2), #4
        OpDesc(opClass="FloatMultAcc", opLat=3), #5
        OpDesc(opClass="FloatDiv", opLat=12, pipelined=False),
        OpDesc(opClass="FloatSqrt", opLat=24, pipelined=False),
        OpDesc(opClass="FloatMisc", opLat=3),
    ]
    count = 2

class CustomLdStU(FUDesc):
    opList = [
        OpDesc(opClass="MemRead"),
        OpDesc(opClass="MemWrite"),
        OpDesc(opClass="FloatMemRead"),
        OpDesc(opClass="FloatMemWrite"),
    ]
    count = 2

class ALUMul(FUDesc):
    opList = [
        OpDesc(opClass="IntAlu"),
        OpDesc(opClass="IntMult", opLat=3),
        OpDesc(opClass="IntDiv", opLat=20, pipelined=False),
    ]
    count = 1

class AluBr(FUDesc):
    opList = [
        OpDesc(opClass="IntAlu"),
    ]

    count = 1

class CustomFUPool(FUPool):
    def __init__(self, opts: Namespace = None):
        super().__init__()

        self.FUList = [
            ALUMul(),
            AluBr(),
            CustomFPU(),
            CustomLdStU(),
            CustomSIMD(opts),
            CustomSIMDLdSt(opts),
            IprPort(),
        ]