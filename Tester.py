from SBH import DeBrujin
from Shotgun import Shotgun
import C_Shotgun
import C_DeBrujin

import GenBank
from Sequencer import Sequencer
from os import listdir
from os.path import isfile, join
from timeit import Timer
from Main import levenshteinDistance


class Stats():
    def __init__(self, name, sbh_time, sbhfast_time, s_time, sfast_time, length, s_acc, sfast_acc, sbh_acc, sbhfast_acc):
        self.name = name
        self.sbh_time = sbh_time
        self.s_time = s_time
        self.length = length
        self.s_acc = s_acc
        self.sbh_acc = sbh_acc
        self.sbhfast_time = sbhfast_time
        self.sfast_time = sfast_time
        self.sfast_acc = sfast_acc
        self.sbhfast_acc = sbhfast_acc


def get_list_of_genbank_files(mypath):
    gbks = [mypath + f for f in listdir(mypath) if isfile(join(mypath, f))]
    return gbks


def sbh(k=20):
    graph = DeBrujin.DeBruijnGraph([dna], k)
    res = graph.eulerianWalkOrCycle()
    answer = res[0] + ''.join(map(lambda x: x[-1], res[1:]))

    ld = float(levenshteinDistance(dna, answer))
    l = float(len(dna))
    s[-1].sbh_acc += ((l - ld) / l) * 100


def sbhfast(k=20):
    graph = C_DeBrujin.DeBruijnGraph([dna], k)
    res = graph.eulerianWalkOrCycle()
    answer = res[0] + ''.join(map(lambda x: x[-1], res[1:]))

    ld = float(levenshteinDistance(dna, answer))
    l = float(len(dna))
    s[-1].sbhfast_acc += ((l - ld) / l) * 100


def shotgun():
    sequencer = Sequencer(dna)  # Sequencer for the clone library
    dna_slices = sequencer.sequence_all()
    reads = []
    for block in dna_slices:  # QUICK'N'DIRTY
        for kmer in block:
            reads.append(kmer)
    ld = float(levenshteinDistance(dna, Shotgun.perform(reads)))
    l = float(len(dna))
    s[-1].s_acc += ((l - ld) / l) * 100


def shotgunfast():
    sequencer = Sequencer(dna)  # Sequencer for the clone library
    dna_slices = sequencer.sequence_all()
    reads = []
    for block in dna_slices:  # QUICK'N'DIRTY
        for kmer in block:
            reads.append(kmer)
    ld = float(levenshteinDistance(dna, C_Shotgun.perform(reads)))
    l = float(len(dna))
    s[-1].sfast_acc += ((l - ld) / l) * 100


# short = get_list_of_genbank_files("gbk/range/")
# short = get_list_of_genbank_files("gbk/huge/")
short = get_list_of_genbank_files("gbk/all_gb/")
print "time/accuracy"
s = []
for f in short:
    stats = Stats("", 0, 0, 0, 0, 0, 0, 0, 0, 0)
    s.append(stats)
    name, l, dna = GenBank.get_dna_from_file(f)
    print l + "\t" + name
    stats.length = len(dna)
    #
    # t = Timer('Tester.sbh()', 'import Tester')
    # reps = 5
    # s[-1].sbh_time = sum(t.repeat(repeat=reps, number=1)) / reps
    # s[-1].sbh_acc /= reps
    # print str(s[-1].sbh_time) + "\t" + str(s[-1].sbh_acc) + "\t" + str(l)

    # t = Timer('Tester.sbhfast()', 'import Tester')
    # reps = 5
    # s[-1].sbhfast_time = sum(t.repeat(repeat=reps, number=1)) / reps
    # s[-1].sbhfast_acc /= reps
    # print str(s[-1].sbhfast_time) + "\t" + str(s[-1].sbhfast_acc) + "\t" + str(l)

    # t = Timer('Tester.shotgun()', 'import Tester')
    # reps = 5
    # s[-1].s_time = sum(t.repeat(repeat=reps, number=1)) / reps
    # s[-1].s_acc /= reps
    # print str(s[-1].s_time) + "\t" + str(s[-1].s_acc) + "\t" + str(l)
    #
    # t = Timer('Tester.shotgunfast()', 'import Tester')
    # reps = 5
    # s[-1].sfast_time = sum(t.repeat(repeat=reps, number=1)) / reps
    # s[-1].sfast_acc /= reps
    # print str(s[-1].sfast_time) + "\t" + str(s[-1].sfast_acc) + "\t" + str(l)
