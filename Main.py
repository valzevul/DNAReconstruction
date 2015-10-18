from SBH import DeBrujin
from Sequencer import Sequencer
# from Shotgun import Shotgun

import sys
import GenBank
import time
import C_Shotgun


sys.setrecursionlimit(100000)  # You see, De Bruijn method is kinda very recursive


def levenshteinDistance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1
    distances = range(len(s1) + 1)
    for index2, char2 in enumerate(s2):
        newDistances = [index2+1]
        for index1,char1 in enumerate(s1):
            if char1 == char2:
                newDistances.append(distances[index1])
            else:
                newDistances.append(1 + min((distances[index1],
                                             distances[index1+1],
                                             newDistances[-1])))
        distances = newDistances
    return distances[-1]


def test_shotgun(dna_slices):
    reads = []
    for block in dna_slices:  # QUICK'N'DIRTY
        for kmer in block:
            reads.append(kmer)
    answer = C_Shotgun.perform(reads)
    print "Test Shotgun method:"
    print answer
    return answer


def test_debruijn(dna, k=20):
    graph = DeBrujin.DeBruijnGraph([dna], k)
    res = graph.eulerianWalkOrCycle()
    answer = res[0] + ''.join(map(lambda x: x[-1], res[1:]))
    print "Test De Bruijn graphs method:"
    print answer
    return answer


def test_all(dna):
    start_dbj = time.clock()
    debruijn = test_debruijn(dna)  # De Bruijn Method
    end_dbj = time.clock()
    print "time elapsed: %d" % (end_dbj - start_dbj)
    sequencer = Sequencer(dna)  # Sequencer for the clone library
    dna_slices = sequencer.sequence_all()  # Clone library
    start_s = time.clock()
    shotgun = test_shotgun(dna_slices)  # Shotgun sequencing method
    end_s = time.clock()
    print "time elapsed: %d" % (end_s - start_s)

    distance = levenshteinDistance(debruijn, dna)
    print "levenshteinDistance(DB, DNA) = %s" % distance
    distance = levenshteinDistance(debruijn, shotgun)
    print "levenshteinDistance(DB, Shotgun) = %s" % distance
    distance = levenshteinDistance(dna, shotgun)
    print "levenshteinDistance(DNA, Shotgun) = %s" % distance


if __name__ == "__main__":
    choice = input("Select DNA length:\n\t1: 300 symb\n\t2: 3000 symb\n\t3: GenBank file\n")
    if choice == 1:
        dna = open("other/resources/dna.txt", 'r').readlines()[0]  # 300 symbols
        test_all(dna)
    elif choice == 2:
        dna_big = open("other/resources/dna_big.txt", 'r').readlines()[0]  # 3000 symbols
        test_all(dna_big)
    elif choice == 3:
        dna_genbank = GenBank.get_dna_from_file("U49845.gbk")
        test_all(dna_genbank)