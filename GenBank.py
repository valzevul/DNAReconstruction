from Bio import SeqIO


def get_dna_from_file(gb_file):
    gb_record = SeqIO.read(open(gb_file, "r"), "genbank")
    # print "Name %s, %i features" % (gb_record.name, len(gb_record.features))
    s = gb_record.seq
    # print str(len(s)) + " symbols"
    return str(gb_record.name), str(len(s)), str(s)
