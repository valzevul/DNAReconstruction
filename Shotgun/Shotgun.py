import sys

def combine(s1, s2, overlap):
    return ''.join([s1, s2[overlap:]])


def remove_contained_reads(reads):
    if len(reads) > 1:
        for i in range(len(reads) - 1):
            for j in range(i + 1, len(reads)):
                if reads[i] in reads[j]:
                    reads[i] = "JUNK"
                elif reads[j] in reads[i]:
                    reads[j] = "JUNK"
    return [read for read in reads if read != "JUNK"]


def get_max_overlap(s1, s2):
    overlap = ""
    S1thenS2 = None
    for i in range(1, (len(s2) + 1)):
        if s1[:i] == s2[-i:] and len(s1[:i]) > len(overlap):
            overlap = s1[:i]
            S1thenS2 = False

    for j in range(1, (len(s1) - len(s2) + 1)):
        subS1 = s1[j:len(s2) + j]
        if subS1 == s2 and len(subS1) > len(overlap):
            overlap = subS1
            S1thenS2 = True

    for k in range(1, (len(s2) + 1)):
        subS2 = s2[:-k]
        subS1 = s1[-(len(s2) - k):]
        if subS1 == subS2 and len(subS1) > len(overlap):
            overlap = subS1
            S1thenS2 = True

    return (overlap, S1thenS2)


def reconstructDNA(reads, threshold):
    overlap_dict = {}
    combined_reads = []
    closed_set = set()

    for i in range(len(reads) - 1):
        for j in range(i + 1, len(reads)):
            max_overlap = get_max_overlap(reads[i], reads[j])
            if len(max_overlap[0]) in overlap_dict:
                overlap_dict[len(max_overlap[0])] \
                    .append((reads[i], reads[j]) if max_overlap[1] else (reads[j], reads[i]))
            else:
                overlap_dict[len(max_overlap[0])] = \
                    [(reads[i], reads[j])] if max_overlap[1] else [(reads[j], reads[i])]

    overlaps = list(overlap_dict.keys())
    overlaps.sort(reverse=True)

    for overlap in overlaps:
        for s1, s2 in overlap_dict[overlap]:
            if overlap >= threshold and (s1 not in closed_set) and (s2 not in closed_set):
                closed_set.add(s1)
                closed_set.add(s2)
                combined_reads.append(combine(s1, s2, overlap))

    for item in reads:
        if item not in closed_set:
            combined_reads.append(item)

    return combined_reads


def perform(reads):
    prev, threshold = -1, 10
    while len(reads) != 1:
        if len(reads) == prev:
            threshold -= 1
        prev = len(reads)
        reads = reconstructDNA(remove_contained_reads(reads), threshold)
    return reads[0]


if __name__ == "__main__":
    pass
