from random import randint


class Helpers:
    def __init__(self):
        pass

    @staticmethod
    def random_step(size):
        return randint(0, int(size * 0.20))  # More than 20% could lead to errors

class Sequencer:
    def __init__(self, dna, k=10, i=500):
        self.dna = dna
        self.dna_len = len(dna)
        self.i = i
        self.k = k

    def __str__(self):
        return "%s: i = %i, k = %i" % (self.dna, self.i, self.k)

    def sequence_once(self):
        sls = []
        count = 0
        count += Helpers.random_step(self.dna_len)
        while count + self.k <= self.dna_len:
            sls.append(self.dna[count:count + self.k])
            count += self.k + Helpers.random_step(self.dna_len)
        return sls

    def sequence_all(self):
        results = []
        for _ in xrange(self.i):
            results.append(self.sequence_once())
        return results

if __name__ == "__main__":
    pass
