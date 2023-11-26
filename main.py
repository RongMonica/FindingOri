from docx import Document
import pandas as pd
# Functions below are to find the hidden information in the Replication Origin of bacterial genome

# find the number of specific pattern in a given genome
def PatternCount(Text,Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern) +1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count

# find the specific pattern and the number it shows up a given genome
def FrequencyMap(Text,k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
    for j in range(n-k+1):
        Pattern = Text[j:j+k]
        freq[Pattern] += 1
    return freq

# find the pattern that show up the most in a given genome
def FrequentWords(Text,k):
    words = []
    freq = FrequencyMap(Text,k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words

# generate a reversed sequence of a given sequence
def Reverse(Pattern):
    rePattern = ""
    for char in Pattern:
        rePattern = char + rePattern
    return rePattern

# generate a complement sequence of a given sequence
def Complement(Pattern):
    comPattern = ""
    for char in Pattern:
        if char == 'A':
            comPattern = comPattern + 'T'
        if char == 'T':
            comPattern = comPattern + 'A'
        if char == 'G':
            comPattern = comPattern + 'C'
        if char == 'C':
            comPattern = comPattern + 'G'
    return comPattern

# generate a reversed complement sequence of a given sequence
def ReverseComplement(Pattern):
    reComPattern = ""
    reComPattern  = Reverse(Pattern)
    reComPattern = Complement(reComPattern)
    return reComPattern

# find the location of a given pattern in a given genome and return the starting site of the positions
def PatternMaching(Pattern, Genome):
    positions = []
    j = len(Genome)
    k = len(Pattern)
    for i in range(j - k + 1):
        if Genome[i:i+k] == Pattern:
            positions.append(i)
    return positions

# Below function is to calculate the number of symbol appears at position i of Genome
def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    extended_genome = Genome + Genome[0:n//2]
    array[0] = PatternCount(extended_genome[0:n//2], symbol)
    for i in range(1,n):
        array[i] = array[i-1]
        if extended_genome[i-1] == symbol:
            array[i] = array[i] - 1
        if extended_genome[i+(n//2)-1] == symbol:
            array[i] = array[i] + 1
    return array

# Below function is to calculate the the difference between the total number of occurrences of G and the total number of occurrences of C that we have encountered so far in Genome by using a skew array. This array, denoted Skew, is defined by setting Skew[i] equal to the number of occurrences of G minus the number of occurrences of C in the first i nucleotides of Genome
# this is a way to locate ori (it should be found where the skew array attains a minimum.)
def SkewArray(Genome):
    skew = [0]*(len(Genome) + 1)
    i = 0
    for char in Genome:
        if char == 'G':
            skew[i+1] = skew[i] + 1
        elif char == 'C':
            skew[i+1] = skew[i] - 1
        else:
            skew[i+1] = skew[i]
        i += 1
    return skew

# This function is to find the position of minimum G-C value, which is the location of ori
def MinimumSkew(Genome):
    skewArray = SkewArray(Genome)
    n = len(skewArray)
    miniValue = min(skewArray)
    positions = []
    for i in range(n):
        if skewArray[i] == miniValue:
            positions.append(i)
    return positions

# This function is to compare the number of differences between two sequences
def HammingDistance(p,q):
    distance = 0
    i = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance

# Find all approximate occurrences of a pattern in a string.
# DnaA box may appear with slight variations
def ApproximatePatternMatching(Text,Pattern,d):
    positions = []
    n = len(Text)
    j = len(Pattern)
    for i in range(n-j+1):
        if HammingDistance(Text[i:i+j],Pattern) <= d:
            positions.append(i)
    return positions

# The function is to find the number of times Pattern appears in Text with at most d mismatches
def ApproxiatePatternCount(Pattern,Text,d):
    count = 0
    n = len(Text)
    j = len(Pattern)
    for i in range(n-j+1):
        if HammingDistance(Text[i:i+j],Pattern) <= d:
            count += 1
    return count


#doc = Document("/Users/huangqingrong/PycharmProjects/Bioinformatics/E_coli_genome.docx")
#Genome = ""
#for paragraph in doc.paragraphs:
#    Genome += paragraph.text
#symbol = 'C'
#print(Genome[0:1000])
#genome2 = "AAAAGGGG"
#symbol2 = 'A'
#data = SymbolArray(Genome, symbol)
#df = pd.DataFrame(data, index = ['Value'])
#excel_path = "/Users/huangqingrong/PycharmProjects/Bioinformatics/findingOri.xlsx"
#df.to_excel(excel_path, index = False)
#print(type(data))
#Genome = "GATACACTTCCCGAGTAGGTACTG"
#print(SkewArray(Genome))
#p = "CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG"
#q = "ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT"
#d = 2
#print(HammingDistance(p,q))
#Genome = "atgaccgggatactgatAAAAAAAAGGGGGGGggcgtacacattagataaacgtatgaagtacgttagactcggcgccgccgacccctattttttgagcagatttagtgacctggaaaaaaaatttgagtacaaaacttttccgaataAAAAAAAAGGGGGGGatgagtatccctgggatgacttAAAAAAAAGGGGGGGtgctctcccgatttttgaatatgtaggatcattcgccagggtccgagctgagaattggatgAAAAAAAAGGGGGGGtccacgcaatcgcgaaccaacgcggacccaaaggcaagaccgataaaggagatcccttttgcggtaatgtgccgggaggctggttacgtagggaagccctaacggacttaatAAAAAAAAGGGGGGGcttataggtcaatcatgttcttgtgaatggatttAAAAAAAAGGGGGGGgaccgcttggcgcacccaaattcagtgtgggcgagcgcaacggttttggcccttgttagaggcccccgtAAAAAAAAGGGGGGGcaattatgagagagctaatctatcgcgtgcgtgttcataacttgagttAAAAAAAAGGGGGGGctggggcacatacaagaggagtcttccttatcagttaatgctgtatgacactatgtattggcccattggctaaaagcccaacttgacaaatggaagatagaatccttgcatAAAAAAAAGGGGGGGaccgaaagggaag ctggtgagcaacgacagattcttacgtgcattagctcgcttccggggatctaatagcacgaagcttAAAAAAAAGGGGGGGa"
#k = 15
#print(FrequentWords(Genome, k))

print(ReverseComplement("AAAACCCGGT"))

print(PatternMaching("ATAT","GATATATGCATATACTT"))