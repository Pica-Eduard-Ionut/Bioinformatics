# implement an app that converts the coding region of a gene into an amino acid sequence. Use the genetic code table from moodle.

RNA = "AACUCCCUAACGGGAGUGGCCAUGCAUGCUUGUCCAAGGUAGUGGAUAGAAAACGUUAAUUCGGCGGAUAUGCAUCCACCCGCAGUGUUGAUUGCUUCAUCUGAGUCAAGCUCGCAUCCCACCGCUUACUGCUUUGCAUUCCCUGAGUAGGGCGUCCACUUGAGCUUCAUCGGAGAGUAACACGACGUAAUUAUGGAAUAAACCCCUGUGGUUGAUUCGCCGCGUUGGGUAUGCGUUAACCUGCAAGACAUUAAUCGUAGCCGUACAGUGAGAGACAGAAGUUGAUAAUGGUGUCAUGCUCUGUCCACUGCGCGAGGCCAACUGGCGCGACUGCUGCGUCAUGUGAGCUCGUGCACCGUCUGUGAUCCAGCACAUCAACCGUCUCAAGUGUGGCAACUGAGCCAGUAGUGGCUGACCUGAUCGCACACCGGCCAAUCCUUCAAUAUGUCUUGGUUGGCCGCGCAACGCUACAGUGGACUUAGUUAGACAGACACCCAACGAUUUCACUUGCAGCUUAAGAACGUUAAAUGCCGAUUUUCCGCUUGCUCAGAAAGCGUAUGCCCCCAUAACUCCAAAAGACUAAGAGGGCGCACGGCGCUGCGUAGCGGACGAUAAUGUCUAAUGCUUCUGGCCGCCUGCCUAGUAAUGCCUGGCCGCAAGCACACUUCCAGAUCGUAAAAAGACACGGUGCGUUCAGACCAGCCGGCUUCCCAUUGGGAGCGUUCAACCUAUCUGUUCUUCAGAUGCGGUGGUAGUCAGCACAACUCCCGACGAAAACCCUAGACAGUCGUACAGCUUUUCGGCCGGUGGUGCUGCGCCACAUAAUAUCAACUGUGUCAGGAAUCCAACCGCAAGGCUAUAGACUUCUGGUAGGGCGAGUUAUCCGUCACCGGCACUAACGGGCUGGUCGUUGACCAAAUCAAGCGAGCAGAGAUUGUUCCUGGAGGCGCCCGCGCCAUUGUCGGCGCCCAAGUAGAUAAGCAAUUCAAAAUGGUCGCUU"

def select(index=0,length=3):
    part = S[index:index+length]
    return part

lookupTable = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",

    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",

    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",

    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
}

def translate(rna):
    start_pos = rna.find("AUG")

    protein = []
    for i in range(start_pos, len(rna) - 2, 3):
        codon = rna[i:i+3]
        amino = lookupTable.get(codon, None)
        if amino is None or amino == "Stop":
            break
        protein.append(amino)
    return protein


protein = translate(RNA)

print("Amino acid sequence:")
print("-".join(protein))
