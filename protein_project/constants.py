AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

AA3_TO_AA1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Sec": "U",
    "Pyl": "O",
    "Ter": "*",
}

PATHOGENIC_LABELS = {
    "pathogenic",
    "likely pathogenic",
    "pathogenic/likely pathogenic",
    "likely pathogenic/pathogenic",
}

BENIGN_LABELS = {
    "benign",
    "likely benign",
    "benign/likely benign",
    "likely benign/benign",
}
