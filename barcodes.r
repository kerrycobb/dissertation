library(DNABarcodes)

i7 <- c("CGATAGAG", "TTCGTTGG", "TGGAGAGT", "TCAGACGA", "GACGAATG", "CATGAGGA", "CGGTTGTT", "TCCGTATG", "TGTGGTAC", "AGAACGAG", "CTTCGTTC", "CCAATAGG")
nhei <- c("CCGAATG", "TTAGGCAG", "AACTCGTCG", "GGTCTACGTG", "GATACCG", "AGCGTTGG", "CTGCAACTG", "TCATGGTCAG")
ecori <- c("CTAACGT", "TCGGTACT", "GATCGTTGT", "AGCTACACTT", "ACGCATT", "GTATGCAT", "CACATGTCT", "TGTGCACGAT", "GCATCAT", "ATGCTGTT", "CATGACCTT", "TGCAGTGAGT")

barcode.set.distances(i7, metric="hamming")
barcode.set.distances(i7, metric="levenshtein")
                      
