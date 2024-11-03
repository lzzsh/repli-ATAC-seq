import pysam
infile = pysam.AlignmentFile("-", "rb")
for sam in infile:
    if sam.is_reverse :
        print(sam.reference_name +'\t'+str(sam.aend-4))
    else:
        pos2 = sam.pos + 5
        print(sam.reference_name +'\t'+str(pos2))
