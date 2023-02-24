genome="ACGCGATAG"


def mRNATranse(data):
    mRNA=[]
    for index in range(len(data)):
        if data[index]=='A':
            mRNA.append('U')
        elif data[index]=='T':
            mRNA.append('A')
        elif data[index]=='G':
            mRNA.append('C')
        elif data[index]=='C':
            mRNA.append('G')
        total_mRNA=''.join(acid for acid in mRNA)
    print(total_mRNA)

mRNATranse(genome)

