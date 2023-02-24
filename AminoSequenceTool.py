genome="ATCGGTGGCTCCGCCGCATCAGATCGAGAAACGTCTGCTCGAACTCTTGAGAGTAACCATTAACAACACGAGTCGGATTCTGACCCAAAACTTGCATCTGACGCTGGTGAGATTCGCTCATACAATGACACTTGAATCCATTCTCGTCTCGGCATTGTTTCTCACACATCTGAAATACCATCGAAGCTTTTGAAGTCCTTTCGCTTTTATTCGATTCGACATCGCCTTCGGGGTTAACATATCGTTCTTGCCCATCGATGATGAAGATTCCATTAGGGATTGAAGAGGAACTCATAACAAACATTGTTATGGTGAATAAGGCCCACTAATTAAAGCCCATACGGTACTGGGCCTAAGATTTAGATATATATTGGGCTATATGATTATTGTTGGGCCATTCTCTAGTATTCTGCATTCGTACCATAATTTACTTTTTAAGTTTTGACTTCTTTTTATTTTCTTTGGGAGCGTTTTCCTCTTTTAAATGACCGAATTTTTAGGGATCATGAAACTAGCGAGGACTATAATATTACTCTCTCATCTAAACTGAAGCAAAGGAATGACAAAATAAAGCTATCAAGCTATTGTTACCAACGCCAACATGTTGGTTCCTGAATTCTCGATATAATTTATCGAGCTTAAGGTTCAGGCTTGGTTTCTTGCACATACCCCACATTTTGATATGGAAACCGATACATCTCCGCCGATTTCATTAATGTGTTAAAGCTAGATAAGGACTCGTAAGGCATTTTCATCCATGGTAGAACTAGCATCATTAAATCTATCTCATGTCAAGACAAATGATTACGAAGGCTTTGTAGAATAACACACTAATAAAAAACAGCTTTCGGGTAAATTACGGAATCCCAAAATGCCTATATAGATAGTGAAAATTCATATGCCCATAAGCGCTAGCCAAAGCAAACTCTAGCTATTAATGTTATCACATTATAAGCATGCATATAATATACTCATGATTAAATATTCTCCATAAAGCAGCCATATAGATTTTCTCCCCCACTTTTTCTTCTTTTTTGACACCCTTTGCTACTTTCCTAAGATACCGATTACAATATGTGTATAAGATTTAAAGCCGAAATTTATCTTCTTCTTTTTTTTGTCAATTTATCTAAACCCAAATAAATTCAAAAGTAAATCCAAAACAAAATATTGTTGATATGATATCAACT"
global total_Amino_list
total_Amino_list=[]
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
    return total_mRNA

def AminoTranse(k,data):
    start=False
    startPoint=0
    endPoint=0
    Amino=[]

    try:
        for index in range(k,len(data)):
            if data[index]=='A' and data[index+1]=='U' and data[index+2]=='G':
                Amino.append('M')
                start=True
                startPoint=index+3
                break
            else:
                continue
    except:
        return total_Amino_list
    if start:
        try:
            for index in range(startPoint,len(data),3):
                if data[index]=='A':
                    if data[index+1]=='G':
                        if data[index+2]=='G' or data[index+2]=='A':
                            Amino.append('R')
                        else:
                            Amino.append('S')
                    elif data[index+1]=='A':
                        if data[index+2]=='G' or data[index+2]=='A':
                            Amino.append('K')
                        else:
                            Amino.append('N')
                    elif data[index+1]=='C':
                        Amino.append('T')
                    else:
                        if data[index+2]=='G':
                            Amino.append('M')
                        else:
                            Amino.append('I')
                elif data[index]=='C':
                    if data[index+1]=='G':
                        Amino.append('R')
                    elif data[index+1]=='A':
                        if data[index+2]=='G' or data[index+2]=='A':
                            Amino.append('Q')
                        else:
                            Amino.append('H')
                    elif data[index+1]=='C':
                        Amino.append('P')
                    elif data[index+1]=='U':
                        Amino.append('L')
                elif data[index]=='U':
                    if data[index+1]=='G':
                        if data[index+2]=='G':
                            Amino.append('W')
                        elif data[index+2]=='C' or data[index+2]=='U':
                            Amino.append('C')
                        else:
                            endPoint=index+3
                            break
                    elif data[index+1]=='A':
                        if data[index+2]=='G' or data[index+2]=='A':
                            endPoint=index+3
                            break
                        else:
                            Amino.append('Y')
                    elif data[index+1]=='C':
                        Amino.append('S')
                    else:
                        if data[index+2]=='A'or data[index+2]=='G':
                            Amino.append('L')
                        else:
                            Amino.append('F')
                else:
                    if data[index+1]=='G':
                        Amino.append('G')
                    elif data[index+1]=='A':
                        if data[index+2]=='A'or data[index+2]=='G':
                            Amino.append('E')
                        else:
                            Amino.append('D')
                    elif data[index+1]=='C':
                        Amino.append('A')
                    else:
                        Amino.append('V')
        except:
            pass
        total_Amino=''.join(Acid for Acid in Amino)
        total_Amino_list.append(total_Amino)
        AminoTranse(endPoint,data)
                

mRNA=mRNATranse(genome)
print("mRNA Sequence")
print(mRNA)
print("\n")
AminoTranse(0,mRNA)
print("Amino_Acid Sequence")
print(total_Amino_list)
