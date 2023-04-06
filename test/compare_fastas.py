from Bio import SeqIO
import sys
rev = str.maketrans("ACGTacgt", "TGCAtgca")
def reverse_complement(seq: str) -> str:
    return seq.translate(rev)[::-1]
def canonical(sequence: str):
    """Returns the smallest value between a sequence and its reverse complement

    Args:
        sequence (str): sequence

    Returns:
        (str): mallest value between a sequence and its reverse complement
    """
    return min(reverse_complement(sequence), sequence)

def read_canonical(path_fasta):
    canonical_l=[]
    with open(path_fasta,"r") as f:
        for record in SeqIO.parse(f,"fasta"):
            canonical_l.append(canonical(str(record.seq)))
    return set(canonical_l)

def compare_two_fasta(path_fasta1,path_fasta2):
    return read_canonical(path_fasta1)==read_canonical(path_fasta2)
def help():
    print("python compare_fastas.py path_ccdbgupdater_fasta path_bifrost_fasta")
    exit()
if __name__=="__main__":
    if len(sys.argv)==1:
        help()
    path1=sys.argv[1]
    path2=sys.argv[2]
    print("SAME CANONICAL SEQUENCES:",compare_two_fasta(path1,path2))
