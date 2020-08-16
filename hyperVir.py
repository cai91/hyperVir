import numpy as np
import matplotlib.pyplot as plt
import subprocess
import scipy.stats
import sys
from scipy.signal import savgol_filter


file=sys.argv[1]

def mutPlot_noGap(seqs):
    sym=['A','I','L','M','V','F','W','Y','N','C','Q','S','T','D','E','R','H','K','G','P','X']
    ent=[]
    mL=0
    # Getting longest sequence
    for s in seqs:
        if len(s)>mL:
            mL=len(s)
    # Counting different characters per position
    for idx in range(mL):
        loc=[0]*len(sym)
        for s in seqs:
            try:
                loc[sym.index(s[idx])]+=1
            except:
                loc[sym.index('X')]+=1  # Any other symbol is pooled into an 'X' (including indels)
        loc=np.array(loc)/len(seqs)
        ent.append(scipy.stats.entropy(loc))
    return ent

def getAA_variation(file):
    # Generating MSA with MAFFT

    print('Generating alignment...')
    with open(file.split('.fa')[0]+'.out','w') as outFile:
        subprocess.call(['mafft','--auto',file],stdout=outFile)

    # Collecting MAFFT aligned protein sequences
    seqs_fibre=[]
    seq=''
    with open(file.split('.fa')[0]+'.out') as inFile: 
        for line in inFile:
            if line[0]=='>':
                if seq!='':
                    seqs_fibre.append(seq)
                    seq=''
            else:
                seq+=line.strip()
        if seq!='':
            seqs_fibre.append(seq)

    # Computing amino acid diversity per position (Shannon's entropy)
    y=mutPlot_noGap(seqs_fibre)
    
    # Smoothing the signal
    savF_w=51 # filter size
    if len(y)<savF_w:
        print(-1) # Proteins are too short for filter size

    yhat = savgol_filter(y, savF_w, 3)
    if sum(yhat)==0.0:
        print(-1)   # There is no diversity (all sequences are identical)

    with open(file.split('.')[0]+'_Aminoacid_diversity.txt','w') as outFile:
        for v in yhat:
            outFile.write(str(v)+'\n')
        
    plt.plot(yhat)
    plt.xlabel('Alignment position')
    plt.ylabel('Amino acid diversity')
    plt.savefig(file.split('.')[0]+'_aminoacid_diversity.png')

def main():
	getAA_variation(file)

main()