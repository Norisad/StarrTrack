#create large DHS region
import sys
import pandas as pd
import numpy as np

# input=sys.argv[1]
# output=sys.argv[2]


#input="/home/sadouni/mount/sacapus/silencer_pipeline/pgk_from_bam/clone_list.1bp.PGK.r1.DHS.bed"
def enlarge_overlapping_region(input,output):
    fi=open(input,"r")
    fo=open(output,"w")
    df = pd.read_table(fi, delimiter='\t',header=None,names=["chr","start","end","point","score","strand","cdna_count","lib_count","region_type","region_id"])
    df1 = (df.groupby('region_id', as_index=False)
         .agg({'chr':'first', 'start':'min', 'end':'max','region_type':'first'})
         [['chr','start','end','region_type','region_id']])
    df1 = df1[df1.region_id != "."]
    df1.to_csv(fo,index=False, sep='\t')

    return(df1)

def make_bins(bins,input,output):
    fi=open(input,"r")
    fo=open(output,"w")
    first_lines=fi.readlines()
    lines = first_lines[1:] #Remove the header
    fo.write("#chr\tstart\tend\tregion_type\tregion_id\tsubregion_id\n") #write new header
    for line in lines: # by line
        elem=line.split("\n")[0].split("\t") #separate the tabulated element
        start=int(elem[1]) # assign start end and region
        end=int(elem[2])
        region=elem[4]
        i=1 #subregion ID
        while start < end : #to dont have greater start than end
            fo.write(('%s\t%i\t') %(elem[0],start)) #write the first start
            start = start + int(bins) #add N bins
            subregion=region+"_"+str(i) #create subregion ID
            if start > end: # to don't have bigger regions than original
                start=end
            i=i+1
            fo.write(('%i\t%s\t%s\t%s\n') %(start,elem[3],elem[4],subregion))
