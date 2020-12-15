#!/usr/bin/python

import random

def mir_location(df,rna_name,seed):        
    if isinstance(seed, int):
        random.seed(seed)
    df1 = df.loc[[rna_name],:]
    if df1.shape[0] == 1:
        return [df1.iloc[0,0],df1.iloc[0,1],df1.iloc[0,2]]
    else:
        random_idx = random.randint(1,df1.shape[0]-1)
        return [df1.iloc[random_idx,0],df1.iloc[random_idx,1],df1.iloc[random_idx,2]]