#!/home/wtk23/anaconda3/bin/python
import numpy as np
import pandas as pd
from functools import reduce
import argparse
import re

def s_all(top):
    return reduce(lambda a,b : a+b, top[1],'')

def s(strand,top):
    return reduce(lambda a,b : a+b, top[top[0] == strand][1],'')

def get_parser():
    parser = argparse.ArgumentParser(description = "sets up a pdPAINT simulation")
    parser.add_argument("--forcefile",type=str)
    parser.add_argument("--opfile",type=str)
    parser.add_argument("--top",type=str)
    return parser

class Force:
    def __init__(self,index,x,y,z,stiff):
        self.index = index
        self.x = x
        self.y = y
        self.z = z
        self.stiff = stiff

        self.force_string = "{\ntype=trap\nparticle=%s\npos0=%s,%s,%s\nstiff=%s\nrate=0\ndir=1.,0.,0.\n}\n"%(self.index,self.x,self.y,self.z,self.stiff)


class Domain:
    def __init__(self,pairs,domain_name):
        self.pairs = list(pairs)
        self.domain_name = domain_name
    def get_dist_string(self):
        dist_string = "{\norder_parameter = mindistance \nname = dist_%s\n"%(self.domain_name,)
        for index,(i,j) in enumerate(self.pairs):
            dist_string += "pair%s = %s, %s\n"%(index+1,i,j)
        dist_string += "interfaces = 1,2,4,8\n"
        dist_string += '}\n'
        return dist_string

    def get_bonds_string(self):
        dist_string = "{\norder_parameter = bond \nname = bonds_%s\n"%(self.domain_name,)
        for index,(i,j) in enumerate(self.pairs):
            dist_string += "pair%s = %s, %s\n"%(index+1,i,j)
        dist_string += '}\n'
        return dist_string


parser = get_parser()
args = parser.parse_args()

top = pd.read_csv(args.top,skiprows=1,header=None,sep=' ')

forces = pd.read_csv(args.forcefile,sep=',')
op = pd.read_csv(args.opfile,sep=',')

def apply_forces(forces,top):

    Forces = []

    for i in range(len(forces)):
        seq = forces.iloc[i]['sequence'].replace(' ','')[::-1] # flip to 3' to 5' direction
        oxdna_seq = s_all(top)

        if forces.iloc[i]['force_side'] == 5:
            index = re.search(seq,oxdna_seq).end()
        elif forces.iloc[i]['force_side'] == 3:
            index = re.search(seq,oxdna_seq).start()
            #ah, there's a bit problem here!
        breakpoint()

        Forces.append(Force(index,forces.iloc[i]['x'],
                                forces.iloc[i]['y'],
                                forces.iloc[i]['z'],
                                forces.iloc[i]['stiff']
                    ))

    with open('ext','w+') as f:
        for i in Forces:
            f.write (i.force_string)


def apply_opfile(op,top,toehold=0):
    oxdna_seq = s_all(top)
    Domains = []
    for i in range(len(op)):
        seq_A = op.iloc[i]['sequence_A'].replace(' ','')[::-1]
        seq_B = op.iloc[i]['sequence_B'].replace(' ','')[::-1]
        is_toehold = op.iloc[i]['is_toehold'] #here we need to actually work out the toehold length, but let's say its zero for now... 

        if not is_toehold: toehold = 0 #fucking hack


        #this slice on the right doesn't work
        if toehold != 0:
            seq_A_true = seq_A[:-toehold]#first we look for intrastrand bonding.
        elif toehold == 0:
            seq_A_true = seq_A
        seq_B_true = seq_B[toehold:]


        seq_A_m_vals = [m for m in re.finditer(seq_A_true,oxdna_seq)]
        
        if len(seq_A_m_vals) == 1:
            seq_A_m = seq_A_m_vals[0]

        elif len(seq_A_m_vals) == 2:
            seq_A_m = seq_A_m_vals[i]

        seq_B_m = re.search(seq_B_true,oxdna_seq) # there's only ever one B...

        pairs = zip(range(seq_A_m.start(),seq_A_m.end()),
                    range(seq_B_m.start(),  seq_B_m.end())[::-1]
        )
        
        assert(len(range(seq_A_m.start(),seq_A_m.end())) == len(range(seq_B_m.start(),  seq_B_m.end())[::-1]))

        Domains.append(Domain(pairs,op.iloc[i]['name']))

    with open("op","w+") as f:
        for domain in Domains:
            f.write(domain.get_bonds_string())
        for domain in Domains:
            f.write(domain.get_dist_string())

    return Domains


def apply_wfile(domains,top):
    #N.B: we have five distance interfaces...
    #is there a way to combine GPU and CPU simualtion here?
    wfile_string = ""
    for b1 in range(len(domains[0].pairs)+1):
        for b2 in range(len(domains[1].pairs)+1):
                for d1 in range(5):
                    for d2 in range(5):
                        wfile_string += "%s %s %s %s %s\n"%(b1,b2,d1,d2,1)
    with open("wfile","w+") as f: f.write(wfile_string)

toehold = 3
apply_forces(forces,top)
domains = apply_opfile(op,top,toehold)
apply_wfile(domains,top)

#now to define the order parameter...
