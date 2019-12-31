#! /home/wtk23/anaconda3/bin/python
import mrdna
import pandas as pd
import numpy as np 
import argparse
from mrdna_parser import make_parser
from mrdna.simulate import multiresolution_simulation as simulate

parser = make_parser()
args = parser.parse_args()

#first we set up the different segments

#The segments are: 
# D1, D1_overhang, D1_polyT, D1_toehold, S1_stem, S1_loop
# D2, D2_overhang, D2_polyT, D2_toehold, S2_stem

#first, D1, here everything is colinear

class Section:
    def __init__(self,name, start_position, stop_position,is_dsdna=True,num_bp=1):
        self.name = name
        self.start_position = start_position
        self.end_position = stop_position
        self.is_dsdna = is_dsdna
        self.num_bp = num_bp
        
    def create_helix(self):
        if self.is_dsdna:
            self.dna = mrdna.DoubleStrandedSegment(name=self.name,
                                                    start_position = self.start_position,
                                                    end_position=self.end_position,
                                                    num_bp = self.num_bp)
        else:
            self.dna = mrdna.SingleStrandedSegment(name=self.name,
                                                    start_position = self.start_position,
                                                    end_position=self.end_position,
                                                    num_nt = self.num_bp)
            #note orientation from 5' to 3'

#D1 segment
d1_num_bp = 32
DNA_LENGTH = 3.4 # Angstroms per base rise
d1_ssDNA = 15
d1_toehold_length = 3
d1_stem_length = 20 - d1_toehold_length
d1_loop_length = 6

d1 = Section('d1',np.array([0,0,0]), np.array([0,0,d1_num_bp*DNA_LENGTH]),num_bp = d1_num_bp)

d1_overhang = Section('d1_ssdna',np.array([0,0,(d1_num_bp+d1_ssDNA)*DNA_LENGTH]), np.array([0,0,d1_num_bp*DNA_LENGTH]),num_bp=d1_ssDNA,is_dsdna=False)

d1_toehold = Section('s1_toehold',np.array([0,0,(d1_num_bp+d1_ssDNA+d1_toehold_length)*DNA_LENGTH]),
                                  np.array([0,0,(d1_num_bp+d1_ssDNA)*DNA_LENGTH]),num_bp = d1_toehold_length, is_dsdna=False)

d1_stem = Section('s1_stem',np.array([0,0,(d1_num_bp+d1_ssDNA+d1_toehold_length+d1_stem_length)*DNA_LENGTH]),
                            np.array([0,0,(d1_num_bp+d1_ssDNA+d1_toehold_length)*DNA_LENGTH]),num_bp = d1_stem_length)

d1_loop = Section('s1_loop', np.array([-d1_loop_length / 2.*DNA_LENGTH,0,(d1_num_bp+d1_ssDNA+d1_toehold_length+d1_stem_length)*DNA_LENGTH]),
                             np.array([d1_loop_length / 2. *DNA_LENGTH,0,(d1_num_bp+d1_ssDNA+d1_toehold_length+d1_stem_length)*DNA_LENGTH]),
                             num_bp = d1_loop_length)

#D5 segment
Y = 10
d5_num_bp = 28

d5 = Section('d5',np.array([0,Y,0]), np.array([0,Y,d5_num_bp*DNA_LENGTH]),num_bp = d5_num_bp)

d5_overhang = Section('d5_ssdna',np.array([0,Y,(d5_num_bp+d1_ssDNA)*DNA_LENGTH]), np.array([0,Y,d5_num_bp*DNA_LENGTH]),num_bp=d1_ssDNA,is_dsdna=False)

d5_toehold = Section('s5_toehold',np.array([0,Y,(d5_num_bp+d1_ssDNA+d1_toehold_length)*DNA_LENGTH]),
                                  np.array([0,Y,(d5_num_bp+d1_ssDNA)*DNA_LENGTH]),num_bp = d1_toehold_length, is_dsdna=False)

d5_stem = Section('s5_stem',np.array([0,Y,(d5_num_bp+d1_ssDNA+d1_toehold_length+d1_stem_length)*DNA_LENGTH]),
                            np.array([0,Y,(d5_num_bp+d1_ssDNA+d1_toehold_length)*DNA_LENGTH]),num_bp = d1_stem_length)

all_segs = [d1,d1_overhang,d1_toehold,d1_stem,d1_loop,d5,d5_overhang,d5_toehold,d5_stem]
for i in all_segs: i.create_helix()

d1.dna.connect_end5(d1_overhang.dna.end3)
d1_toehold.dna.connect_end3(d1_overhang.dna.start5)
d1_stem.dna.connect_start3(d1_toehold.dna.start5)
d1_loop.dna.connect_start5(d1_stem.dna.end3)
d1_loop.dna.connect_end3(d1_stem.dna.end5)

segs_list = [x.dna for x in all_segs]

#breakpoint()

model = mrdna.SegmentModel(
                    segs_list,
                    local_twist = True,
                    dimensions=(5000,5000,5000),
                )


model.set_sequence('T'*100000)
prefix = "DNA" 

run_args = dict(
            model = model,
            output_name = prefix,
            job_id = "job-" + prefix,
            directory = args.directory,
            gpu = args.gpu,
            minimization_output_period = int(args.output_period),
            coarse_local_twist = args.coarse_local_twist,
            fix_linking_number = args.fix_linking_number,
            coarse_output_period = int(args.output_period),
            fine_output_period = int(args.output_period),
            minimization_steps = 0, # int(args.minimization_steps),
            coarse_steps = int(args.coarse_steps),
            fine_steps = int(args.fine_steps),
            backbone_scale = args.backbone_scale,
            oxdna_steps = args.oxdna_steps,
            oxdna_output_period = args.oxdna_output_period
        )


simulate( **run_args )
