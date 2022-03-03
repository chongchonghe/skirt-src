#!/usr/bin/env python
"""run_amr2cube.py
Create grid data from RAMSES outputs using amr2cube.f90

"""
import os
from academicpython.tools import betterRun

def runcmd(cmd):
    betterRun(cmd, prt=True, check=True)

def run_amr2cube(jobid, outs, center='c',width=0.8, lma=9, suffix='',
                 out_dir="../data/data_amr2cube",
                ):
    amr2cube = "/startrek2nb/chongchong/Sam/test_amr2cube/amr2cube_mod"
    sam_dir = "/startrek2nb/chongchong/Sam"
    os.system(f"mkdir -p {out_dir}/Job{jobid}")
    if isinstantce(center, str):
        if center == 'c':
            center = [.5, .5, .5]
    left = [c - width / 2 for c in center]
    right = [c + width / 2 for c in center]
    params = (f"-xmi {left[0]} -xma {right[0]} -ymi {left[1]} -yma {right[1]} "
              f"-zmi {left[1]} -zma {right[2]}")
    for i in outs:
        print("run_amr2cube(), doing out", i)
        for field, typ in zip(['den', 'xHII'], [1, 12]):
            denname = f"{out_dir}/Job{jobid}/out{i:02d}_{field}_l{lma}{suffix}.dat"
            if not os.path.isfile(denname):
                cmd = (f"{amr2cube} -inp {sam_dir}/Job{jobid}/output_000{i:02d}"
                       f" -out {denname} -typ {typ} -lma {lma} {params}")
                runcmd(cmd)
                print(f"{denname} created")
            else:
                print(f"{denname} exists. Skipping")
        # or ['den', 'xHII', 'pre'], [1, 12, 11]

if __name__ == "__main__":

    run_amr2cube('2.2.2', range(15, 49+1))
    run_amr2cube('3.2.2', range(14, 44+1, 2))
    run_amr2cube('4.2.1', range(14, 48+1, 1))

    #run_amr2cube('3.2.2', [38], width=0.98, suffix="_w")
    #run_amr2cube('3.2.2', [26], width=0.98, suffix="_w")
