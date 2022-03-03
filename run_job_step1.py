#!/startrek/chongchong/anaconda3/envs/yt/bin/python -u
"""run_job_step1.py

Should run on startrek.
1. Create part_CK for all outputs of a job
2. Run RAMSKI to create hydro for all outputs

"""

import os
import sys
import f90nml
from academicpython.tools import betterRun

# sys.path.append("/startrek2nb/chongchong/Sam/coding")
sys.path.append("/startrek/chongchong/Sam/coding")
from to_skirt import to_skirt
from pkgs import ramses

def run_222():
    """Run RAMSKI in all outputs of a job"""

    #jobid = '2.2.2'
    #width = 0.8
    #exe = "/startrek/chongchong/Academic/SKIRT/RAMSKI_test/ramski"
    #skirt_job_dir = "../data/yorp07/run_v6"
    ## outs = r.get_all_outputs()
    #outs = range(15, 50)
    #for out in outs:
    #    print(f"\nDoing out {out}")

    #    # make part_CK
    #    out_dir = f"{skirt_job_dir}/out{out:02d}"
    #    os.system(f"mkdir -p {out_dir}")
    #    to_skirt(jobid=jobid,
    #             output=out,
    #             fn_out=f"{out_dir}/part_CK",
    #             family='CK',
    #             width_boxlen=width,
    #             letdie=True,
    #             skip_exist=True,
    #    )
    #    # cmd = "python /startrek2nb/chongchong/Sam/coding/to_skirt.py " +\
    #    #     f"{jobid} {out} {out_dir}/part_CK -s CK -width_boxlen {width} -die"
    #    # run_o, run_e = betterRun(cmd, check=0)
    #    # if run_e != '':
    #    #     print('Error:', run_e)
    #    #     break
    #    # print(run_o)

    #    # make hydro using RAMSKI
    #    inp = f"/startrek2nb/chongchong/Sam/Job2.2.2/output_000{out:02d}"

    #    # if not os.path.isfile(os.path.join(out_dir, 'hydro')):
    #    #     if os.path.isfile("../data/run_v4/sam_{jobid}_out{out:02d}/hydro"):
    #    #         cmd = f"ln -s ../../../run_v4/sam_{jobid}_out{out:02d}/hydro {out_dir}/hydro"
    #    #         print(cmd)
    #    #         os.system(cmd)
    #    #     else:
    #    #         cmd = f"{exe} -inp {inp} -nmlpath {skirt_job_dir}/main.nml -outdir {out_dir}"
    #    #         print(cmd)
    #    #         run_o, run_e = betterRun(cmd, check=0)
    #    #         if run_e != '':
    #    #             print(f"This is out {out}")
    #    #             print(run_e)
    #    #             break
    #    #         print(run_o)
    #    # else:
    #    #     print(f"{out_dir}/hydro exists. Skipped")

    #    nml = "main.nml" if jobid[0] != '4' else "main_l15.nml"

    #    if not os.path.isfile(os.path.join(out_dir, 'hydro')):
    #        if os.path.isfile("../data/run_v4/sam_{jobid}_out{out:02d}/hydro"):
    #            cmd = f"ln -s ../../../run_v4/sam_{jobid}_out{out:02d}/hydro {out_dir}/hydro"
    #            print(cmd)
    #            os.system(cmd)
    #        else:
    #            cmd = f"{exe} -inp {inp} -nmlpath {skirt_job_dir}/{nml} -outdir {out_dir}"
    #            print(cmd)
    #            run_o, run_e = betterRun(cmd, check=0)
    #            if run_e != '':
    #                print(f"This is out {out}")
    #                print(run_e)
    #                break
    #            print(run_o)
    #    #else:
    #    #    if os.path.islink("../data/run_v4/sam_{jobid}_out{out:02d}/hydro"):
    #    #        print(f"skip ../data/run_v4/sam_{jobid}_out{out:02d}/hydro")
    #    #    else:
    #    #        cmd = f"{exe} -inp {inp} -nmlpath {skirt_job_dir}/{nml} -outdir {out_dir}"
    #    #        print(cmd)
    #    #        run_o, run_e = betterRun(cmd, check=0)
    #    #        if run_e != '':
    #    #            print(f"This is out {out}")
    #    #            print(run_e)
    #    #            break
    #    #        print(run_o)
    return

def run_job_out(jobid, out, fn_sink, nml=None,
                skirt_job_dir="../data/yorp07/run_v6",
                part="part_CK", letdie=True,
               ):
    """Run RAMSKI in all outputs of a job. Will do the following for the given output:
        - Produce particle data for SKIRT run
        - Produce hydro data for SKIRT run using ramski

    Args:
        nml (string): Default: main-jobi{jobid}.nml
    """

    exe = "/startrek/chongchong/Academic/SKIRT/RAMSKI_test/ramski"

    if nml is None:
        # nml = "main.nml" if jobid[0] != '4' else "main_l15.nml"
        nml = f"main-job{jobid}.nml"
    fn_nml = f"{skirt_job_dir}/{nml}"
    the_nml = f90nml.read(fn_nml)

    # get width
    params = the_nml["PARAMS"]
    msg = ("{xi}max and {xi}min does not center around 0.5. This is not"
           "consistant with my to_skirt.py file, so I will stop here.")
    iscenter = True
    for xi in ['x', 'y', 'z']:
        iscenter = iscenter and (abs(float(params[xi + "max"]) + float(params[xi + "min"]) - 1) \
                < 1e-10)
    if iscenter:
        center = None   # center is the box center
    else:
        center = [(params["xmax"] + params["xmin"])/2,
                  (params["ymax"] + params["ymin"])/2,
                  (params["zmax"] + params["zmin"])/2]
    width = float(params[xi + "max"]) - float(params[xi + "min"])

    print(f"\nDoing out {out}")
    if os.stat(fn_sink).st_size == 0:
        print(fn_sink, 'is empty. Skipped.')
        return

    # make part_CK
    out_dir = f"{skirt_job_dir}/Job{jobid}/out{out:02d}"
    os.system(f"mkdir -p {out_dir}")
    # name_prt = params["name_prt"]
    # assert name_prt != "trash", "name_prt is trash. Please change it!"
    name_prt = part
    fn_out = os.path.join(out_dir, name_prt)
    if not __debug__:
        fn_out += "_debug"
    to_skirt(jobid=jobid,
             output=out,
             fn_out=fn_out,
             center=center,
             family='CK',
             width_boxlen=width,
             letdie=letdie,
             skip_exist=True,
    )
    if not __debug__:
        return

    # make hydro using RAMSKI
    hydro_fp = os.path.join(out_dir, the_nml["PARAMS"]["name_hdr"])
    inp = f"/startrek2nb/chongchong/Sam/Job{jobid}/output_000{out:02d}"
    if not os.path.isfile(hydro_fp):
        cmd = f"{exe} -inp {inp} -nmlpath {fn_nml} -outdir {out_dir}"
        print(cmd)
        run_o, run_e = betterRun(cmd, check=0)
        if run_e != '':
            print(f"This is out {out}")
            print(run_e)
            return
        print(run_o)
    else:
        print(f"{out_dir}/hydro exists. Skipped")

def run_job(jobid, outs=None, skip=1, **kwargs):
    """
    Example:
        >>> run_job(jobid, outs, skirt_job_dir=..., part=..., letdie=0)
    """

    r = ramses.Ramses(jobid, jobdir='/startrek2nb/chongchong/Sam')
    if outs is None:
        outs = r.get_all_outputs()
    for out in outs:
        if not out % skip == 0:
            continue
        fn_sink = r.get_sink_path(out)
        run_job_out(jobid, out, fn_sink, **kwargs)

if __name__ == "__main__":

    if len(sys.argv) >= 2:
        run_job(sys.argv[1])
    else:
        print(f"Usage: {sys.argv[0]} jobid")
