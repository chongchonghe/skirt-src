# run the whole skirt project for a given jobid

import sys, os
from academicpython.tools import betterRun
from run_job_step1 import run_job_out
#from mod_ski import mod_ski

sys.path.append("/startrek2nb/chongchong/Sam/coding")
from pkgs import ramses

# constants
out_fp_fmt = "/startrek2nb/chongchong/Sam/Job{jobid}/output_00001/info_00001.txt"

jobid = '3.2.2'

# decide outputs to run
r = ramses.Ramses(jobid, jobdir='/startrek2nb/chongchong/Sam')
outs_all = r.get_all_outputs()
nml="main_w.nml"
part = "part_w"
outs = outs_all

outs = [44]
outs = [26]
for out in outs:

    # make part_CK and hydro
    fn_sink = r.get_sink_path(out)
    skirt_job_dir = "../data/yorp07/run_v6"
    run_job_out(jobid, out, fn_sink, nml=nml, skirt_job_dir=skirt_job_dir,)

    # mode ski
    ski_fn = nml.replace("nml", "ski")
    ski_mod_fn = ski_fn.replace(".ski", "_1e6ph.ski")
    skirt_job_dir_injob = os.path.join(skirt_job_dir, "Job" + jobid)
    if os.path.exists(os.path.join(skirt_job_dir_injob, ski_mod_fn)):
        print(f"{skirt_job_dir_injob}/{ski_mod_fn} exists. Skipped modifying ski file")
    else:
        # old
        #mod_ski(...)
        cmd = "{cwd}/mod_ski.py {skirt_job_dir_injob}/{fi} {skirt_job_dir_injob}/{fo} {nph} {lmax} {fn_par} {skirt_job_dir}/{fn_nml} {fn_ramses_info}".format(
            skirt_job_dir=skirt_job_dir,
            skirt_job_dir_injob=skirt_job_dir_injob,
            cwd='.',
            fi=ski_fn,
            fo=ski_mod_fn,
            nph='1e6',
            lmax=9,
            fn_par=part,
            fn_nml=nml,
            fn_ramses_info=out_fp_fmt.format(jobid=jobid),
        )
        print(cmd)
        os.system(cmd)
