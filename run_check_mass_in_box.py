#!/startrek/chongchong/anaconda3/envs/yt/bin/python -u
"""

"""

import os
import sys
import numpy as np

sys.path.append("/startrek2nb/chongchong/Sam/coding")
from to_skirt import to_skirt
from pkgs import ramses

jobid = '2.2.2'
r = ramses.Ramses(jobid, jobdir="/startrek2nb/chongchong/Sam")
outs = r.get_all_outputs()
masses = []
n = []
skirt_job_dir = f"../data/yorp07/run_v6/Job{jobid}"
m_CK = []
n_CK = []
for out in outs:
    #masses.append(np.sum(r.get_sink_masses()))
    n.append(np.sum(r.is_sink_alive(out)))
    fn = f"{skirt_job_dir}/out{out:02d}/part_CK_corrected"
    if not os.path.isfile(fn):
        #m_CK.append(np.nan)
        n_CK.append(np.nan)
        continue
    m = np.loadtxt(fn)
    #m_CK.append(np.sum(m[:, 0]))
    n_CK.append(m.shape[0])
#np.savetxt("t.txt", np.stack((masses, np.array(m_CK)), axis=1))
np.savetxt("t.txt", np.stack((n, np.array(n_CK)), axis=1), )
