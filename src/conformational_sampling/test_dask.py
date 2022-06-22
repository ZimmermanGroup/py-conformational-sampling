# %%
from IPython.display import display

from dask.distributed import Client, progress
from dask_jobqueue import SLURMCluster

import time
import random

def inc(x):
    time.sleep(random.random())
    return x + 1

def double(x):
    time.sleep(1)
    return 2 * x

def add(x, y):
    time.sleep(random.random())
    return x + y

def test_local():
    with Client(n_workers=2) as client:
        future = client.submit(double, 3)
        display(client)
        display(future.result(timeout=5))

def test_slurm():
    with SLURMCluster(
        queue='zimintel',
        cores=1,
        memory='1GB',
        processes=1,
        walltime='00:10:00',
    ) as cluster:
        print(cluster.job_script())
        cluster.scale(jobs=2)
        with Client(cluster) as client:
            display(client)
            future = client.submit(double, 3)
            display(future.result(timeout=5))
if __name__ == '__main__':
    test_slurm()
    
# %%
# from dask_jobqueue import SLURMCluster
# from dask.distributed import Client
# with SLURMCluster(
#     queue='zimintel',
#     cores=1,
#     memory='1GB',
#     processes=1
# ) as cluster, Client(cluster) as client:
#     cluster.scale(jobs=2)
    # futures = client.map(execute_xtb, range(len(complexes)), complexes)
    # mol = futures[0].result(timeout=10).to_rdkit_mol()
# Chem.SanitizeMol(mol)
# display(nglview.show_rdkit(mol))
