import os
os.environ['OMP_NUM_THREADS'] = '1'
print(f'{os.environ["OMP_NUM_THREADS"] = }')