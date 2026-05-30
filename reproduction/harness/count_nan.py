import sys
import pandas as pd
import numpy as np
d = pd.read_csv(sys.argv[1])
ens = d.drop(columns=['binAges']).to_numpy()
all_nan = np.all(~np.isfinite(ens), axis=0)
print(f'{sys.argv[1]}: {all_nan.sum()}/{ens.shape[1]} members are fully NaN; '
      f'{(~all_nan).sum()} survive')
