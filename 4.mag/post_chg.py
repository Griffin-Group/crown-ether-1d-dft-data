#!/usr/bin/env python3
print("Starting Python")
import numpy as np
import sys
# Indicies
Ba = np.r_[0:3]
Mn = np.r_[3:6]
H = np.r_[6:78]
C = np.r_[78:114]
Br = np.r_[114:126]
O = np.r_[126:144]
indices = dict(Ba=Ba, Mn=Mn, H=H, C=C, Br=Br, O=O)
neutral = dict(Ba=10, Mn=7, H=1, C=4, Br=7, O=6)
# Load file
chg = np.genfromtxt(sys.argv[1], skip_header=2, usecols=(4,), max_rows=144)
# Charges
for k in indices:
    if sum(chg) < 450:
        # Already a neutral system (e.g. spin density)
        print(f"{k} charge: {sum(chg[indices[k]])}")
    else:
        # Reference point is zero electrons, not zero charge.
        print(f"{k} charge: {neutral[k] * len(indices[k]) - sum(chg[indices[k]])}")
print("Finished Python")
