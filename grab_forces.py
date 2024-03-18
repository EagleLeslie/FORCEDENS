import numpy as np
import os
from xdatcar import Xdatcar

def grab_forces(atoms):

    os.system(f"grep -A{atoms+2} 'TOTAL-FORCE (eV/Angst)' OUTCAR > forces")

if __name__ == "__main__":
    xdat = Xdatcar()

    grab_forces(xdat.totatoms)