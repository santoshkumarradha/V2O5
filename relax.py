from gpaw import restart
from ase.parallel import paropen as open
from ase.io import read
from gpaw import GPAW, PW, FermiDirac, restart
from ase.optimize import BFGS
from ase.parallel import world, parprint
from ase.io import read, write, Trajectory

atom = read('liv2o5.cif')
calc = GPAW(
    mode='lcao',
    basis='dzp',
    xc='PBE',
    kpts=(8, 8, 2),
    nbands=-10,
    random=True,  # random guess (needed if many empty bands required)
    occupations=FermiDirac(0.01),
    txt='v2o5_gs.txt')
atom.calc = calc
dyn = BFGS(atom, trajectory='v2o5.traj', logfile='log', restart='v2o5.pckl')
dyn.run(fmax=0.05)
Trajectory('v2o5.traj')[-1].write("relaxed_liv2o5.cif")