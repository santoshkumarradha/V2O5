from ase import Atoms
from gpaw import GPAW, FermiDirac, Davidson, Mixer, restart
from gpaw.cdft.cdft import CDFT
from gpaw.cdft.cdft_coupling import CouplingParameters
from gpaw import restart
from ase.parallel import world, parprint
from ase.io import read
import numpy as np
start_fresh = False
if start_fresh:
    atom = read('relaxed_liv2o5.cif')
    mom = np.zeros(atom.get_global_number_of_atoms())
    mom[np.where(atom.get_atomic_numbers() == 23)] = 1
    atom.set_initial_magnetic_moments(magmoms=mom)
    calc = GPAW(
        # mode='fd',
        basis='dzp',
        xc='PBE',
        kpts=(8, 8, 2),
        nbands=-30,
        random=True,  # random guess (needed if many empty bands required)
        occupations=FermiDirac(0.01),
        txt='v2o5_gs.txt',
        symmetry='off',
        spinpol=True,
    )

    do_cdft = False
    if do_cdft:
        c = 1
        cdft = CDFT(
            calc=calc,
            atoms=atom,
            charge_regions=[[0]],  # choose atom 0 as the constrained region
            charges=[c],  # constrain +1 charge
            charge_coefs=[2.7],  # initial guess for Vc
            method='L-BFGS-B',  # Vc optimization method
            txt='v2o5_c={}.cdft'.format(c),  # cDFT output file
            minimizer_options={'gtol': 0.01})  # tolerance for cdft
        atom.calc = cdft

    else:
        atom.calc = calc

    energy = atom.get_potential_energy()
    calc.write('v2o5_gs.gpw', mode="all")

else:
    atom, calc = restart('v2o5_gs.gpw')

if False:
    calc = GPAW('v2o5_gs.gpw',
                fixdensity=True,
                symmetry='off',
                kpts={
                    'path': 'GXSYG',
                    'npoints': 80
                },
                txt='v2o5_bands.txt',
                convergence={'bands': 'CBM+2.5'})
    atom.calc = calc
    calc.get_potential_energy()
    calc.write('v2o5_bands.gpw', mode="all")

atom, calc = restart('v2o5_bands.gpw')
bs = calc.band_structure()
bs.plot(filename='bandstructure.png', show=False, emax=10.0, emin=-10.0)