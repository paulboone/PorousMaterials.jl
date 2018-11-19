using PorousMaterials

framework = Framework("SBMOF-1.cif")
molecule = Molecule("Xe")
forcefield = LJForceField("UFF.csv")
grid = energy_grid(framework, molecule, forcefield, n_pts=(10, 10, 10))
bgrid = binarize(grid, 0.0)
