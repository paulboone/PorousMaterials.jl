# cof_gcmc.jl
# Adrian Henle 2019
using PorousMaterials
using PyPlot
using Printf
using CSV
using DataFrames

# Simulation Parameters
crystals = ["COF102", "COF103"]   # List of Host Structures
xtal_extension = Dict("COF102" => ".cif", "COF103" => ".cssr")
xtal_densities = Dict("COF102" => 420.54, "COF103" => 388.72)
temp = 298.0                                # Temperature [K]

molecule = Molecule("CH4")
ljff = LJForceField("Dreiding.csv", cutoffradius=12.5, mixing_rules="Lorentz-Berthelot")

for crystal in crystals
    frame = Framework(crystal * xtal_extension[crystal])
    strip_numbers_from_atom_labels!(frame)
    @assert isapprox(crystal_density(frame), xtal_densities[crystal], atol=0.1)
    
    # Load results from published simulations
    # source: https://pubs.acs.org/doi/abs/10.1021/jp507152j
    lit_data = CSV.read(crystal * "_dreiding_298K.txt", delim = "\t", header = 4, datarow = 5);
    names!(lit_data, [:P, :N]); # [Pa] vs. [mol/kg]
    lit_data.P = lit_data.P ./ 100000; # Pressure now in [bar]
    # Molar Mass of methane
    MW = 16.04;
    lit_data.N = lit_data.N .* MW; # Adsorption now in [g / kg]

    results = adsorption_isotherm(
        frame, molecule, temp, lit_data[:, :P], ljff, eos=:PengRobinson,
        n_burn_cycles=10000, n_sample_cycles=10000
    )
    n_pred = [result["⟨N⟩ (mmol/g)"] for result in results]
    
    # plot result
    figure()
    scatter(lit_data.P, n_pred)
    scatter(lit_data.P, lit_data.N)
    title(@sprintf("Isotherms for Methane in %s at %0.1f K", crystal, temp))
    xlabel("Pressure (bar)")
    ylabel("⟨N⟩ (mg Methane / g COF)")
    legend(["PorousMaterials.jl", "J. Phys. Chem. C"])
    savefig("$(crystal)_ch4_lit_vs_porousmaterials.png", format="png")
end
