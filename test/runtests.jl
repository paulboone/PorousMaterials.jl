#!/usr/bin/env julia

# Details from http://www.stochasticlifestyle.com/finalizing-julia-package-documentation-testing-coverage-publishing/
# Start Test Script
using PorousMaterials
using Base.Test

# Run Tests

@printf("------------------------------\nTesting Crystal.jl\n\n")
frame = read_crystal_structure_file("test_structure.cif")
strip_numbers_from_atom_labels!(frame)
@testset "Crystal Tests" begin
	@test frame.box.f_to_c * frame.box.c_to_f ≈ eye(3)
	@test frame.box.a == frame.box.b == frame.box.c == 1.0
	@test frame.box.α == frame.box.β == frame.box.γ == 90 * (pi / 180)
	@test frame.box.Ω ≈ 1
	@test frame.xf == [0.0, 0.0, 0.0][:, :]
	@test frame.atoms == ["Zn"]
    @test frame.charges[1] ≈ 1.082578
    @test chemical_formula(frame) == Dict("Zn" => 1)
end;

@printf("------------------------------\nTesting Forcefield.jl\n\n")
const ljforcefield = read_forcefield_file("test_forcefield.csv", cutoffradius=12.5, mixing_rules="Lorentz-Berthelot")
rep_factors = replication_factors(frame.box, ljforcefield)
@testset "Forcefield Tests" begin
	@test ljforcefield.pure_sigmas["He"] == 1.0
	@test ljforcefield.pure_epsilons["Zn"] == 12.0
	@test ljforcefield.sigmas_squared["Zn"]["He"] == ((1.0 + 3.0) / 2) ^ 2
	@test ljforcefield.epsilons["He"]["Zn"] == sqrt(12.0 * 3.0)
	@test ljforcefield.cutoffradius_squared == 12.5 ^ 2
	@test rep_factors == (25, 25, 25)
end;

@printf("------------------------------\nTesting Energetics.jl\n\n")
molecule1 = Molecule(1, ["He"], [0.5, 0.5, 0.5][:,:], [0.0])
molecule2 = Molecule(1, ["He"], [0.5 + rep_factors[1], 0.5 + rep_factors[2], 0.5 + rep_factors[3]][:,:], [0.0])
frame2 = read_crystal_structure_file("SBMOF-1.cif")
rep_factors2 = replication_factors(frame2.box, ljforcefield)
molecule3 = Molecule(1,["Xe"], zeros(3,1), [0.0])
energy1 = vdw_energy(frame2, molecule3, ljforcefield, rep_factors2)
molecule3.x[1] = 0.494265; molecule3.x[2] = 2.22668; molecule3.x[3] = 0.450354;
energy2 = vdw_energy(frame2, molecule3, ljforcefield, rep_factors2)
@testset "Energetics Tests" begin
	@test vdw_energy(frame, molecule1, ljforcefield, rep_factors) ≈ vdw_energy(frame, molecule2, ljforcefield, rep_factors)
	@test vdw_energy(frame, molecule1, ljforcefield, (1,1,1)) ≈ 4 * ljforcefield.epsilons["He"]["Zn"] * ((ljforcefield.sigmas_squared["Zn"]["He"] / 0.75) ^ 6 - (ljforcefield.sigmas_squared["Zn"]["He"] / 0.75) ^ 3 )
	@test isapprox(energy1, -5041.58, atol = 0.005) # comparing to RASPA
	@test isapprox(energy2, 12945.838, atol = 0.005)
end;

@printf("------------------------------\n")

#
#GCMC Tests
#

"""

    use_nearest_image = completely_outside_box(molecule, sim_box)

returns true if the given molecule is completely outside of a given box
"""
function completely_outside_box(molecule::Molecule, box::Box)
    xf = box.c_to_f * molecule.x
    for coords = 1:3
        #if none of the coords are less than 1 it must be outside of the box
        if sum(xf[coords, :] .<= 1.0) == 0
            return true
        #if none of the coords are greater than 0 it must be outside of the box
        elseif sum(xf[coords, :] .>= 0.0) == 0
            return true
        end
    end
    return false
end

#
#INSERTION TESTS
#
molecules = Array{Molecule}(0)
repfactors = replication_factors(frame.box, ljforcefield)
sim_box = replicate_box(frame.box, repfactors)

for i = 1:100
    insert_molecule!(molecules, sim_box, "C")
    @assert(!completely_outside_box(molecules[i], sim_box), "Molecule outside of simulation box")
    @assert(length(molecules) == i, "Molecules array not modified")
    if i > 1
        @assert(sum(molecules[i - 1].x .≈ molecules[i].x) == 0, "Molecules not being inserted at random coordinates")
    end
end

#
#DELETION TESTS
#
for i = 1:100
    delete_molecule!(rand(1:length(molecules)), molecules)
    @assert(length(molecules) == 100 - i, "Molecule not deleted")
end

#
#TRANSLATION TESTS
#
molecules = [Molecule(1, ["C"], reshape(sim_box.f_to_c * [0.99, 0.99, 0.01], 3, :), [0.0])]
for i = 1:10000
    old_coords = molecules[1].x
    translate_molecule!(1, molecules, sim_box)
    translate_molecule!(1, molecules, sim_box)
    @assert(sum(molecules[1].x .== old_coords) == 0, "Molecule not moved")
    @assert(!completely_outside_box(molecules[1], sim_box),
        "Molecule completely outside of simulation box")
end

#
#VAN DER WAALS GUEST-GUEST ENERGY TESTS
#
sim_box = construct_box(25.0, 25.0, 25.0, π/2, π/2, π/2)
#TODO add CH4 to UFF file to have sigmas and epsilons ready
molecules = [Molecule(1, ["C"], reshape([5.0, 12.0, 12.0], 3, :), [0.0]),
    Molecule(1, ["C"], reshape([11.0, 12.0, 12.0], 3, :), [0.0])]

#calculate the energy between the two molecules using lennard_jones
dx1 = molecules[1].x .- molecules[2].x
energy_1 = lennard_jones(dot(dx1, dx1),
    ljforcefield.sigmas_squared[molecules[1].atoms[1]][molecules[2].atoms[1]],
    ljforcefield.epsilons[molecules[1].atoms[1]][molecules[2].atoms[1]])
#    get(ljforcefield.sigmas_squared[molecules[1].atoms[1]], molecules[2].atoms[1], -1),
#    get(ljforcefield.epsilons[molecules[1].atoms[1]], molecules[2].atoms[1], -1))

#making sure the assertion does what it is suppsoed to
@printf("\nEnergy calculated using lennard jones: %f \nEnergy calculated using gg vdw: %f\n", energy_1,
    guest_guest_vdw_energy(1, molecules, ljforcefield, sim_box))

#calculate the energy between two molecules to make sure it works correctly
@assert(guest_guest_vdw_energy(1, molecules, ljforcefield, sim_box) ≈ energy_1,
    "Did not calculate energy correctly between two molecules")

#calculate the energy between the same two molecules with the order reversed
#   to make sure that the energy between them stays the same
@assert(guest_guest_vdw_energy(1, molecules, ljforcefield, sim_box) ≈
    guest_guest_vdw_energy(2, molecules, ljforcefield, sim_box),
    "Energy is not the same between two molecules when reversing ")

#this removes the second molecule and replaces it with another
delete_molecule!(2, molecules)
push!(molecules, Molecule(1, ["C"], reshape([24.0, 12.0, 12.0], 3, :), [0.0]))
#calculate the energy between the two molecules using lennard_jones
dx2 = molecules[1].x .- molecules[2].x
energy_2 = lennard_jones(dot(dx2, dx2),
    ljforcefield.sigmas_squared[molecules[1].atoms[1]][molecules[2].atoms[1]],
    ljforcefield.epsilons[molecules[1].atoms[1]][molecules[2].atoms[1]])

#these molecules are outside of the cut off radius in the unit cell, but
#   they should be able to interact using nearest image convention
@assert(guest_guest_vdw_energy(1, molecules, ljforcefield, sim_box) ≈ energy_1,
    "Did not use nearest image convention")
push!(molecules, Molecule(1, ["C"], reshape([11.0, 12.0, 12.0], 3, :), [0.0]))

#making sure that the molecules I want are in the array
println(molecules)

#making sure the assertion does what it is suppsoed to
@printf("\nEnergy calculated using lennard jones: %f \nEnergy calculated using gg vdw: %f\n", 2 * energy_1,
    guest_guest_vdw_energy(1, molecules, ljforcefield, sim_box))

#inserts the original molecule and tests energy again. should be the sum of
#   the original two calculations
@assert(guest_guest_vdw_energy(1, molecules, ljforcefield, sim_box) ≈
    2 * energy_1,
    "Did not calculate guest-guest energy correctly for more than two molecules")