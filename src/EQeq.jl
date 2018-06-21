using PorousMaterials
using OffsetArrays
using CSV
using DataFrames

# todo delete these constants
const ϵ₀ = 4.7622424954949676e-7  # \epsilon_0 vacuum permittivity units: electron charge^2 /(A - K)
const ϵᵣ = 1.67 # Global dielectric strength to prevent infinite charge separation
#const K = 14.4
const R = 8.314 # J/ (mol K)
const kb = 8.6173303e-5 # eV / K
"""

"""
struct IonizationData
    """
    Array of ionization potentials for the atoms (units: K).
        entry i is ionization potential for atom i
    """
    ionization_potentials::Dict{Symbol, Array{Float64, 1}}
    "Electron affinity (units: K)"
    electron_affinity::Dict{Symbol, Float64}
    """
    Taylor expansion of energy of atom i is:
        E_i(Q_i) = E_i(Q_i*) + ∂E_i (Q_i - Q_i*) ...
        this holds Q_i^*, center about which Taylor expansion is performed.
    """
    charge_center::Dict{Symbol, Int64}
end

function read_ionization_data()
    df = CSV.read("data/ionizationData.csv")
    
    ionization_potentials = Dict{Symbol, Array{Float64, 1}}()
    electron_affinity = Dict{Symbol, Float64}()
    charge_center = Dict{Symbol, Int64}()
    
    # TODO put units in column names of data.
    for row in eachrow(df)
        atom = Symbol(row[:atom])
        charge_center[atom] = row[:chargeCenter]
        electron_affinity[atom] = row[:Electron_affinity]
        ionization_potentials[atom] = Float64[]
        for i = 1:8
            ie = row[Symbol("Ionization_energy_$i")]
            if (ie == "na") | (ie == "np") # not possible/ not available?
                break
            end
            push!(ionization_potentials[atom], parse(Float64, row[ie]))
        end
    end
    return IonizationData(ionization_potentials, electron_affinity, charge_centers)
end

function old_read_ionization_data()
    # Data is read in as eV. This is changed to K by dividing by the Boltzmann constant, `kb`
    f = open("data/ionizationData.csv")
    lines = readlines(f)
    close(f)

    ionization_potential = Dict{Symbol, Array{Float64, 1}}()
    n_ionization = Dict{Symbol, Int64}()
    status = Dict{Symbol, AbstractString}()
    chargeCenter = Dict{Symbol, Int64}()

    for line in lines[2:end]
        line = split(line, ",")
        atom = Symbol(line[1])
        status[atom] = line[2]

        chargeCenter[atom] = parse(Int64, line[3])

        if line[4] == "<0.5"
            electron_affinity = 0.5 / kb
        else
            electron_affinity = parse(Float64, line[4]) / kb
        end
        if atom == :H
            electron_affinity = -2 / kb
        end

        io_potential = tryparse(Float64, line[5])
        potential = Array{Float64, 1}()
        push!(potential, electron_affinity)
        n_io = 0
        while !isnull(io_potential)
            n_io += 1
            push!(potential, get(io_potential) / kb) 
            if n_io + 5 > length(line)
                break
            end
            io_potential = tryparse(Float64, line[5 + n_io])
        end

        ionization_potential[atom] = potential
        n_ionization[atom] = n_io
    end

    return IonizationData(ionization_potential, n_ionization, status, chargeCenter)
end

"""
chemical hardness J:= = I₁ - A = ∂²E_i/∂Q_i², 
where I₁ is ionization potential and A is electron affinity
"""
function chemical_hardness(ion_data::IonizationData, atom::Symbol)
    return ion_data.ionization_potentials[atom][1] - ion_data.electron_affinity[atom]
end

"""
electronegativity Χ := (I₁ + A) / 2
where I₁ is ionization potential and A is electron affinity
"""
function electronegativity(ion_data::IonizationData, atom::Symbol)
    return (ion_data.ionization_potentials[atom][1] + ion_data.electron_affinity[atom]) / 2.0
end

function ∂E_i∂_Q_i(ion_data::IonizationData, atom::Symbol)
    charge_center = ion_data.charge_center[atom]
    if charge_center == 0
        return electronegativity(ion_data, atom)
    else
        return (ion_data.ionization_potentials[charge_center + 1] + ion_data.ionization_potentials[charge_center]) / 2.0
    end
end

function ∂²E_i∂_Q_i²(ion_data::IonizationData, atom::Symbol)
    charge_center = ion_data.charge_center[atom]
    if charge_center == 0
        return chemical_hardness(ion_data, atom)
    else
        return ion_data.ionization_potentials[charge_center + 1] - ion_data.ionization_potentials[charge_center]
    end
end

"""
Entries in matrix
"""
function compute_ϕ_ij(framework::Framework, i::Int, j::Int, ion_data::IonizationData; sr_cutoff::Float64=12.5)
    rep_factors = replication_factors(framework, sr_cutoff_r)
    eparams, kvecs, eikar, eikbr, eikcr = setup_Ewald_sum(sr_cutoff_r, framework.box)

    # build molecules out of framework atoms in question.
    x_i = framework.box.f_to_c * framework.xf[:, i]
    x_j = framework.box.f_to_c * framework.xf[:, j]
    atom_i = Molecule(:atom_i, 
                      LennardJonesSphere[],
                      PointCharge(1.0, x_i)
                      x_i)
    atom_j = Molecule(:atom_j, 
                      LennardJonesSphere[],
                      PointCharge(1.0, x_j)
                      x_j)

    ϕ = PorousMaterials.electrostatic_potential_energy([atom_i, atom_j], eparams, kvecs, eikar, eikbr, eikcr)

    return ϕ.sr + ϕ.lr + ϕ.self * 2.0 + 

 #     K = 1 / (4 * π * ϵ₀ * ϵᵣ)
 #     J = compute_pure_J(ioData)
 # 
 #     α = compute_α(frame, repfactors, eparams, ind)
 #     β = compute_β(frame, repfactors, eparams, eikar, eikbr, eikcr, ind, kvecs)
 #     if ind[1] == ind[2]
 #         return K * α - eparams.α * K / sqrt(π) + β + J[frame.atoms[ind[1]]]
 #     else
 #         E0 = compute_E0(frame, repfactors, J, ioData, ind)
 #         return K / 2 * (α + E0) + β / (4 * π)
 #     end
end

function assign_charges(frame::Framework, ioData::IonizationData, repfactors::Tuple{Int, Int, Int}, sr_cutoff_r::Float64, ϵ::Float64)
    A_J = ones(frame.n_atoms, frame.n_atoms)
    for i = 1:frame.n_atoms-1
        for j = 1:frame.n_atoms
            A_J[i,j] = compute_JPE(frame, ioData, repfactors, sr_cutoff_r, ϵ, (i, j)) - compute_JPE(frame, ioData, repfactors, sr_cutoff_r, ϵ, (i+1, j))
        end
    end
    return A_J
end

function compute_α(frame::Framework, repfactors::Tuple{Int, Int, Int}, eparams::EwaldParams, ind::Tuple{Int, Int})
    α = 0.0
    for ra = 0:(repfactors[1] - 1), rb = 0:(repfactors[2] - 1), rc = 0:(repfactors[3] - 1)
        dxf = frame.xf[:,ind[1]] - (frame.xf[:,ind[2]] + [ra, rb, rc])
        if norm(dxf) == 0
            continue
        end
        nearest_image!(dxf, repfactors)
        dx = frame.box.f_to_c * dxf
        @inbounds @fastmath r = sqrt(dx[1,1] * dx[1,1] + dx[2,1] * dx[2,1] + dx[3,1] * dx[3,1])
        if r < eparams.sr_cutoff_r
            @inbounds @fastmath α += erfc(r * eparams.α) / r
        end
    end
    return α
    # Returns units of Å^(-1)
end

function compute_β(frame::Framework, repfactors::Tuple{Int, Int, Int}, eparams::EwaldParams, eikar::OffsetArray{Complex{Float64}}, eikbr::OffsetArray{Complex{Float64}}, eikcr::OffsetArray{Complex{Float64}}, ind::Tuple{Int, Int}, kvectors::Array{Kvector, 1})
    β = 0.0
    for ra = 0:(repfactors[1] - 1), rb = 0:(repfactors[2] - 1), rc = 0:(repfactors[3] - 1)
        @inbounds dxf = frame.xf[:,ind[1]] - (frame.xf[:,ind[2]] + [ra, rb, rc])
        if norm(dxf) == 0
            continue
        end
        @inbounds dx = frame.box.f_to_c * dxf
        @inbounds k_dot_dx = transpose(eparams.box.reciprocal_lattice) * dx

        fill_eikr!(eikar, k_dot_dx[1, 1], eparams.kreps[1], false)
        fill_eikr!(eikbr, k_dot_dx[2, 1], eparams.kreps[2], true)
        fill_eikr!(eikcr, k_dot_dx[3, 1], eparams.kreps[3], true)

        for kvector in kvectors
            @unsafe @inbounds β += kvector.wt * real(eikar[kvector.ka] * eikbr[kvector.kb] * eikcr[kvector.kc]) # [β] ≡ Å^3 * K * e^(-2)
        end
    end
    return β / eparams.box.Ω # This has already been divided by ϵ₀, so we'll have to be careful when multiplying by K (according to the paper)
    # Returns units of K * e^(-2)
end

function compute_E0(frame::Framework, repfactors::Tuple{Int, Int, Int}, pure_J::Dict{Symbol, Float64}, ioData::IonizationData, ind::Tuple{Int, Int})
    E0 = 0.0
    K = 1 / (4 * π * ϵ₀ * ϵᵣ)
    for ra = 0:(repfactors[1] - 1), rb = 0:(repfactors[2] - 1), rc = 0:(repfactors[3] - 1)
        @inbounds dxf = frame.xf[:, ind[1]] - (frame.xf[:, ind[2]] + [ra, rb, rc])
        if norm(dxf) == 0
            continue
        end
        nearest_image!(dxf, repfactors)
        @inbounds dx = frame.box.f_to_c * dxf
        r = norm(dx)
        J = (pure_J[frame.atoms[ind[1]]] + pure_J[frame.atoms[ind[2]]]) / 2
        E0 += exp(-(J * r / K)^2) * (J/K - J^2*r/K^2 - 1/r) # [E0] ≡ Å^(-1)
    end
    return E0
    #Returns units of Å^(-1)
end


function compute_b_vector(frame::Framework, pure_J::Dict{Symbol, Float64}, ioData::IonizationData)
    b = zeros(frame.n_atoms)
    for i = 1:length(b)-1
        atom1 = frame.atoms[i]
        atom2 = frame.atoms[i+1]
        #ionization_potential[:Atom][1] is Electron Affinity
        #ionization_potential[:Atom][2...] are 1st, 2nd, ... ionization energies
        b[i] = ioData.ionization_potential[atom2][1] - ioData.ionization_potential[atom1][1] #- pure_J[atom2] * ioData.chargeCenter[atom2] + pure_J[atom1] * ioData.chargeCenter[atom1]
    end
    return b
end
