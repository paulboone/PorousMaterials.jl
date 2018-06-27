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
 #  TODO delete this
 # function old_read_ionization_data()
 #     # Data is read in as eV. This is changed to K by dividing by the Boltzmann constant, `kb`
 #     f = open("data/ionizationData.csv")
 #     lines = readlines(f)
 #     close(f)
 # 
 #     ionization_potential = Dict{Symbol, Array{Float64, 1}}()
 #     n_ionization = Dict{Symbol, Int64}()
 #     status = Dict{Symbol, AbstractString}()
 #     chargeCenter = Dict{Symbol, Int64}()
 # 
 #     for line in lines[2:end]
 #         line = split(line, ",")
 #         atom = Symbol(line[1])
 #         status[atom] = line[2]
 # 
 #         chargeCenter[atom] = parse(Int64, line[3])
 # 
 #         if line[4] == "<0.5"
 #             electron_affinity = 0.5 / kb
 #         else
 #             electron_affinity = parse(Float64, line[4]) / kb
 #         end
 #         if atom == :H
 #             electron_affinity = -2 / kb
 #         end
 # 
 #         io_potential = tryparse(Float64, line[5])
 #         potential = Array{Float64, 1}()
 #         push!(potential, electron_affinity)
 #         n_io = 0
 #         while !isnull(io_potential)
 #             n_io += 1
 #             push!(potential, get(io_potential) / kb) 
 #             if n_io + 5 > length(line)
 #                 break
 #             end
 #             io_potential = tryparse(Float64, line[5 + n_io])
 #         end
 # 
 #         ionization_potential[atom] = potential
 #         n_ionization[atom] = n_io
 #     end
 # 
 #     return IonizationData(ionization_potential, n_ionization, status, chargeCenter)
 # end

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
function ϕ_ij(framework::Framework, i::Int, j::Int, ion_data::IonizationData; sr_cutoff::Float64=12.5)
    rep_factors = replication_factors(framework, sr_cutoff_r)
    eparams, kvecs, eikar, eikbr, eikcr = setup_Ewald_sum(sr_cutoff_r, framework.box)
    
    ϕ = 0.0
    # build molecules out of framework atoms in question.
    charge_i = PointCharge(1.0, framework.box.f_to_c * framework.xf[:, i])
    
    # short-range contribution
    for ra = 0:(repfactors[1] - 1), rb = 0:(repfactors[2] - 1), rc = 0:(repfactors[3] - 1)
        # disallow self-interaction
        if (i == j) && ([ra, rb, rc] == [0, 0, 0])
            continue
        end
        # fractional coordinate of all atom j's within short-range cutoff
        xf = framework.box.f_to_c * (framework.xf[:, j] + 1.0 * [ra, rb, rc])
        charge_j = PointCharge(1.0, xf)
        # TODO divide by two or not?
        if i == j
            ϕ += _ϕ_sr(charge_i, charge_j, eparams) / FOUR_PI_ϵ₀
        else
            ϕ += _ϕ_sr(charge_i, charge_j, eparams) / FOUR_PI_ϵ₀ / 2.0
        end
    end

    # long-range contribution
    charge_j = PointCharge(1.0, framework.box.f_to_c * framework.xf[:, j])
    # TODO divide by two or not?
    if i == j
        ϕ += _ϕ_lr(charge_i, charge_j, eparams, kvecs, eikar, eikbr, eikcr)
    else
        ϕ += _ϕ_lr(charge_i, charge_j, eparams, kvecs, eikar, eikbr, eikcr) / 2.0
    end

    # spurious self-interaction contribution
    if i == j
        ϕ -= 2.0 * eparams.α / sqrt(π) / FOUR_PI_ϵ₀
    end

    # ionization energy   
    if i == j
        ϕ += ∂²E_i∂_Q_i²(ion_data, framework.atoms[i])
    end
    
    # orbital overlap energy
    # TODO check Chris Wilmer's code: is this chemical hardness or double deriviave
    if i != j
        # geometric mean of chemical hardni divided by 4 π ϵ₀
        J_mean_x_four_pi_ϵ = sqrt(∂²E_i∂_Q_i²(ion_data, framework.atoms[i]) * ∂²E_i∂_Q_i²(ion_data, framework.atoms[j])) * FOUR_PI_ϵ₀
        # distance between atom i and j
        dxf = framework.xf[:, i] - framework.xf[:, j]
        nearest_image!(dxf, (1, 1, 1))
        dx = framework.box.f_to_c * dxf
        r = norm(dx)

        ϕ += exp(- (J_mean_x_four_pi_ϵ * r) ^ 2) * (J_mean_x_four_pi_ϵ - J_mean_x_four_pi_ϵ ^ 2 * r - 1.0 / r) / FOUR_PI_ϵ₀ / 2.0
    end
    return ϕ

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
    A = ones(frame.n_atoms, frame.n_atoms)
    for i = 1:frame.n_atoms-1
        for j = 1:frame.n_atoms
            A[i, j] = ϕ
        end
    end
    return A
end
