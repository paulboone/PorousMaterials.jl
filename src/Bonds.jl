"""
    bonding_rule = BondingRule(:Ca, :O, 0.4, 2.0)
    bonding_rules = [BondingRule(:H, :*, 0.4, 1.2),
                     BondingRule(:*, :*, 0.4, 1.9)]

A rule for determining if two atoms within a framework are bonded. 

# Attributes
-`species_i::Symbol`: One of the atoms types for this bond rule
-`species_j::Symbol`: The other atom type for this bond rule
-`min_dist`: The minimum distance between the atoms for bonding to occur
-`max_dist`: The maximum distance between the atoms for bonding to occur
"""
struct BondingRule
    species_i::Symbol
    species_j::Symbol
    min_dist::Float64
    max_dist::Float64
end

"""
    atom_to_radius = cordero_covalent_atomic_radii()

Read in a list of covalent atomic radii as proposed by Cordero et al. DOI 10.1039/B801115J.

Caveats:
* Csp, Csp2, Csp3 are specified in Table 2 of Cordero et al. and included here. C was not specified by Cordero et al. We assigned radius to C as the maximum radius of (Csp, Csp2, Csp3).
* Mn, Fe, and Co were given separate radii for high spin and low spin in Table 2 of Cordero et al. We assigned radii to Mn, Fe, and Co based on the maximum of the radii in the low and high spin states.

# Returns
* `atom_to_radius::Dict{Symbol, Float64}`: dictionary that maps atomic symbol to covalent radius. `atom_to_radius[:H]` gives covalent atomic radius of H.
"""
function cordero_covalent_atomic_radii()
    wherez_radii_file = normpath(joinpath(pathof(PorousMaterials), "..", "covalent_radii.csv"))
    df = CSV.read(wherez_radii_file, comment="#")
    atom_to_radius = Dict{Symbol, Float64}()
    for atom in eachrow(df)
        atom_to_radius[Symbol(atom[:atom])] = atom[:covalent_radius_A]
    end
    return atom_to_radius
end

"""
    default_bondingrules = default_bondingrules()

Returns the default bonding rules. Using `append!` and/or `prepend!` to add to the default bonding rules:

# Example
```
bond_rules = default_bondingrules()
prepend!(bond_rules, BondingRule(:Cu, :*, 0.1, 2.6))
```

# Returns
-`default_bondingrules::Array{BondingRule, 1}`: The default bonding rules: `[BondingRule(:*, :*, 0.4, 1.2), BondingRule(:*, :*, 0.4, 1.9)]`
"""
default_bondingrules() = [BondingRule(:H, :*, 0.4, 1.2), BondingRule(:*, :*, 0.4, 1.9)]

"""
    remove_bonds!(framework)

Remove all bonds from a framework structure.

# Arguments
-`framework::Framework`: the framework that bonds wil be removed from
"""
function remove_bonds!(framework::Framework)
    while ne(framework.bonds) > 0
        rem_edge!(framework.bonds, collect(edges(framework.bonds))[1].src, collect(edges(framework.bonds))[1].dst)
    end
end

"""
    are_atoms_bonded = is_bonded(framework, i, j, bonding_rules=[BondingRule(:H, :*, 0.4, 1.2), BondingRule(:*, :*, 0.4, 1.9)],
                                 include_bonds_across_periodic_boundaries=true)

Checks to see if atoms `i` and `j` in `framework` are bonded according to the `bonding_rules`.

# Arguments
-`framework::Framework`: The framework that bonds will be added to
-`i::Int`: Index of the first atom
-`j::Int`: Index of the second atom
-`bonding_rules::Array{BondingRule, 1}`: The array of bonding rules that will
    be used to fill the bonding information. They are applied in the order that
    they appear.
-`include_bonds_across_periodic_boundaries::Bool`: Whether to check across the
    periodic boundary when calculating bonds

# Returns
-`are_atoms_bonded::Bool`: Whether atoms `i` and `j` are bonded according to `bonding_rules`

"""
function is_bonded(framework::Framework, i::Int64, j::Int64,
                   bonding_rules::Array{BondingRule, 1}=[BondingRule(:H, :*, 0.4, 1.2), BondingRule(:*, :*, 0.4, 1.9)];
                   include_bonds_across_periodic_boundaries::Bool=true)
    species_i = framework.atoms.species[i]
    species_j = framework.atoms.species[j]

    cartesian_dist_between_atoms = distance(framework, i, j, include_bonds_across_periodic_boundaries)

    # loop over possible bonding rules
    for br in bonding_rules
        # determine if the atom species correspond to the species in `bonding_rules`
        species_match = false
        if br.species_i == :* && br.species_j == :*
            species_match = true
        elseif br.species_i == :* && (species_i == br.species_j || species_j == br.species_j)
            species_match = true
        elseif br.species_j == :* && (species_i == br.species_i || species_j == br.species_j)
            species_match = true
        elseif (species_i == br.species_i && species_j == br.species_j) || (species_j == br.species_i && species_i == br.species_j)
            species_match = true
        end

        if species_match
            # determine if the atoms are close enough to bond
            if br.min_dist < cartesian_dist_between_atoms && br.max_dist > cartesian_dist_between_atoms
                return true
            end
        end
    end
    return false
end

"""
    infer_bonds!(framework, include_bonds_across_periodic_boundaries, 
                    bonding_rules=[BondingRule(:H, :*, 0.4, 1.2), BondingRule(:*, :*, 0.4, 1.9)])

Populate the bonds in the framework object based on the bonding rules. If a
pair doesn't have a suitable rule then they will not be considered bonded. 

`:*` is considered a wildcard and can be substituted for any species. It is a
good idea to include a bonding rule between two `:*` to allow any atoms to bond
as long as they are close enough.

The bonding rules are hierarchical, i.e. the first bonding rule takes precedence over the latter ones.

# Arguments
-`framework::Framework`: The framework that bonds will be added to
-`include_bonds_across_periodic_boundaries::Bool`: Whether to check across the
    periodic boundary when calculating bonds
-`bonding_rules::Array{BondingRule, 1}`: The array of bonding rules that will
    be used to fill the bonding information. They are applied in the order that
    they appear.
"""
function infer_bonds!(framework::Framework, include_bonds_across_periodic_boundaries::Bool,
                      bonding_rules::Array{BondingRule, 1}=
                      [BondingRule(:H, :*, 0.4, 1.2), BondingRule(:*, :*, 0.4, 1.9)])
    @assert ne(framework.bonds) == 0 @sprintf("The framework %s already has bonds. Remove them with the `remove_bonds!` function before inferring new ones.", framework.name)

    # loop over every atom
    for i in 1:framework.atoms.n_atoms
        # loop over every unique pair of atoms
        for j in i+1:framework.atoms.n_atoms
            if is_bonded(framework, i, j, bonding_rules; include_bonds_across_periodic_boundaries=include_bonds_across_periodic_boundaries)
                add_edge!(framework.bonds, i, j)
            end
        end
    end
end

"""
    sane_bonds = bond_sanity_check(framework)

Run sanity checks on `framework.bonds`.
* is the bond graph fully connected? i.e. does every vertex (=atom) in the bond graph have at least one edge?
* each hydrogen can have only one bond
* each carbon can have a maximum of four bonds

if sanity checks fail, refer to [`write_bond_information`](@ref) to write a .vtk to visualize the bonds.

Print warnings when sanity checks fail.
Return `true` if sanity checks pass, `false` otherwise.
"""
function bond_sanity_check(framework::Framework)
    sane_bonds = true
    for a = 1:framework.atoms.n_atoms
        ns = neighbors(framework.bonds, a)
        # is the graph fully connected?
        if length(ns) == 0
            @warn "atom $a = $(framework.atoms.species[a]) in $(framework.name) is not bonded to any other atom."
            sane_bonds = false
        end
        # does hydrogen have only one bond?
        if (framework.atoms.species[a] == :H) && (length(ns) > 1)
            @warn "hydrogen atom $a in $(framework.name) is bonded to more than one atom!"
        end
        # does carbon have greater than four bonds?
        if (framework.atoms.species[a] == :C) && (length(ns) > 4)
            @warn "carbon atom $a in $(framework.name) is bonded to more than four atoms!"
        end
    end
    return sane_bonds
end

"""
    bonds_equal = compare_bonds_in_framework(framework1, framework2, atol=0.0)

Returns whether the bonds defined in framework1 are the same as the bonds
defined in framework2. It checks whether the atoms in the same positions
have the same bonds.

# Arguments
-`framework1::Framework`: The first framework
-`framework2::Framework`: The second framework
-`atol::Float64`: absolute tolerance for the comparison of coordinates in the framework

# Returns
-`bonds_equal::Bool`: Wether the bonds in framework1 and framework2 are equal
"""
function compare_bonds_in_framework(fi::Framework, fj::Framework; atol::Float64=0.0)
    if ne(fi.bonds) != ne(fj.bonds)
        return false
    end

    num_in_common = 0
    for edge_i in collect(edges(fi.bonds))
        for edge_j in collect(edges(fj.bonds))
            # either the bond matches going src-src dst-dst
            if  (fi.atoms.species[edge_i.src] == fj.atoms.species[edge_j.src] &&
                 fi.atoms.species[edge_i.dst] == fj.atoms.species[edge_j.dst] &&
                 isapprox(fi.atoms.xf[:, edge_i.src], fj.atoms.xf[:, edge_j.src]; atol=atol) &&
                 isapprox(fi.atoms.xf[:, edge_i.dst], fj.atoms.xf[:, edge_j.dst]; atol=atol)) ||
                # or the bond matches going src-dst dst-src
                (fi.atoms.species[edge_i.src] == fj.atoms.species[edge_j.dst] &&
                 fi.atoms.species[edge_i.dst] == fj.atoms.species[edge_j.src] &&
                 isapprox(fi.atoms.xf[:, edge_i.src], fj.atoms.xf[:, edge_j.dst]; atol=atol) &&
                 isapprox(fi.atoms.xf[:, edge_i.dst], fj.atoms.xf[:, edge_j.src]; atol=atol))
                num_in_common += 1
                break
            end
        end
    end
    return num_in_common == ne(fi.bonds) && num_in_common == ne(fj.bonds)
end
