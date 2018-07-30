using CSV
using DataFrames


#create dictionary for bonding rules

#reminder of syntax, should be deleted eventually
#framework = read_crystal_structure_file("SBMOF-1_cory.cif")


#function find_bonds(framework::Framework, bonding_rules)
function find_bonds(framework::Framework)

    #read in  atom_properties
    #will be useful for both the total number of atoms
    #as well as for comparing atoms and being able to convert from atom_id to atom name (maybe)
    atoms = CSV.read(PATH_TO_DATA * "atom_properties.csv")
    n = length(atoms[:atom])

    #initializes feat_array
    feat_array = zeros(Float64, framework.n_atoms, 2 * n)

    #creates a matrix just for use in the second iterator
    #ALMOST CERTAINTLY A BETTER WAY TO DO THIS!!
    atoms_list = atoms[:atom]

    #modifies feature vector matrix with atom identification for each atom id
    for i = 1:framework.n_atoms
        for j = 1:n
            #checks if name of atom is the same as the corresponding posiiton in the feature matrix
            if string(framework.atoms[i]) == atoms_list[j]
                feat_array[i, j] = 1
                break
            end
        end
    end

    for i = 1:framework.n_atoms
        #defines atom_1 to be i
        atom_1_id = i

        for j = 1:framework.n_atoms
            #defines atom_2 to be j
            atom_2_id = j

            #find distance between atom we care about and current atom
            #dist = distance(atom_1, atom_2)

            #find characteristic bond length
            charac_bond_length = bonding_rules(atom_1_id, atom_2_id, framework)

            bond_length = 1
            #creates bond in feature array
            if bond_length < charac_bond_length && bond_length > 0.4
                #finds index of atom that is bonded
                k = find(atoms_list .== string(framework.atoms[atom_2_id]))
                #add bond to feature array
                feat_array[atom_1_id, n + k] += 1
            end
        end
    end
    return feat_array
end

function distance(atom_1::Int64, atom_2::Int64)
    cart = framework.box.f_to_c * framework.xf

end

function bonding_rules(atom_1_id::Int64, atom_2_id::Int64, framework::Framework)
    #finds which atoms are being bonded
    atom_1 = framework.atoms[atom_1_id]
    atom_2 = framework.atoms[atom_2_id]

    #bonding rules for Hydrogen
    if atom_1 == :H || atom_2 == :H
        charac_bond_length = 1.2

    #Bonding rules for Calcium
elseif atom_1 == :Ca || atom_2 == :Ca
        charac_bond_length = 2.5

    #Bonding rules for the rest
    else
        charac_bond_length = 1.9
    end

    return charac_bond_length
end
