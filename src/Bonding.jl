using CSV
using DataFrames


#create dictionary for bonding rules

#reminder of syntax, should be deleted eventually
#framework = read_crystal_structure_file("SBMOF-1_cory.cif")


#function find_bonds(framework::Framework, bonding_rules)
function find_bonds(framework::Framework)

    #read in  atom_properties
    #will be useful for both the total number of n_atoms
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
        atom_1 = i
        for j = 1:framework.n_atoms
            #should ignore the bonds betweeen the same atom
            #distance would be 0 and could confuse bonding rules

            #defines atom_2 to be j
            atom_2 = j

            #find distance between atom we care about and current atom
            #dist = distance(atom_1, atom_2)

            #find characteristic bond length
            #charac_bond_length = bonding_rules(atom_1, atom_2)

            bond_length = 1
            charac_bond_length = 2
            #creates bond in feature array
            if bond_length < charac_bond_length
                #finds index of atom that is bonded 
                k = find(atoms_list .== string(framework.atoms[j]))
                feat_array[atom_1, n + k] += 1
            end
        end
    end
    return feat_array
end




function distance(atom_1, atom_2)
    cart = framework.box.f_to_c *framework.xf

end

function bonding_rules()

end
