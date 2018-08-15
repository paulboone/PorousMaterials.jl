using CSV
using DataFrames

#Easy copy and paste for framework, should be deleted eventually
#framework = read_crystal_structure_file("SBMOF-1_cory.cif")

#read in  atom_properties
atoms = CSV.read(PATH_TO_DATA * "atom_properties.csv")

#Creates the dictionary of covalent radii for bonding rules
function create_bonding_rules()
    bonding_rules = Dict{Symbol,Float64}()
    for i = 1:length(atoms[:atom])
        if ! ismissing(atoms[Symbol("singlebondcovalentradius[Angstrom]")][i])
            bonding_rules[Symbol(atoms[:atom][i])] = Float64(atoms[Symbol("singlebondcovalentradius[Angstrom]")][i])
        end
    end
    return bonding_rules
end

#function find_bonds(framework::Framework, bonding_rules)
function find_bonds(framework::Framework)

    bonding_rules = create_bonding_rules()

    n = length(atoms[:atom])

    #initializes feat_array
    feat_array = zeros(Float64, framework.n_atoms, 2 * n)

    #modifies feature vector matrix with atom identification for each atom id
    for i = 1:framework.n_atoms
        for j = 1:n
            #checks if name of atom is the same as the corresponding posiiton in the feature matrix
            if string(framework.atoms[i]) == atoms[:atom][j]
                feat_array[i, j] = 1
                break
            end
        end
    end

    #actually finds the bonds
    for atom_1_id = 1:framework.n_atoms

        #changes framework.xf to be a vector of distance from the atom in question
        #to all other atoms in the framework
        atom_1_vector = framework.xf[:, atom_1_id]
        distance = framework.xf .- atom_1_vector

        #couldn't get the element wise nearest_image to work
        #so the fix is to use a for loop
        for l = 1:framework.n_atoms
            #NEEDS TO BE RECHECKED TO MAKE SURE IT IS ACTUALLY DOING THE NEAREST IMAGE
            nearest_image!(distance[:, l])
        end

        #converts to cartesian distance
        distance = framework.box.f_to_c * distance

        for atom_2_id = 1:framework.n_atoms

            #find magnitude of vector between atom we care about and current atom
            bond_length = norm(distance[:, atom_2_id])

            #find characteristic bond length
            charac_bond_length = bonding_rules[framework.atoms[atom_1_id]] + bonding_rules[framework.atoms[atom_2_id]]

            #this is the percentage representing the furthest and shortest distance the bond
            #can be from the average covalent radii sum and still be condidered bonded
            tol = 0.1

            #creates bond in feature array
            if ((charac_bond_length * (1 - tol)) < bond_length < charac_bond_length * (1 + tol))
                #finds index of atom that is bonded
                k = find(atoms[:atom] .== string(framework.atoms[atom_2_id]))
                #add bond to feature array
                feat_array[atom_1_id, n + k] += 1

            end
        end
    end
    return feat_array
end
