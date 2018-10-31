
function find_MOF(MOFname::String)
    if isfile(joinpath(PATH_TO_DATA, "crystals", MOFname * "_clean_min_charges.cif"))
        framework = Framework(MOFname * "_clean_min_charges.cif")
    elseif isfile(joinpath(PATH_TO_DATA, "crystals", MOFname * ".cif"))
        framework = Framework(MOFname * ".cif")
    else
        return "No MOF of that name exists in Crystal folder"
    return framework
    end
end

function read_pdb(MOFname::String)
    file = open(joinpath(PATH_TO_DATA, MOFname * ".pdb"))
    lines = readlines(file)
    close(file)
    return lines
end

function identify_atoms_in_pdb(lines::Array{String, 1}, atoms::DataFrames.DataFrame, atoms_in_MOF::Int64)
    count = 1
    #creates an array of string to hold atom identifiers
    id = Array{String, 1}(undef, atoms_in_MOF)

    #gets array of atomic names from PDB file
    for (i, line) in enumerate(lines)
        line = split(line)

        if line[1] == "HETATM"
            id[count] = line[end]
            count += 1
        end
    end

    return id
end

function feature_array(MOFname::String)
    #gets info from pdb file
    lines = read_pdb(MOFname)

    #gets atom information to find name and corresponding atomic number
    atoms = CSV.read(joinpath(PATH_TO_DATA, "atom_properties.csv"))

    #creates framework of the MOF in question
    MOF = find_MOF(MOFname)

    atoms_in_MOF = MOF.atoms.n_atoms

    #creates a nice array where each row corresponds to the number of the atom
    #in the MOF and the entry is the name of the atom
    atom_ids = identify_atoms_in_pdb(lines, atoms, atoms_in_MOF)

    #initializes feat_array
    feat_array = zeros(Float64, atoms_in_MOF, 2 * length(atoms[:atom]))


    #NOTE haven't touched anything below here yet

#=

    #modifies feature vector matrix with atom identification for each atom id
    for i = 1:framework.atoms.n_atoms
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

        #changed from element wise (for loop) to new nearest image
        nearest_image!(distance)

        #can be deleted later
        #for l = 1:framework.n_atoms
        #    nearest_image!(distance[:, l])
        #end

        #converts to cartesian distance
        distance = framework.box.f_to_c * distance

        for atom_2_id = 1:framework.atoms.n_atoms

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
    =#
    return atom_ids
end
