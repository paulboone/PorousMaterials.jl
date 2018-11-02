
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

function extract_from_pdb(lines::Array{String, 1}, atoms::DataFrames.DataFrame, MOF::Framework)

    #NOTE right now the code can only handle 10 bonds, which seems to be enough to me

    count = 1
    count_2 = 1
    #creates an array of string to hold atom identifiers and coords from the pdb
    info = Array{String, 2}(undef, MOF.atoms.n_atoms, 14)

    #array to hold the cart coords once theyve been translated to floats from strings
    pdb_cart_coords = zeros(Float64, MOF.atoms.n_atoms, 3)

    #arrays used to store atomic number, and bonding info
    unordered_MOF_info = zeros(Float64, MOF.atoms.n_atoms, 11)
    ordered_MOF_info = zeros(Float64, MOF.atoms.n_atoms, 11)

    #stores array of atomic names from PDB file
    #NOTE This doens't appear to be working
    #       all entries appear as undef
    for (i, line) in enumerate(lines)
        line = split(line)

        if line[1] == "HETATM"
            #Stores Atom name
            info[count, 1] = line[end]
            #stores atom cartesian coordinate in angstroms
            info[count, 2:4] = line[6:8]
            count += 1

        end

        if line[1] == "CONECT"
            for j = (1:length(line) - 2)
                #this makes sense if you read the pdb format
                info[count_2, j + 1] = line[j + 1]
            end
            count_2 += 1
        end



    end
#=
    #converts the string representations of the cartesian coords from the pdb
    #into floats in a new array
    pdb_cart_coords = map(x->parse(Float64, x), info[:, 2:4])


    #stores atomic numbers to the corresponding atom in the MOF
    for i = 1:length(info[1,:])
        for j = 1:length(atoms[:atom])
            #checks if name of atom is the same as the corresponding posiiton in the feature matrix
            if info[i, 1] == atoms[:atom][j]
                unordered_MOF_info[i, 1] =  j
                break
            end
        end
    end



    #corrects atom order by comparing the cartesian coordinates of each atom in the pdb vs the framework
    #first convert framework from frac coords to cartesian
    MOF_cart_coords = (MOF.box.f_to_c * MOF.atoms.xf)'
    #Then compare those coords to the coords from the pdb to unscramble
    for (i, atom_1) in enumerate(MOF_cart_coords[1, :])
        for (j, atom_2) in enumerate(MOF_cart_coords[1, :])
            if pdb_cart_coords[i, 1:3] == MOF_cart_coords[j, 1:3]
                ordered_MOF_info[j, 1:4] = unordered_MOF_info[i, 1:4]
            end
        end
    end

    =#
    return info
end




function feature_array(MOFname::String)
    #gets info from pdb file
    lines = read_pdb(MOFname)

    #gets atom information to find name and corresponding atomic number
    atoms = CSV.read(joinpath(PATH_TO_DATA, "atom_properties.csv"))

    #creates framework of the MOF in question
    MOF = find_MOF(MOFname)

    #creates a nice array where each row corresponds to the number of the atom
    #in the MOF and the entry is the name of the atom
    MOF_info = extract_from_pdb(lines, atoms, MOF)

    #initializes feat_array
    feat_array = zeros(Float64, MOF.atoms.n_atoms, 2 * length(atoms[:atom]))

    #modifies feature vector matrix with atom identification for each atom id

    #=
    for i = 1:length(atom_ids)
        for j = 1:length(atoms[:atom])
            #checks if name of atom is the same as the corresponding posiiton in the feature matrix
            if atom_ids[i] == atoms[:atom][j]
                feat_array[i, j] = 1
                break
            end
        end
    end

    #actually finds the bonds
    for atom_1_id = 1:length(atom_ids)
        for atom_2_id = 1:length(atom_ids)

            #creates bond in feature array
            #if
            #    feat_array[atom_1_id, n + k] += 1
            #end
        end
    end
    =#
    return MOF_info
end
