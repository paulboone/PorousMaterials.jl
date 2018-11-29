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

#NOTE This code could be made stand alone for use in reading PDB files if the
#bonding bit was removed to it's own code
#Another option would be to integrate it into the framework function so
#Framework() can read pdb files.
function extract_bonds_from_pdb(lines::Array{String, 1}, atoms::DataFrames.DataFrame, MOF::Framework)

    #creates an array of string to hold atom identifiers and coords from the pdb
    info_String = fill("O", (MOF.atoms.n_atoms, 4))

    #arrays used to store atomic number, and bonding info
    unordered_MOF_info = []
    bond_dict = Dict{Int64, Array{Int64, 1}}()

    for (i, line) in enumerate(lines)
        line = split(line)

        if line[1] == "HETATM"
        atom_id = parse(Int64, line[2])
        bond_dict[atom_id] = []
            #Stores Atom name
            info_String[atom_id, 1] = line[end]
            #stores atom cartesian coordinate in angstroms
            info_String[atom_id, 2:4] = line[6:8]
        end

        if line[1] == "CONECT"
            atom_id = parse(Int64, line[2])
            for j = 3:length(line)
                append!(bond_dict[atom_id], parse(Int64, line[j]))
            end
        end
    end

    #stores atomic numbers to the corresponding atom in the MOF
    for i = 1:length(info_String[:, 1])
        #changes atom name to atomic number in info matrix
        for j = 1:length(atoms[:atom])
            if info_String[i, 1] == uppercase(String(atoms[:atom][j]))
                info_String[i, 1] =  string(j)
                break
            end
        end
    end

    #finds the max bonds an atom has in the pdb file
    #useful to set the number of columns in bonding array
    max_bonds = 0
    for i = 1:length(info_String[:, 1])
        if max_bonds < length(bond_dict[i])
            max_bonds = length(bond_dict[i])
        end
    end

    #adds bond info to an array with other info from info_String
    #done this way because it is impossible (to me) to append multidimensional arrays
    for atom_id = 1:length(info_String[:, 1])
        temp_array = []
        push!(temp_array, parse(Int64, info_String[atom_id, 1]))

        if haskey(bond_dict, atom_id)
            for k = 1:length(bond_dict[atom_id])
                push!(temp_array, (bond_dict[atom_id][k]))
            end
        end
        #needed to concatenate arrays of different sizes (adds zeros)
        while length(temp_array) != (max_bonds + 1)
            push!(temp_array, 0.0)
        end

        temp_array = temp_array'
        #only starts vertical concatenation on the second loop
        if atom_id == 1
            unordered_MOF_info = temp_array
        else
            unordered_MOF_info = vcat(unordered_MOF_info, temp_array)
        end
    end

    #converts the string representations of the cartesian coords from the pdb
    #into Floats in a new array
    info_Float64 = map(x->parse(Float64, x), info_String)

    #Corrects row order in pdb to match the order in the Framework
    MOF_cart_coords = (MOF.box.f_to_c * MOF.atoms.xf)'
    new_row = Dict{Int64, Int64}()
    ordered_MOF_info = zeros(Int64, MOF.atoms.n_atoms, max_bonds + 1)
    for i = 1:length(MOF_cart_coords[:, 1])
        for j = 1:length(MOF_cart_coords[:, 1])
            #corrects row order
            if isapprox(info_Float64[i, 2:4], MOF_cart_coords[j, 1:3]; atol = 0.001)
                ordered_MOF_info[j, :] = unordered_MOF_info[i, :]
                new_row[i] = j
            end
        end
    end

    #corrects bonds because they point to rows, and rows have changed
    for i = 1:length(ordered_MOF_info[:, 1])
        for j = 2:length(ordered_MOF_info[1, :])
            if ordered_MOF_info[i,j] != 0
                ordered_MOF_info[i,j] = new_row[ordered_MOF_info[i,j]]
            end
        end
    end

    return ordered_MOF_info
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
    ordered_MOF_info = extract_bonds_from_pdb(lines, atoms, MOF)

    #initializes feat_array
    feat_array = zeros(Int64, MOF.atoms.n_atoms, 2 * length(atoms[:atom]))

    #creates Feature array
    for i = 1:length(ordered_MOF_info[:, 1])
        #adds one to location corresponding to atomic number for each row
        feat_array[i, Int64(ordered_MOF_info[i, 1])] = 1

        #adds bonds in feature array
        for j = 2:length(ordered_MOF_info[i, :])
            #Stops loop when all bonds are added
            if ordered_MOF_info[i, j] == 0
                break
            end
            feat_array[i, Int64(ordered_MOF_info[Int64(ordered_MOF_info[i, j]), 1]) + length(atoms[:atom])] += 1
        end

    end

    return feat_array
end
