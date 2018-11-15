
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

    #NOTE Arrays should be changed so they are sized to exactly the number of bonds
    #This can be done with the function "insert!"

    zero_location = 0

    #creates an array of string to hold atom identifiers and coords from the pdb
    info_String = Array{String, 2}(undef, MOF.atoms.n_atoms, 14)

    #NOTE This turns array of undefs for info_String into array of strings of "0"
    #This is necessary because I want to know the length of the array later
    #as well I would like to know where information stops and undefs start in a row
    #There's got to be a better way to do this...
    for i = 1:MOF.atoms.n_atoms
        for j = 1:14
            info_String[i,j] = "0"
        end
    end

    #array to hold the cart coords once theyve been translated to floats from strings
    info_Float64 = zeros(Float64, MOF.atoms.n_atoms, 14)

    #arrays used to store atomic number, and bonding info
    unordered_MOF_info = zeros(Float64, MOF.atoms.n_atoms, 11)
    ordered_MOF_info = zeros(Float64, MOF.atoms.n_atoms, 11)

    for (i, line) in enumerate(lines)
        line = split(line)

        if line[1] == "HETATM"
        atom_id = parse(Int64, line[2])
            #Stores Atom name
            info_String[atom_id, 1] = line[end]
            #stores atom cartesian coordinate in angstroms
            info_String[atom_id, 2:4] = line[6:8]
        end

        if line[1] == "CONECT"
        atom_id = parse(Int64, line[2])

        #quick and dirty way to find the first entry with a zero in info_String[]
        #necessary because pdb fiels split bonding into two lines, and
        #the info_String matrix needs to have all bonding info in one line
        for k = 1:length(info_String[atom_id, :])
            if info_String[atom_id, k] == "0"
                zero_location = k - 1
                break
            end
        end

            for j = 1:(length(line) - 2)

                #stores the bonding info as the atoms bonded to it
                #row in info_String is the atom_id from the pdb, while row is the first
                #entry with a zero in the row
                info_String[atom_id, zero_location + j] = line[j + 2]
            end
        end
    end



    #stores atomic numbers to the corresponding atom in the MOF
    for i = 1:length(info_String[:, 1])
        for j = 1:length(atoms[:atom])
            #changes atom name to atomic number in new unordered matrix
            if info_String[i, 1] == uppercase(String(atoms[:atom][j]))
                unordered_MOF_info[i, 1] =  j
                break
            end
        end
    end


    #converts the string representations of the cartesian coords from the pdb
    #into Floats in a new array
    info_Float64 = map(x->parse(Float64, x), info_String[:, 2:14])

    unordered_MOF_info[:, 2:11] = info_Float64[:, 4:13]
    #corrects atom order by comparing the cartesian coordinates of each atom in the pdb vs the framework
    #first convert framework from frac coords to cartesian

#NOTE needs to be changed to use isapprox()

    MOF_cart_coords = (MOF.box.f_to_c * MOF.atoms.xf)'
    #Then compare those coords to the coords from the pdb to unscramble
    for (i, atom_1) in enumerate(MOF_cart_coords[:, 1])
        for (j, atom_2) in enumerate(MOF_cart_coords[:, 1])
            if isapprox(info_Float64[i, 1:3], MOF_cart_coords[j, 1:3]; atol = 0.001)
                ordered_MOF_info[j, :] = unordered_MOF_info[i, :]
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
    feat_array = zeros(Float64, MOF.atoms.n_atoms, 2 * length(atoms[:atom]))

    #creates Feature array
    for i = 1:length(ordered_MOF_info[:,1])
        #adds one to location corresponding to atomic number for each row
        feat_array[i, Int64(ordered_MOF_info[i, 1])] = 1

        #adds bonds in feature array
        for j = 2:length(ordered_MOF_info[i, :])
            #Stops loop when all bonds are added
            if ordered_MOF_info[i, j] == 0
                break
            end
            #NOTE bond added at correct position minus 1??
            feat_array[i, Int64(ordered_MOF_info[Int64(ordered_MOF_info[i, j]), 1]) + length(atoms[:atom])] += 1
        end

    end
    return feat_array
end
