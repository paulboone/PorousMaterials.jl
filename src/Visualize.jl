using Makie, GeometryTypes, Colors

struct FrameImage
    frame::Framework
    atomic_radii::Dict{Symbol, Float64}
    atom_colors::Dict{Symbol, Tuple{Int, Int, Int}}
    scene::Makie.Scene{:makie}
end

function construct_frame_image(frame::Framework, resolution = (500, 500))
    atom_colors = read_cpk_colors()
    atomic_radii = read_atomic_radii()
    scene = Scene(resolution = resolution)
    return frameimg = FrameImage(frame, atomic_radii, atom_colors, scene)
end

function plot_atoms(frimg::FrameImage)
    coords = frimg.frame.box.f_to_c * frimg.frame.xf
    for atom in unique(frimg.frame.atoms)
        meshscatter(coords[1, frimg.frame.atoms .== atom],
                    coords[2, frimg.frame.atoms .== atom],
                    coords[3, frimg.frame.atoms .== atom],
                    color = RGB(frimg.atom_colors[atom][1]/255, 
                                frimg.atom_colors[atom][2]/255, 
                                frimg.atom_colors[atom][3]/255),
                    markersize = 2*frimg.atomic_radii[atom])
    end
    center!(frimg.scene)
end

function add_bonds(frimg::FrameImage)
    N = frimg.frame.n_atoms
    add_bond = falses(N,N)
    dist = fill(Inf32, (N,N))
    for i in 1:N
        for j in i+1:N
            dist[i,j] = norm(frimg.frame.box.f_to_c * (frimg.frame.xf[:,i] - frimg.frame.xf[:,j]))
            if dist[i,j] < 2
                add_bond[i,j] = true
            end
        end
    end
    print(typeof(dist))
    meshC = GeometryTypes.GLNormalMesh(GeometryTypes.Cylinder{3, Float32}(
                                       GeometryTypes.Point3f0(0., 0., 0.),
                                       GeometryTypes.Point3f0(0., 0., 1.),
                                       Float32(1)), 30)

    for i = 1:N
        atoms_coord = frimg.frame.box.f_to_c * frimg.frame.xf[:,add_bond[i,:]]
        origin = frimg.frame.box.f_to_c * frimg.frame.xf[:,i]
        n_bond = size(atoms_coord, 2)
        if n_bond < 1
            continue
        end
        sizesC = [Vec3f0(0.1, 0.1, dist[i, add_bond[i,:]][j]) for j = 1:n_bond]
        Qlist = zeros(n_bond, 4)
        for k = 1:n_bond
            ct = GeometryTypes.Cylinder{3, Float32}(
                         GeometryTypes.Point3f0(origin[1], origin[2], origin[3]),
                         GeometryTypes.Point3f0(atoms_coord[1,k], atoms_coord[2,k], atoms_coord[3,k]), Float32(1))
            Q = GeometryTypes.rotation(ct)
            #r = sqrt(1 + Q[1, 1] + Q[2, 2] + Q[3, 3])
            r = dist[i, add_bond[i,:]][k]
            Qlist[k, 4] = r
            Qlist[k, 1] = (Q[3,2] - Q[2,3]) / (4 * r)
            Qlist[k, 2] = (Q[1,3] - Q[3,1]) / (4 * r)
            Qlist[k, 3] = (Q[2,1] - Q[1,2]) / (4 * r)
        end
        rotationC = AbstractVector[Vec4f0(Qlist[i, 1], Qlist[i, 2], Qlist[i, 3], Qlist[i, 4]) for i = 1:n_bond]
        pG = [GeometryTypes.Point3f0(atoms_coord[1,i], atoms_coord[2,i], atoms_coord[3,i]) for i = 1:n_bond]
        println(pG)
        @printf("\n\n\n")
        meshscatter(pG, marker = meshC, color = :grey,  markersize = sizesC, rotations = rotationC)

    end
    return add_bond, dist
end
