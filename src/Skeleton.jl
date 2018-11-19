"""
    bgrid = binarize(grid, pot_energy_threshold)

Binarize the a potential energy grid so that values in `grid.data` below the 
`pot_energy_threshold` are assigned `true` for accessible and `false` otherwise.
Returns a new `Grid{Bool}`.

# Attributes
- `grid::Grid{Float64}`: potential energy grid to binarize
- `pot_energy_threshold::Float64`: energy value below which a grid point is considered
accessible (`true`).
# Returns
- `bgrid::Grid{Bool}`: binarized energy grid, whose `data` contain `true` for accessible, `false` for inaccessible.
"""
function binarize(grid::Grid{Float64}, pot_energy_threshold::Float64)
    return Grid(grid.box, grid.n_pts, 
                convert(Array{Bool}, grid.data .< pot_energy_threshold), 
                grid.units, grid.origin)
end

function _distance_grid(bgrid::Grid{Bool})
    # create new grid that contains the distance to the closest grid point
    # the distance is an integer that measures the distance in terms of the nb of voxels.
    dist = -1 * ones(Int, grid.n_pts...)
    for i = 1:bgrid.n_pts[1], j = 1:bgrid.n_pts[2], k = 1:bgrid.n_pts[3]
        # if not accessible, make distance zero
        if ! bgrid.data[i, j, k]
            dist[i, j, k] = 0
        end
    end
end
