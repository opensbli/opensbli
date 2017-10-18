## Removed from metrics.py
def create_metric_bc_kernel(metrics, detJ, direction, side, block):
    # if block.ndim == 1:
    #     transformed_metrics = -1*metrics
    #     transformed_det = -1*detJ
    # else:
    #     transformation_matrix = eye(block.ndim)
    #     transformation_matrix[direction, direction] = -1
    #     transformed_metrics = metrics*transformation_matrix.adjugate()
    #     transformed_det = transformation_matrix.det()*detJ

    kernel = Kernel(block, computation_name="Metric boundary dir%d side%d" % (direction, side))
    kernel.set_boundary_plane_range(block, direction, side)
    kernel.ranges[direction][side] = block.ranges[direction][side]
    halos = kernel.get_plane_halos(block)
    # Add the halos to the kernel in directions not equal to boundary direction
    for i in [x for x in range(block.ndim) if x != direction]:
        kernel.halo_ranges[i][0] = block.boundary_halos[i][0]
        kernel.halo_ranges[i][1] = block.boundary_halos[i][1]
    # # Add halos across the boundary side only
    # kernel.halo_ranges[direction][side] = block.boundary_halos[direction][side]

    loc = [0 for i in range(block.ndim)]
    new_loc = loc[:]
    if side == 0:
        from_side_factor = -1
        to_side_factor = 1
    elif side == 1:
        from_side_factor = 1
        to_side_factor = -1

    lhs = [x for x in metrics if x not in [-1, 0, 1]] + [detJ]
    rhs = lhs[:]
    # rhs = [x for x in transformed_metrics if x not in [-1, 0, 1]] + [transformed_det]
    lhs = block.dataobjects_to_datasets_on_block(lhs)
    rhs = block.dataobjects_to_datasets_on_block(rhs)

    new_loc[direction] += to_side_factor

    transfer_indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, abs(halos[direction][side]) + 1)]

    final_equations = []
    for index in transfer_indices:
        array_equations = []
        loc_lhs, loc_rhs = loc[:], loc[:]
        loc_lhs[direction] += index[0]
        loc_rhs[direction] += index[1]
        for left, right in zip(lhs, rhs):
            left = convert_dataset_base_expr_to_datasets(left, loc_lhs)
            right = convert_dataset_base_expr_to_datasets(right, loc_rhs)
            array_equations += [Eq(left, right, evaluate=False)]
        final_equations += array_equations
    kernel.add_equation(final_equations)
    pprint(kernel.equations)
    kernel.update_block_datasets(block)
    return kernel

## Removed from metrics.py
    def apply_boundary_conditions(cls, block):
        """No global boundary conditions for the metric terms"""
        # ndim = block.ndim
        # directions = [i for i in range(ndim)]
        # metrics = Matrix(ndim, ndim, cls.FD_metrics[:])
        # detJ = cls.detJ
        # kernels = []
        # arrays = metrics[:] + [detJ]
        # arrays = block.dataobjects_to_datasets_on_block(arrays)
        # pprint(arrays)
        # kernels = block.apply_boundary_conditions(arrays)
        # for direction in directions:
        # for side in [0,1]:
        # kernels += [create_metric_bc_kernel(metrics, detJ, direction, side, block)]
        # cls.Kernels += kernels
        # exit()
        return