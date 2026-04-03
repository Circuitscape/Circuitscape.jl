# Computational Limitations, Speed, and Landscape Size

We have tested this code on landscapes with up to 437 million cells. Increasing numbers of connections using diagonal (eight neighbor) connections will decrease the size of landscapes that can be analyzed. Also, increasing landscape size or numbers of focal nodes will increase computation time. Note that due to the matrix algebra involved with solving many pairs of focal nodes, Circuitscape will run much faster when focal points (each focal node falls within only one grid cell), rather than focal regions (at least one focal node occupies multiple grid cells), are used.

## Memory Limitations

There are several ways to increase the solvable grid size:

- Set impermeable areas of your resistance map to NODATA
- Use focal points instead of regions in pairwise mode
- Connect cells to their four neighbors only (`connect_four_neighbors_only = True`)
- Disable current and voltage maps (`write_cur_maps = False`, `write_volt_maps = False`)
- Use the one-to-all or all-to-one modes, which typically use less memory and run more quickly than pairwise mode
- Use the `cg+amg` solver instead of `cholmod` for large grids (CHOLMOD uses significantly more memory)
- Coarsen your grids (use larger cell sizes) -- this often produces qualitatively similar results (see McRae et al. 2008)

The all-to-one mode can be an alternative to pairwise mode when the goal is to produce a cumulative map of important connectivity areas among multiple source/target patches.
