# material_properties

Simulation results and analysis codes for the computation of rheological properties in https://doi.org/10.1038/s41567-024-02558-1

The `source_data` folder contains all experimental and computational data for the final figures.

The `trajectories` folder contains the three trajectory files for the WT^+NLS system at the simulation temperature, T = 53. For the A1-LCD, this temperature corresponds to 296.8 K.

Please also see the .m files in the `codes` folder. Note the codes now include the use of connected, rather than biconnected, graphs. The previous version of the code, as was used in the paper, used biconnected graphs, meaning that the removal of any node will not split the graph into multiple disconnected parts. Using biconnected graphs may be more robust if there are any dangling chains that may or may not actually be part of the condensate. Nonetheless, it results only in minor differences.
