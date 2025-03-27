# spatial_maturation_preterm_brain
Code related to publication on spatial maturation in the preterm brain

Four key sections

1. Run Source Space Reconstruction on a list of EDF files and generate SSR mat files - go from 128 to 50 channels
2. Extract Features from a list of SSR mat files and store features in a mat file.
     -> this process is deatiled in extract_all_features_from_file_list.m
4. Build regionally specific age predictions across all available regions
5. Analysis & Plotting
