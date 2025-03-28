# spatial_maturation_preterm_brain
Code related to publication on spatial maturation in the preterm brain

Four key sections

1. Read EDF files and perform artefact removal and store in mat format
2. Do source space reconstruction (go from 125 to 58 channels) and extract features from a list of mat files.

     -> this process is detailed in extract_all_features_from_file_list.m
   
3. Build regionally specific age predictions across all available regions

     -> this process is details in do_all_iteration_for_github.m and will use the inputs in the all_data.mat file soted in the supporting_mat_files directory, this contains EEG features and PMA, not raw EEG recordings.
     -> to do type the following commands into matlab >> load all_data; [r, ca, cb, outs1, outs2, err0, err1, err2, err3, err4, val_ref, ppx, out1, out2, out4, df] = do_iteration_all_for_github(fv2, pma2, id2, demos, MyAtlas, 0, fref);
  
4. Analysis & Plotting
