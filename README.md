'''Scripts to produce the ascot input files starting from BioSaw files '''<br><br>
The Toroidal magnetic field must include both the assialsymmetric component and the non-assialsymmetric component. <br>
Br and Bz are only the perturbation <br>

1) produce B field using biosaw
2) use biosaw2h5.m to convert the *.bpert file into an h5 file (let's say output_biosaw.h5). Here you have the R,phi,Z grid, the coils position and the B field.
  2a) if needed, use function merge_biosaw from pybiosaw.plot_biosawoutput to merge EFCC and TF
  2b) from pybiosaw.plot_biosawoutput you can plot TF and EFCC coils location

in biosaw2ascot folder: <br>
3) use pybiosaw.ascotinput_from_biosawoutput files to read output_biosaw.h5 and produce the file in the correct ascot5 structure (input_toascot.h5). Here you have R,phi,Z grid and B field<br>
4) you can check the resulting fields into input_toascot.h5 using pybiosaw.check_ascot_input functions<br>
5) now you can copy the "bfield" field of input_toascot.h5 into the ascot.h5 file you want to use<br>
