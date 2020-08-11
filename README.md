# NTuple

A program that sorts the Geant4 NTuple output into histograms, matrices, etc. This program accepts a settings file that includes thresholds, resolutions, etc, for each detection system. This program was written by Vinzenz Bildstein, and later modified and edited by Evan Rand. This branch of the code is for modifications to do with TI-STAR, and will rely heavily on the older [TRexGeant4Analysis](https://github.com/VinzenzBildstein/TRexGeant4Analysis) code.

-----------------------------------------
 Installation
-----------------------------------------
The NTuple code requires the CommandLineInterface library (https://github.com/GRIFFINCollaboration/CommandLineInterface.git), please clone and compile this library first.
The Makefile assumes that the code of the CommandLineInterface library is in ~/CommandLineInterface, if this is not the case, edit the Makefile to point COMM_DIR to the directory where the code is.
The Makefile also assumes that you have a ~/lib directory where the shared-object libraries from the CommandLineInterface are installed. Again, if this is not the case, modify the Makefile so that LIB_DIR points to where these shared-object libraries are.
For the TI-STAR branch of the analysis code, the libCustomClasses shared-object library (which contains the TistarSettings class) must also be set so that the SIM_DIR variable points to the GEANT4 simulation directory.
Once the Makefile has been adjusted, the code can be compiled with a simple make.

Things to remember:
(1) set the COM_DIR variable in the Makefile to point to your CommandLineInterface installation
(2) set the SIM_DIR variable in the Makefile to point to your Griffinv10 simulation directory (where the build/include dir's are)
(3) include the directory containing the libCommandLineInterface.so library in your LD_LIBRARY_PATH
(4) include the directory containing the libCustomClasses.so library in your LD_LIBRARY_PATH (the simulation build dir)
