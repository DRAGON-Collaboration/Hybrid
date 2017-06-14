# DRAGON Hybrid Detector
GEANT4 source code for an application simulating a hybrid ionization chamber / double-sided-silicon
strip detector to operate as a dE vs E detector for heavy ions within astrophysical energy ranges.

Simulating a Hybrid Ionization Chamber / Double-Sided-Silicon-Strip Detector
Author: Devin Burke  Email: burkeds@mcmaster.ca
McMaster University
Department of Physics and Astronomy

**Geant4 Compatibility: Geant4.10**

The application is written in C++ and uses the GEANT4 simulation toolkit as well as standard
C++ libraries. This application builds from previous work done by Dario Gigliotti and from the
GEANT4 example application TESTEM7. 
The source code can be compiled using cmake in the standard way most GEANT4 applications can be compiled (i.e. - performing an "out of source" build as is described [here](http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/InstallationGuide/html/ch03s02.html)):

```
mkdir build
cd build
cmake -DGEANT4_DIR=$INSTALL_DIR$ $PATH_TO_SOURCE$ make -jn
```

where `n` is the number of cores you wish to use.

Commands can be entered from the UI or from batch. To run a macro from batch enter:

`./Si_Ion_Chamber_v8.2 $MACRO$ > $G4cout OUTPUT FILE$`

To run a macro from the UI use the command:
`/control/execute $MACRO$`
			
The application opens in the pre-initialization stage. Here the user can set detector geometries and material
types. 
Type `ls` to see the list of commands.

App-specific commands are located in the `testem` command directory. Associated help commands can be entered
to explain their function and a list can be seen in the source code file DetectorMessenger.cc. All UI 
commands are defined within their associated messenger source code files.

All materials and geometry can be defined from the UI commands provided the materials are defined in
DetectorConstruction.cc with the exception of the DSSSD dead layer which must currently be altered 
in the source code.

Isobutane pressure can be set from 1 to 31 torr by entering the command:
`/testem/det/setGasMat isobutane$Ptorr`
where `$P` is an integer from 1 to 31.
There is nothing stopping the user from defining other pressures within the source code. Just add the
material in DetectorConstruction.cc.

The length of the gas chamber along the beam axis can be set by:
`/testem/det/setGasThick`

The sensitive regions of the gas can be set by three commands:
```
/testem/det/setAnodePosition	//sets the origin of the sensitive region
/testem/det/setAnodeLength		//sets the region length along the beam axis
/testem/det/setSegmentLength	//divides the length into equal segments with dead regions
								//if the lengths do not divide evenly into the total length.
```

Once detector parameters are set the user should enter the initialized stage.
`/run/initialize`

From here the user can specify the parameters of the general particle source (GPS). The GPS is unaltered
from the GEANT4 tool and the user can see the GEANT4 documentation for instruction in its use however
I have included an example below.

Specifying a beam of carbon ions
```
/gps/particle ion    	//specify an ion type projectile
/gps/ion 6 12 0			//The ion has Z=6	A=12 Q=0
/gps/ene/type Gauss		//The beam energy spectrum is of gaussian type, each ion is sampled from a Gaussian profile.
/gps/ene/mono 10 MeV	//The gaussian profile peaks at 10 MeV
/gps/ene/sigma 1 MeV 	//The gaussian profile has a STDEV of 1 MeV
/gps/pos/type Beam		//Specifies a beam type position distribution
/gps/pos/sigma_r 2 mm	//Set STDEV in radial of the beam positional profile
/gps/pos/rot1 0 0 1		//Set the 1st vector defining the rotation matrix
/gps/ang/type beam1d	//Specifies a 1D beam angular distribution
/gps/ang/sigma 10 mrad	//Set STDEV in the angular distribution 
/gps/ang/rot1 0 0 1
/run/beamOn 1000		//Begins a run with 1000 events
```

Simulation results are displayed in the GEANT4 window at the end of the run. Results are also recorded in text
files within the application directory as well as a ROOT tree. 
Efficiency is calculated and displayed in the GEANT4 window. It is the ratio of event hits to total events. 
In this application an event is considered a hit when the primary event stops within the active region of the DSSSD
volume.

If you have any questions you can email me. Visualization should work however it crashes sometimes and not others 
and I don't know why. Visualization has a big performance impact and I do not recommend you use it when running
long simulations.