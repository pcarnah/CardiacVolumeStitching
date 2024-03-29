# CardiacVolumeStitching

Slicer modules used for generated multi-view compounded volumes from TEE images. Requires 1 or more mid-esophageal views and 4 or more transgastric views with approximately 80% overlap between volumes.

Complete package can be found in Releases.

Citation: Carnahan, P.; Moore, J.;
Bainbridge, D.; Chen, E.C.S.; Peters,
T.M. Multi-View 3D Transesophageal
Echocardiography Registration and
Volume Compounding for Mitral
Valve Procedure Planning. Appl. Sci.
2022, 12, 4562. https://doi.org/10.3390/app12094562

This module is configured as a Superbuild which includes dependencies Elastix and SimpleElastix. Tested on Slicer version 4.13.0 2022-04-24.

### Usage

Load the acquired sequences into Slicer (Transgastric and mid-esophogeal as described in paper). From Cardiac Volume Stitching module, select input sequences using the dropdown and + and - buttons. The position of the acquisitions can be set in the volume list. Select the sequence to use as the master sequence for the purpose of the timing of the sequence frames. Select or create a sequence to hold the output, and click apply to proceed with proccessing.

<div align="center">
    <img src="/Screenshot.PNG" width="800px"</img> 
</div>

### Build instructions

Create a build folder to configure cmake project into. Do not build in-source. Configure CMake from the build folder with Slicer_DIR configured to location of build of 3D Slicer.

```
cmake -G "Visual Studio 17 2022" -A x64 -DSlicer_DIR:PATH="C:/D/Slcr/S4R/Slicer-build" ../CardiacVolumeStitching
cmake --build . --config Release
```

To install, navigate into the inner-build directory, optionally configure CMake with the desired install location, and build the INSTALL target.
```
cmake -DCMAKE_INSTALL_PREFIX="%BUILD_FOLDER%/SER-install" .
cmake --build . --config Release --target INSTALL
```

A compressed package that can be installed via the Slicer extension manager can also be created using the PACKAGE target
```
cmake --build . --config Release --target PACKAGE
```
