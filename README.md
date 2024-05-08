# FrameFieldTracing
This project uses cross-field information, that are stored in the input meshes in the .om file format, to trace a discrete tmesh on the mesh surface.

# Building
mkdir build
cd build
cmake ..
make

# Usage
./FrameFieldTracer infile [optional: angle_treshold]

e.g. ./FrameFiledTracer ../data/spot.om 0.2

## License
This project is licensed under the BSD-3 License