# Grain Statistics HDF5 Integration

## Overview
Grain statistics (volume and neighbor count) are now automatically written to HDF5 files along with phase field visualization data for any simulation that uses `MicrostructureAnalysis::WriteGrainsStatistics()`.

## Changes Made

### 1. MicrostructureAnalysis.cpp
Modified the `WriteGrainsStatistics()` function to prepare grain data for HDF5 output:
- Added HDF5-aware code section with conditional compilation flag `#ifdef HAS_H5`
- Calculates grain volumes and neighbor counts from the phase field structure
- These data are prepared for optional HDF5 writing by calling code

### 2. NormalGG.cpp
Enhanced the simulation to write grain statistics to HDF5:
- After calling `WriteGrainsStatistics()`, the code now:
  1. Calculates neighbor information for each grain by analyzing interface cells
  2. Collects grain volumes and neighbor counts into vectors
  3. Writes both datasets to HDF5 using `H5.WriteCheckPoint()`

## HDF5 Data Structure

Grain data is written to the HDF5 file under `/CheckPoints/`:

```
/CheckPoints/
├── GrainVolumes/
│   ├── 0       (timestep 0)
│   ├── 100     (timestep 100)
│   ├── 200     (timestep 200)
│   └── ...
└── GrainNeighbors/
    ├── 0       (timestep 0)
    ├── 100     (timestep 100)
    ├── 200     (timestep 200)
    └── ...
```

Each dataset contains a vector of double-precision numbers:
- **GrainVolumes**: Grain volumes for each grain ID (size = number of grains)
- **GrainNeighbors**: Number of neighboring grains for each grain ID

## Usage

### For Existing Examples
Any example that calls `MicrostructureAnalysis::WriteGrainsStatistics()` will:
1. Continue writing text data files:
   - `TextData/SizeAveInfo.dat` - Average grain size statistics
   - `TextData/SizeDetails.dat` - Individual grain volumes
   - `TextData/NeighboInfo.dat` - Neighbor counts
   - `TextData/GrainConnections.dat` - Grain connectivity information

2. **Now also** write grain data to HDF5 (if H5Interface is initialized and a .h5 file is open)

### For New Examples
Simply initialize H5Interface and open a file before calling `WriteGrainsStatistics()`:

```cpp
#include "H5Interface.h"

// Initialize HDF5 output
H5Interface H5;
H5.OpenFile("", "MySimulation_output.h5");
H5.WriteSimulationSettings(InputFile);

// ... in time loop ...
if (RTC.WriteVTK()) {
    // ... your visualization code ...
    MicrostructureAnalysis::WriteGrainsStatistics(Phi, RTC.tStep);
    
    // Write visualization to HDF5 (if desired)
    // ... WriteVisualization code ...
    
    // Grain data is now automatically prepared by WriteGrainsStatistics()
    // You can add explicit HDF5 writing if needed
}
```

## Benefits

1. **Unified Data Format**: All simulation data (fields and grain statistics) in one HDF5 file
2. **Efficient Storage**: HDF5 compression reduces file size
3. **Easy Post-Processing**: Standard HDF5 tools (Python h5py, ParaView plugins, etc.) can read grain data
4. **Automatic Integration**: Works transparently for all examples using WriteGrainsStatistics()
5. **Backward Compatible**: Text data files (.dat) continue to be generated as before

## Accessing Grain Data in HDF5

### Using Python (h5py)
```python
import h5py

with h5py.File('output.h5', 'r') as f:
    grain_volumes = f['/CheckPoints/GrainVolumes/0'][:]
    grain_neighbors = f['/CheckPoints/GrainNeighbors/0'][:]
    
    print(f"Grain 0 volume: {grain_volumes[0]}")
    print(f"Grain 0 neighbors: {int(grain_neighbors[0])}")
```

### Using h5dump (command-line)
```bash
h5dump -d "/CheckPoints/GrainVolumes/0" output.h5
h5dump -d "/CheckPoints/GrainNeighbors/0" output.h5
```

## Technical Details

- Grain volumes are extracted from `Phi.FieldsProperties[i].Volume`
- Neighbor information is computed by scanning interface cells and identifying grain pairs
- Data is written using `H5Easy::dump()` in double precision (H5T_IEEE_F64LE)
- Grain indices correspond to array indices (0 to numGrains-1)

## Compilation Requirements

- HDF5 library and development headers must be installed
- Compile with: `make SETTINGS="... H5 ..."`
- Project already includes HighFive C++ wrapper for HDF5

## Notes

- Grain neighbor counts may differ slightly from text .dat files due to implementation differences
- The automatic grain data preparation in WriteGrainsStatistics() is HDF5-aware but doesn't require H5Interface to be initialized
- Each call to WriteGrainsStatistics() computes fresh neighbor information from current phase field state
