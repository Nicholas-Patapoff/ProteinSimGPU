# ProteinSimGPU

**ProteinSimGPU** is a molecular dynamics (MD)-style simulation engine for proteins with both **CPU** and **GPU (CUDA)** backends.  
It supports **AMBER-style input files** (`.pdb`, `.parm7`, `.crd`) and produces coordinate trajectories (`.crd`) that can be visualized or analyzed further.

This project demonstrates how GPU acceleration can significantly improve protein simulation performance compared to CPU-only runs.

---

## Features
- ⚡ **GPU acceleration** with CUDA (dynamic parallelism enabled)
- 🖥️ **CPU fallback** for reference and debugging
- 📦 Support for AMBER-style **PDB/parm7/CRD** inputs
- 🧪 Example protein data included (`Data/` folder)
- 📊 Outputs trajectories in `.crd` format

---

## Requirements
- **C++20 compiler** (e.g., `g++`)
- **CUDA Toolkit** (tested with CUDA 11+)
- **NVIDIA GPU** with compute capability **8.6** (Ampere, e.g., RTX 30-series)  
  > To use a different GPU architecture, edit `Makefile`:
  ```make
  NVFLAGS=-arch=sm_86   # adjust for your GPU
  ```

---

## Building
Clone and build with `make`:

```bash
git clone https://github.com/Nicholas-Patapoff/ProteinSimGPU.git
cd ProteinSimGPU
make
```

This produces an executable:

```
myApp
```

---

## Usage

### Example run (CPU + GPU)
The entry point (`main.cpp`) demonstrates running both CPU and GPU simulations:

```cpp
Environment coords("Data/chignolin_ext.pdb");
parm data("Data/chig_ext.parm7");

// CPU simulation
cpu_simulation small_probe(coords, data, 1, "test_prot");
small_probe.update_coord(0.01, 10000, 1000);

// GPU simulation
gpu_sim(10000, 1000, 0.01, ready_data(coords, data));
```

### Parameters
- **Steps**: total number of timesteps (`10000`)
- **Print Frequency**: output interval (`1000`)
- **Timestep**: integration timestep (`0.01`)

---

## Inputs & Outputs
- **Inputs**
  - `.pdb` — Protein structure (atom positions)
  - `.parm7` — Force field / simulation parameters
  - `.crd` (optional) — Starting coordinates

- **Outputs**
  - `.crd` trajectory file (e.g., `gpuOUT.crd`)
  - Printed runtime information (`Completed CPU in X seconds!`)

Example input files are provided in:
- `Data/` (chignolin test protein)
- Top-level test files (`test_prot.pdb`, `test_prot.parm7`, `test_prot.crd`)

---

## Example Run
After building:

```bash
./myApp
```

This will:
1. Run a CPU simulation of `chignolin_ext.pdb` with `chig_ext.parm7`
2. Output timing results to console
3. Run the equivalent GPU simulation
4. Produce a `.crd` trajectory (e.g., `gpuOUT.crd`)

---

## Development Notes
- Source files:
  - `CPU_Simulation.cpp` — CPU implementation
  - `GPU_sim.cu`, `GPU.h`, `GPU_Data.*` — GPU simulation
  - `Vector_Math.*` — math utilities
  - `parmParse.cpp`, `parm.h` — parameter parsing
  - `Env.h`, `Environment.cpp` — environment setup
- Configurable through `main.cpp` (hardcoded demo example included)
- Extendable for different proteins by changing PDB/PARM inputs

---

## License
MIT License (or specify if different).

---

## Author
Developed by **Nicholas Patapoff**.  
For inquiries, please open an issue on the [GitHub repository](https://github.com/Nicholas-Patapoff/ProteinSimGPU).
