# Resummation-Based SSE QMC for Higher Symmetric Continuous 'N' Heisenberg Models

This repository contains the implementation of a **Resummation-Based Stochastic Series Expansion (SSE) Quantum Monte Carlo (QMC)** algorithm for simulating **higher symmetric continuous 'N' Heisenberg models**. The code is designed to efficiently handle systems with SU(N) symmetry, leveraging resummation techniques to enhance computational efficiency and accuracy.

---

## Features

- Implementation of the SSE QMC method adapted for SU(N) Heisenberg models.
- Support for varying lattice geometries, system sizes, and interaction strengths.
- Enhanced performance using resummation techniques for high-order terms in the expansion.
- Calculation of key observables, including energy, correlation functions, and symmetry properties.

---

## Getting Started

### **Prerequisites**

Before using this code, ensure that the following are installed on your system:

- **C++ compiler** supporting C++17 or later (e.g., GCC, Clang)
- **CMake** for build automation (recommended)
- Libraries like **Boost** or **Eigen** (if required in the code).

### **Installation**

Clone the repository to your local machine:

```bash
git clone https://github.com/YourUsername/resummation-sse-qmc.git
cd resummation-sse-qmc
