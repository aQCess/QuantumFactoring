# Quantum factoring algorithm using Grover search

A code for implementing quantum circuits to factor an $n$-bit integer $N$ using $2n-5$ qubits, based on the algorithm described in *S. Whitlock and T. D. Kieu, "Quantum factoring algorithm based on Grover search", arXiv reference to be confirmed.*

The algorithm doesn’t depend on any properties of the number to be factored, has guaranteed convergence, and doesn’t require complex classical pre or post-processing.

**Note** This open source code is built using the free library to compose and manage quantum circuits: [MIMIQ](https://github.com/qperfect-io/MimiqCircuits.jl) by [QPerfect](https://www.qperfect.io/). To run the algorithm on the MIMIQ cloud quantum simulator you can request a short trial or subscription at https://www.qperfect.io/.

## COPYRIGHT

Copyright © 2023 University of Strasbourg

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

## SETUP

These instructions explain the steps involved in getting the code to run on linux.

0. Have Julia installed. The code has been tested in Julia v1.9. If you are installing from scratch I recommend using the [Juliaup](https://github.com/julialang/juliaup#installation) version management tool.

0. Download this repository to a local folder, or if using [git](https://git-scm.com/):

    `git clone https://github.com/aQCess/QuantumFactoring.git`

0. Install MIMIQ circuits library (Julia version)

	```julia
	julia> ]
	(v1.9) pkg> registry update
	(v1.9) pkg> registry add https://github.com/qperfect-io/QPerfectRegistry.git
	(v1.9) pkg> add MimiqCircuits
	```

## USAGE
Tested with `MimiqCircuits.jl v0.8.0` 

### *Full algorithm (requires access to the MIMIQ cloud quantum simulator)*
```Julia
julia> using MimiqCircuits
julia> include("factor_grover_2n-5.jl")

# establish a connection to MIMIQ cloud (requires credentials)
julia> conn = connect(; url=QPERFECT_CLOUD)
julia> factorize(1073, conn; compiled=false)
```
```Julia
Result :			 37 x 29 = 1073
Trials :			 1
Success probability :		 1.0 (100/100 samples)
Simulation time :		 0.798278261 seconds
Size of the Hilbert space :	 2^17
Grover iteration steps :	 8
Number of one-qubit gates :	 7499
Number of two-qubit gates :	 12842
```

### *Step-by-step*
```julia
julia> using MimiqCircuits
julia> include("factor_grover_2n-5.jl")

julia> N = 1073 # a number to factorize
julia> s = 1 # guess for s = +/- 1
julia> K = 8 # number of Grover iteration steps

# obtain the oracle polynomial
julia> f = foracle(N, s);

# define three quantum registers
julia> regs = ([1:3...], [4:7...], [8:17...]);
```

```julia
# build the circuit (decomposed)
julia> circ = grover_loop_decomp(f, K, regs...)
```
```julia
17-qubit circuit with 20341 instructions:
├── H @ q[1]
├── H @ q[2]
├── H @ q[3]
├── H @ q[4]
├── H @ q[5]
├── H @ q[6]
├── H @ q[7]
├── H @ q[17]
├── CRZ(π/2) @ q[16], q[17]
⋮   ⋮
├── H @ q[14]
├── CRZ(-1π/2) @ q[14], q[15]
├── CRZ(-1π/4) @ q[14], q[16]
├── CRZ(-1π/8) @ q[14], q[17]
├── H @ q[15]
├── CRZ(-1π/2) @ q[15], q[16]
├── CRZ(-1π/4) @ q[15], q[17]
├── H @ q[16]
├── CRZ(-1π/2) @ q[16], q[17]
└── H @ q[17]
```

```julia
# establish a connection to MIMIQ cloud (requires credentials)
julia> conn = connect(; url=QPERFECT_CLOUD);

# execute circuit on MIMIQ cloud
julia> res = execute_circuit(conn, circ);

# obtain the factors
julia> factors = get_factors(histsamples(res), f, regs)
julia> println("Factors of $N are $(factors[1][1:2])")
```
```julia
Factors of 1073 are (37, 29)
```
