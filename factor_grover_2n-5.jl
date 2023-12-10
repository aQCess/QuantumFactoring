#
# Copyright © 2023 University of Strasbourg.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

""" 
Quantum factoring algorithm using Grover search

Factor an n-bit integer N using 2n-5 qubits.

Algorithm described in S. Whitlock & T. D. Kieu, 
"Quantum factoring algorithm using Grover search"
preprint arXiv:TBC
"""

# import some useful functions
include("_decompositions.jl")
include("_report.jl")
include("_registers.jl")

""" 
    foracle(N::Int)

returns a dictionary of coefficients of the polynomial f(x,y)
such that f(a,b) == 0, with an entry for each value of s=+/-1
"""
function foracle(N::Int, s::Int)
    @assert abs(s)==1

    S = floor(Int, 1.5-0.5*mod(N,6))
    p1,p0 = 6, 6+s   # 6(a+1)+s
    q1,q0 = 6, 6+s*S # 6(b+1)+s*S 

    # calculate polynomial coefficients
    A = div(p1*q1, 6)   # ab
    B = div(p1*q0, 6)   # a
    C = div(p0*q1, 6)   # b
    D = div(p0*q0-N, 6) # offset
    
    return (A,B,C,D)
end


"""
Quantum factoring algorithm 

Calls grover_loop() or grover_loop_decomp() if compiled = false

### Arguments
* N::Integer Number to be factorized N = p*q

### Keyword arguments
* dn::Integer        Estimate of the bit length difference between the factors p,q. 
                     If not specified the algorithm will perform multiple trials with different dn
* s::Integer         Guess parameter s=+1 or s=-1
* compiled::Bool     Run compiled version of the circuit (default true). 
                     Makes the simulation much slower if false
* K::Integer         Number of Grover iterationsteps (default K_opt)
* nsamples::Integer  Number of samples to take of the final state (default 100) 
* pq::Tuple{Int,Int} Hint for factors (p,q) used to compute optimal register size (default nothing)
"""
function factorize(N::Integer, conn::Connection; 
                    s=[1,-1], dn=0:floor(Int, log2(N))-5, compiled=true, 
                    K=nothing, nsamples=100, pq=nothing)
    
    # check for trivial factors
    if N%2==0
        return FactorReport("p"=>2, "q"=>div(N,2), "alg"=>"N%2", "P"=>1.0)
    elseif N%3==0
        return FactorReport("p"=>3, "q"=>div(N,3), "alg"=>"N%3", "P"=>1.0)
    end

    # otherwise go ahead and run factoring algorithm
    trials = 0

    for d in dn, s_ in s
        # oracle polynomials f(a,b) == 0 for s+/-1
        f = foracle(N,s_) 

        if pq isa Tuple{Int,Int} 
            qregs,s_ = registers(pq)
        else
            qregs = registers(N, d)
        end
        if qregs == nothing
			break
		end
        nx, ny, nz = length.(qregs)

        for nsols in 1:2

            # number of iteration steps for convergence
            if K == nothing
                K_ = floor(Int, pi/4*sqrt(2^(nx+ny)/nsols))
            else
                K_ = K
            end
            
            # build circuit, first assuming nsols = 1
            if compiled
                circ = grover_loop(f, K_, qregs[1:2]...)
            else
                circ = grover_loop_decomp(f, K_, qregs...)
            end

            # run simulation and obtain sample results
            res = execute_circuit(conn, circ; nsamples=nsamples)
            trials += 1

            # get factors p,q
            factors = get_factors(histsamples(res), f, qregs)

            if length(factors) == 0
                break #factors not found
            end
            p,q = factors[1][1:2]
            occ = sum([f[3] for f in factors])
        
            # return results if solution found
            if p*q == N && length(factors) == nsols

                results = ( "N"=>N, "p"=>p, "q"=>q,
                            "P"=>occ/nsamples,
                            "trials"=>trials,
                            "alg"=>"grover",
                            "compiled"=>compiled, 
                            "tapply"=>res.timings["apply"], 
                            "n1gates"=>count(i->numqubits(i)==1, circ),
                            "n2gates"=>count(i->numqubits(i)==2, circ),
                            "steps"=>K_, 
                            "nqubits"=>numqubits(circ), 
                            "occurances"=>occ,
                            "nsamples"=>nsamples,
                            "nreg"=>length.(qregs),
                            "bitstring"=>factors[1][4],
                            "nsols"=>nsols,
                            "s"=>s_ )

                return FactorReport(results...)
            end
        end
    end
    return FactorReport()
end

""" 
    execute_circuit(conn, circ)

execute circuit on MIMIQ cloud and return results

### Arguments
* conn::Connection MIMIQ connection
* circ::Circuit Circuit to simulate
"""
function execute_circuit(conn::Connection, circ::Circuit; nsamples::Integer=100)

    # evaluate the circuit on MIMIQ cloud
    job = execute(conn, circ, nsamples = nsamples)
    res = getresults(conn, job; interval=0.1)

    # return samples
    return res
end
"""
Grover loop using compiled oracle. Does not use the zreg

### Arguments
* foracle::NTuple{4,Int} polynomial coefficients for oracle function
* K::Int Number of interation steps
* xreg::Vector{Integer} Indices of the x register qubits
* yreg::Vector{Integer} Indices of the y register qubits
"""
function grover_loop(foracle::NTuple{4,Int}, K::Int,
                    xreg::Vector{Int}, yreg::Vector{Int})
    
    nx = length(xreg)
    ny = length(yreg)

    # build circuit
    circ = Circuit()

    # initialize registers
    for j in 1:nx+ny
        push!(circ, GateH(), j)
    end

    for k in 1:K
        # oracle
        push!(circ, PolynomialOracle(nx,ny,foracle...), xreg..., yreg...)

        # diffusion operator
        push!(circ, Diffusion(nx+ny), xreg..., yreg...)
    end
    return circ

end

"""
Grover loop decompiled to elementary one and two qubit gates

### Arguments
* foracle::NTuple{4,Int} Polynomial coefficients for oracle
* K::Int Number of interation steps
* xreg::Vector{Integer} Indices of the x register qubits
* yreg::Vector{Integer} Indices of the y register qubits
* zreg::Vector{Integer} Indices of the z register qubits
"""
function grover_loop_decomp(foracle::NTuple{4,Int}, K::Int,
                            xreg::Vector{Int}, yreg::Vector{Int}, zreg::Vector{Int})

    nx = length(xreg)
    ny = length(yreg)
    nz = length(zreg)

    # special operations
    phimult = phimultiply(foracle..., xreg, yreg, zreg)
    inv_phimult = inverse(phimult)

    # build circuit
    circ = Circuit()

    # initialize registers
    if nx+ny>0
        push!(circ, GateH(), vcat(xreg, yreg) )
    end

    append!(circ, quantumfouriertransform(zreg))
    for k in 1:K
        # begin oracle
        append!(circ, inv_phimult)
        append!(circ, diffusion(zreg, [xreg...,yreg...]) ) 
        append!(circ, phimult) 
        # end oracle     

        # diffusion operator
        append!(circ, diffusion([xreg...,yreg...], zreg) ) 
    end
    append!(circ, inverse(quantumfouriertransform(zreg)))
    return circ
end

"""
    get_factors(samples, foracle, registers)

Extracts the factors p and q from the sampled bitstrings

### Arguments
* occurances::Dict{BitString, Int64} Dictionary of sampled bitstrings and occurances
* foracle::NTuple{4,Int} Polynomial coefficients for oracle
* registers::NTuple(3,Vector{Int})

Returns a vector of [(p,q,occurances,bitstring), ...]
"""
function get_factors(occurances, foracle::NTuple{4,Int}, registers)
    xreg, yreg, zreg = registers

    occurances = sort(occurances, byvalue=true)

    factors = []
    for (bs_,o) in occurances

        if length(xreg)==0 a=0 else a=bitstring_to_integer(bs_[xreg]) end
        if length(yreg)==0 b=0 else b=bitstring_to_integer(bs_[yreg]) end

        if foracle[1]*a*b + foracle[2]*a + foracle[3]*b + foracle[4]==0 
            p = 6*a+foracle[3]
            q = 6*b+foracle[2]
            push!(factors, (p,q,o,bs_))
        end

    end

    return factors
end

;