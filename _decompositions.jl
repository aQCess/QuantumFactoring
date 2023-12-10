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
    phimultiply(a, b, c, d, regs...)

Multiply and add operation in Fourier domain for z-> a*xy + b*x + c*y + d

one qubit gates: nz
two qubit gates: 5*nx*ny*(2*nz-nx-ny+2)/2 + nx*(2*nz-nx+1)/2 + ny*(2*nz-ny+1)/2 ~ n^3

#  n2gates < nz * (5*nx*ny + nx + ny) + 6*nx*ny - (5*nx*ny+nx+ny-1)*(nx+ny)/2
#       
"""
function phimultiply(a, b, c, d, regs...)

    # size of the registers
    nx, ny, nz = length.(regs)
    circ = Circuit()

    # multiplication a*x*y 
    for j in 0:nx-1, i in 0:ny-1, k in i+j:nz-1

        angle = a * 2*pi / 2.0^(k - j - i + 1)

        if mod(angle, 2*pi) != 0 # do nothing if lambda is multiple of 2pi

            ctr1, ctr2, target = regs[1][j+1], regs[2][i+1], regs[3][k+1]

            #push!(circ, Control(2,GateP(lambda)), ctr1, ctr2, target)
            # decomposed in terms of elementary two-qubit gates
            push!(circ, GateCP(angle/2), ctr2, target)
            push!(circ, GateCX(), ctr1, ctr2)
            push!(circ, GateCP(-angle/2), ctr2, target)
            push!(circ, GateCX(), ctr1, ctr2)
            push!(circ, GateCP(angle/2), ctr1, target)
        end

    end

    # add b*x
    for j in 0:nx-1, k in j:nz-1
        angle = b* 2*pi / 2.0^(k - j + 1)
        if mod(angle, 2*pi) != 0
            push!(circ, GateCP(angle), regs[1][j+1], regs[3][k+1])
        end
    end

    # add c*y
    for i in 0:ny-1, k in i:nz-1
        angle = c* 2*pi / 2.0^(k - i + 1)
        if mod(angle, 2*pi) != 0
            push!(circ, GateCP(angle), regs[2][i+1], regs[3][k+1])
        end
    end

    # add d
    for k in 0:nz-1
        angle = d * 2*pi / 2.0^(k+1)
        if mod(angle, 2*pi) != 0
            push!(circ, GateP(angle), regs[3][k+1])
        end
    end
    
    return circ
end


""" 
    Quantum Fourier Transform 

one qubit gates: n
two qubit gates: n(n-1)/2
"""
function quantumfouriertransform(reg)

    qc = Circuit()
    q1,q2 = 1,length(reg)
    
     # Generate multiple groups of diminishing angle CRZs and H gate
     for i_qubit in q1:q2
    
         hidx = q2 - i_qubit
    
         # multiple controlled RZs with decreasing angles
          if hidx < q2
              ncrz = i_qubit-1
              for j in q1:ncrz
                  divisor = 1<<(ncrz-j+1)
                  push!(qc, GateCRZ(pi / divisor), reg[hidx+q1], reg[q1+q2-j] )
              end
          end
    
         # H gates applied to all qubits
         push!(qc, GateH(), reg[hidx+q1])
      end
    return qc

end



	
""" 
    multiCX(qreg, qfree)
    
Uses the recursive decomposition described in [Barenco1995], Lemma 7.2

Assuming k control bits and 1 target bit. 
Requires an auxiliary register of k-4 free bits

one qubit gates: 16*n - 48
two qubit gates: 12*n - 36
"""
function multiCX(qreg, qfree)

	"""decomposed circuit for CCX gate modulo pi phase shifts 6.2"""
	function CCX(c1::Int, c2::Int, t::Int)
		cc = Circuit()
		push!(cc, GateRY(pi/4), t)
		push!(cc, GateCX(),  c2, t)
		push!(cc, GateRY(pi/4), t)
		push!(cc, GateCX(),  c1, t)
		push!(cc, GateRY(-pi/4), t)
		push!(cc, GateCX(),  c2, t)
		push!(cc, GateRY(-pi/4), t)
		return cc
	end
	
	c = Circuit()
	controls = qreg[1:end-1]
	target = qreg[end]

	if length(controls)==1
		c= push!(Circuit(), GateX(), controls[1])
	elseif length(controls)==2
		c= CCX(controls[1],controls[2], target)
	else
		zig = Circuit()
		#iterate over control bits, excluding first 2
		j=1
		for ctr1 in reverse(controls[3:end-1])
		    tgt, ctr2 = qfree[j], qfree[j+1]
		    append!(zig, CCX(ctr1, ctr2, tgt))
		    j+=1
		end
		
		zag = inverse(zig)

		append!(c, CCX(controls[end], qfree[1], target))
		append!(c, zig)
		append!(c, CCX(controls[1], controls[2], qfree[j]))
		append!(c, zag)
		append!(c, CCX(controls[end], qfree[1], target))
		append!(c, zig)
		append!(c, CCX(controls[1], controls[2], qfree[j]))
		append!(c, zag)
	end
	return c	
end

"""
one qubit gates: 16n^2 - 112n + 200
two qubit gates: 12n^2 - 82n + 149
"""
function multiCP(a, qreg, qfree)

    c = Circuit()
    nq = length(qreg)
    if nq == 1
        push!(c, GateP(a), qreg...)
    elseif nq == 2
        push!(c, GateCP(a), qreg...)
    elseif nq == 3
        push!(c, GateCP(a/2), qreg[2], qreg[3])
        push!(c, GateCX(), qreg[1], qreg[2])
        push!(c, GateCP(-a/2), qreg[2], qreg[3])
        push!(c, GateCX(), qreg[1], qreg[2])
        push!(c, GateCP(a/2), qreg[1], qreg[3])
    else
        push!(c, GateCP(a/2), qreg[end-1], qreg[end])      
        
        append!(c, multiCX(qreg[1:end-1], qfree))
        
        push!(c, GateCP(-a/2), qreg[end-1], qreg[end])
        
        append!(c, inverse(multiCX(qreg[1:end-1], qfree)))
        
        Q = [qreg[1:end-2]..., qreg[end]]
        append!(c, multiCP(a/2, Q, qfree))
    end

    return c
end


""" diffusion operator from Grovers algorithm 

uses multiCP gate using auxiliary qubits

one qubit gates: 16n^2 - 110n + 200
two qubit gates: 12n^2 - 82n + 149
"""
function diffusion(qreg, qfree)
    c = Circuit()
    if length(qreg)>1
        for q in qreg
            push!(c, GateRY(pi/2), q)
        end
        
        # multicontrol Z gate = C^k P(pi)
        append!(c, multiCP(pi, qreg, qfree))

        for q in qreg
            push!(c, GateRY(-pi/2), q)
        end
    end
    return c
end

nothing
