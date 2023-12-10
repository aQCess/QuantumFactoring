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


""" returns the number of bits in argument j """
function nbits(j)
    if j == 0 return 0 end
    return length(digits(j, base=2))
end

"""
Return registers with a fixed size without knowledge of the factors

### Arguments
* N::Integer Number to be factorized N = p*q
* d::Integer Bit length difference (ny-nx)/2
"""
function registers(N::Integer, d)
    k = nbits(N)-4  # force 2n-5
    nx,ny = floor(Int, k/2-d), ceil(Int, k/2+d)
    nz = nx + ny + 3

    if nx < 0 || ny < 0
        return nothing
    end

    # create registers
    xreg = Int[1:nx...] 
    yreg = Int[nx+1:nx+ny...]
    zreg = Int[nx+ny+1:nx+ny+nz...]

    return (xreg, yreg, zreg)
end


"""
Return registers with the optimal (minimum) size. 
Requires knowledge of the factors p,q

### Arguments
* pq::Tuple(Int,Int) prime factors
"""
function registers(pq::Tuple{Integer,Integer})
    p,q = pq[1],pq[2]

    a = ifelse(p % 6 == 1, div(p-1,6)-1, div(p+1,6)-1)
    b = ifelse(q % 6 == 1, div(q-1,6)-1, div(q+1,6)-1)
    s = Int(1.5-0.5*mod(p,6))

    nx, ny = nbits(a), nbits(b)
    nz = nx + ny + 3

    xreg = Int[1:nx...] 
    yreg = Int[nx+1:nx+ny...]
    zreg = Int[nx+ny+1:nx+ny+nz...]
    return (xreg, yreg, zreg), s
end

nothing
