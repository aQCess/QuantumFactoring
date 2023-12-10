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

struct FactorReport
    data::Dict
    FactorReport(res...) = new(Dict(res))
end
Base.getindex(r::FactorReport,j) = if haskey(r.data,j) r.data[j] else nothing end

function Base.show(io::IO, r::FactorReport)
    if isempty(r.data)
        println(io, "Could not find prime factors. Possibly a prime.")
    else
        p,q,P = r["p"], r["q"], r["P"]
        if r["alg"] == "grover"
            println(io, "Result :\t\t\t $p x $q = $(p*q)")
            println(io, "Trials :\t\t\t $(r["trials"])")
            println(io, "Success probability :\t\t $P ($(r["occurances"])/$(r["nsamples"]) samples)")
            println(io, "Simulation time :\t\t $(r["tapply"]) seconds")
            println(io, "Size of the Hilbert space :\t 2^$(r["nqubits"])")
            println(io, "Grover iteration steps :\t $(r["steps"])")
        elseif r["alg"] == "N%3" || r["alg"] == "N%2"
            println(io, "Result :\t\t\t $p x $q = $(p*q)")
        end
        if r["compiled"]==false
            println(io, "Number of one-qubit gates :\t $(r["n1gates"])")
            println(io, "Number of two-qubit gates :\t $(r["n2gates"])")
        end
    end
end

function Base.show(io::IO, ::MIME"text/html", r::FactorReport)
    if isempty(r.data)
        println(io, "Could not find prime factors. Possibly a prime.")
    else
        N,p,q,P = r["N"], r["p"], r["q"], r["P"]    
        print(io, "<font size=6><center>Factoring result: $(p*q) = $p x $q </center></font>")
        if r["alg"] == "grover"
            print(io, "<center>With success probability $P ($(r["occurances"])/$(r["nsamples"]) samples)</center>")

            print(io, "<table>")
            print(io, "<tr><td> <b>Trials :</b></td><td><b>$(r["trials"]) </b> </td><tr>")
            print(io, "<tr><td> <b>Simulation time :</b></td><td><b>$(r["tapply"]) seconds </b> </td><tr>")
            print(io, "<tr><td> Size of the Hilbert space :</td><td>2^$(r["nqubits"])</td><tr>")
            print(io, "<tr><td> Grover iteration steps :</td><td>$(r["steps"])</td><tr>")
        elseif r["alg"] == "N%3" || r["alg"] == "N%2"
            println(io, "<center>Trivial factor found</center>")
        end


        if r["compiled"]==false
            print(io, "<tr><td> <b>Number of one-qubit gates :</b></td><td><b>$(r["n1gates"])</b></td><tr>")
            print(io, "<tr><td> <b>Number of two-qubit gates :</b></td><td><b>$(r["n2gates"])</b></td><tr>")
        end
        print(io, "</table>")
    end
end

nothing