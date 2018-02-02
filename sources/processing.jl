# function SymKLDiv(J1::Array{Float64,2}, h1::Array{Float64,1}, sample1::Array{Int64,1}, J2::Array{Float64,2}, h2::Array{Float64,1}, sample2::Array{Float64,2})



# end

type graph
    J::Array{Float64,2}
    h::Array{Float64,1}
    sample::Array{Int64,2}    
    L::Int64
end

function ReadGraphs(path::String, steps::Array{Int64,1}, L::Int64)
    q = 21
    S = size(steps)
    graph_list = Array{graph,1}(S)
    pos = 1
    @time for s in steps
        @printf("Reading graph %d ... ", s)
        J = readdlm(@sprintf("%s/decimation_results/plmInf_%d_mat.txt", path, s))
        h = J[end,:]
        J = J[1:(end-1),:]
        sample = readdlm(@sprintf("%s/samples/sample_%d.txt", path, s))
        sample = convert(Array{Int64,2}, sample[2:end,:]+1)
        graph_list[pos] = graph(J,h,sample,L)
        pos += 1
        @printf("done.\n")
    end
    return graph_list
end

function SymKLDiv(graph_list::Array{graph,1})

    S = size(graph_list)[1]
    M = size(graph_list[1].sample)[1]
    pairw_kl =  zeros(Float64,S,S)
    energies = Array{Float64,1}(M)
    for s in 1:S
        @printf("KL div for model %d ...\n ", s)
        # Sample s - energy in s
        ComputeEnergies!(energies, graph_list[s].J, graph_list[s].h, graph_list[s].sample)
        E_s_s = mean(energies)
        @time for t in (s+1):S
            println(t)
            # Sample t - energy in s
            ComputeEnergies!(energies, graph_list[s].J, graph_list[s].h, graph_list[t].sample)
            E_t_s = mean(energies)
            # Sample s - energy in t
            ComputeEnergies!(energies, graph_list[t].J, graph_list[t].h, graph_list[s].sample)
            E_s_t = mean(energies)
            # Sample t - energy in t -- Could be more efficient with this one, ... 
            ComputeEnergies!(energies, graph_list[t].J, graph_list[t].h, graph_list[t].sample)
            E_t_t = mean(energies)
            pairw_kl[s,t] = (E_t_s - E_t_t) + (E_s_t - E_s_s)
            pairw_kl[t,s] = pairw_kl[s,t]
        end
    end
    return pairw_kl
end

function CompareCorrelations(graph_list::Array{graph,1}, C_nat::Array{Float64,2})

    S = size(graph_list)[1]
    q = 21
    L = convert(Int64,size(C_nat)[1]/q) 
    M = size(graph_list[1].sample)[1]
    weights = ones(Int64,M)
    results = zeros(Float64,S,2)
    for s in 1:S
        (f1,f2) = ReturnFreqs(graph_list[s].sample)
        C = f2 - f1*f1'
        results[s,1] = cor(reshape(C_nat, L*q*L*q,1), reshape(C, L*q*L*q,1))[1]
        results[s,2] = vecnorm(C_nat-C)
    end
    return results

end

"""
ComputeEnergies!(energies::Array{Float64,1}, J::Array{Float64,2}, h::Array{Float64,1}, sample::Array{Int64,2})

    Computes energies of all configurations in `sample` with graph `J` and `h`.
"""
function ComputeEnergies!(energies::Array{Float64,1}, J::Array{Float64,2}, h::Array{Float64,1}, sample::Array{Int64,2})

    (M,L) = size(sample)
    q = 21
    for m = 1:M
        energies[m] = 0
        for i = 1:L
            for j = (i+1):L
                energies[m] -= J[(i-1)*q+sample[m,i], (j-1)*q+sample[m,j]]
            end
        energies[m] -= h[(i-1)*q+sample[m,i]]
        end
    end
    return energies
end

"""
ComputeEnergies(J::Array{Float64,2}, h::Array{Float64,1}, sample::Array{Int64,2})

    Computes energies of all configurations in `sample` with graph `J` and `h`.
"""
function ComputeEnergies(J::Array{Float64,2}, h::Array{Float64,1}, sample::Array{Int64,2})

    (M,L) = size(sample)
    q = 21
    energies = zeros(Float64, M)
    for m = 1:M
        for i = 1:L
            for j = (i+1):L
                energies[m] -= J[(i-1)*q+sample[m,i], (j-1)*q+sample[m,j]]
            end
        energies[m] -= h[(i-1)*q+sample[m,i]]
        end
    end
    return energies
end



"""
   ReturnFreqs(msa::Array{Int64,2})

Return single point and pairwise frequencies for input `msa`. 
"""
function ReturnFreqs(msa::Array{Int64,2} ; weights=ones(Int64,size(msa)[1]))

    (M,L) = size(msa)
    q = 21
    f_2 = zeros(Float64,L*q, L*q)
    f_1 = zeros(Float64,L*q)
    for m in 1:M
        for j in 1:L
            f_1[(j-1)*q+msa[m,j]] += weights[m]
            f_2[(j-1)*q+msa[m,j], (j-1)*q+msa[m,j]] += weights[m]
            for i in (j+1):L
                f_2[(i-1)*q+msa[m,i], (j-1)*q+msa[m,j]] += weights[m]
                f_2[(j-1)*q+msa[m,j], (i-1)*q+msa[m,i]] += weights[m]
            end
        end
    end

    Meff = sum(weights)
    @printf("Meff : %f\n",Meff)
    f_2 = f_2./Meff
    f_1 = f_1./Meff

    return (f_1,f_2)

end