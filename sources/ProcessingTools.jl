type graph
    J::Array{Float64,2}
    h::Array{Float64,1}
    L::Int64
    q::Int64
end

"""
   ReturnFreqs(Y::Array{Int64,2}, w::Array{Float64,1}, q::Int64)

Return single point and pairwise frequencies for input array `Y`. Weights have to be provided in vector `w`.  
"""
function ReturnFreqs(Y::Array{Int64,2}, w::Array{Float64,1}, q::Int64)
    Y = Y'
    (L,M) = size(Y)
    f2 = zeros(Float64, L*q, L*q)
    f1 = zeros(Float64,L*q) 

    if size(w)[1]!=M
        error("AlignmentStats.jl - ReturnFreqs: incorrect number of weights\n")
    end 

    Meff = sum(w)
    flush(STDOUT)
    for m in 1:M
        for j in 1:L
            f1[(j-1)*q+Y[j,m]] += w[m]
            f2[(j-1)*q+Y[j,m], (j-1)*q+Y[j,m]] += w[m]
            for i in (j+1):L
                f2[(i-1)*q+Y[i,m], (j-1)*q+Y[j,m]] += w[m]
                f2[(j-1)*q+Y[j,m], (i-1)*q+Y[i,m]] += w[m]
            end
        end
    end

    f2 = f2./Meff
    f1 = f1./Meff

    return (f1,f2)
end

"""
   ReturnFreqs(msa::Array{Int64,2}; q=findmax(msa)[1], reweighting=0, weights::String="")

Return single point and pairwise frequencies for input `msa`. 
"""
function ReturnFreqs(msa::Array{Int64,2}; q=findmax(msa)[1], reweighting=0, weights::String="", saveweights::String="", verbose::Bool=true)

    (M,L) = size(msa)
    f_2 = zeros(Float64,L*q, L*q)
    f_1 = zeros(Float64,L*q)
    w = ones(Float64,M)
    if weights!=""
        verbose?print("AlignmentStat.jl - ReturnFreqs: Reading weights... "):Void
        w = readdlm(weights)
        verbose?println(size(w)):Void
        if size(w)[1]!=M
            error("AlignmentStats.jl - ReturnFreqs: incorrect number of weights\n")
        end
        verbose?println("done."):Void
    elseif reweighting !=0
        verbose?println("AlignmentStat.jl - ReturnFreqs: Computing weights... "):Void
    verbose?begin @time w = ComputeWeights(msa,output=saveweights) end:begin w = ComputeWeights(msa,output=saveweights) end
        verbose?println("done."):Void
    end

    verbose?println("Computing frequencies..."):Void
    verbose?begin println("Meff = ",sum(vec(w))); @time (f_1,f_2) = ReturnFreqs(msa, vec(w), q) end:begin (f_1,f_2) = ReturnFreqs(msa, vec(w), q) end
    return (f_1,f_2)

end

"""
    FittingQuality(f2_1::Array{Float64,2}, f1_1::Array{Float64,1}, f2_2::Array{Float64,2}, f1_2::Array{Float64,1}, q::Int64; withdiag=false)
Compares frequencies `(f1_1,f2_1)` to `(f1_2,f2_2)`. Output is 
1. Pearson correlation between connected correlations. (default: without diagonal elements)
2. Frobenius norm of the difference between connected correlations.
3. Same as 1. for magnetizations.
4. Same as 2. for magnetizations.
"""
function FittingQuality(f2_1::Array{Float64,2}, f1_1::Array{Float64,1}, f2_2::Array{Float64,2}, f1_2::Array{Float64,1}, q::Int64; withdiag=false)
    L::Int64=size(f2_1)[1]/q
    C1 = f2_1 - f1_1*f1_1'
    C2 = f2_2 - f1_2*f1_2'
    id = zeros(Bool,L*q,L*q)
    # if !withdiag
    #     for i = 1:L
    #         C1[(i-1)*q+(1:q), (i-1)*q+(1:q)] = 0
    #         C2[(i-1)*q+(1:q), (i-1)*q+(1:q)] = 0
    #     end
    # end
    for i = 1:L
        for j = (i+1):L
            id[(j-1)*q+(1:q),(i-1)*q+(1:q)]=true
        end
    end
    cc = cor(C1[id],C2[id])
    froc = vecnorm(C1[id] - C2[id])
    cm = cor(f1_1,f1_2)
    from = vecnorm(f1_1-f1_2)
    return (cc,froc,cm,from)
end


"""
    ComputeWeights(msa::String ; theta::Float64 = 0.2)

Compute weights of each sequence in file `msa` using reweighting threshold `theta` (default 0.2).
"""
function ComputeWeights(msa::String; theta::Float64 = 0.2, output::String = "")
    Y = readdlm(msa,Int64);
    if findmin(Y[2:end,:])[1]==0
        Y = Y[2:end,:];
        Y = Y +1;
    end

    return ComputeWeights(Y,theta=theta,output=output)
end

function ComputeWeights(Y::Array{Int64,2};theta::Float64=0.2,output::String="")
    Y = Y';
    M::Int64 = size(Y,2)
    L::Int64 = size(Y,1)
    h::Float64 = L * (1-theta)
    weights = ones(Float64, M,1);
    d::Int64=0
    for m = 1:M
        for l = (m+1):M
            d = 0
            for i = 1:L
                d += Y[i,m]==Y[i,l];
            end
            if d > h
                weights[m]+=1
                weights[l]+=1
            end
        end
    end
    if output!=""
        writedlm(output, 1./weights, ' ')
    end
    return 1./weights
end


"""
    SymKLDiv(graph1::graph, graph2::graph)

KL-div between `graph1` and `graph2`, using `sample1` and `sample2` for averaging. Output syntax: E_a_b --> sample from a, energies from b.
"""

function SymKLDiv(graph1::graph, graph2::graph, sample1::Array{Int64,2}, sample2::Array{Int64,2},)
    M1 = size(sample1)[1]
    M2 = size(sample2)[1]
    energies_s1 = Array{Float64,1}(M1)
    energies_s2 = Array{Float64,1}(M2)
    # Sample 1 - energy in 1
    ComputeEnergies!(energies_s1, graph1.J, graph1.h, sample1, graph1.q)
    E_1_1 = mean(energies_s1)
    # Sample 1 - energy in 2
    ComputeEnergies!(energies_s1, graph2.J, graph2.h, sample1, graph1.q)
    E_1_2 = mean(energies_s1)
    # Sample 2 - energy in 1
    ComputeEnergies!(energies_s2, graph1.J, graph1.h, sample2, graph1.q)
    E_2_1 = mean(energies_s2)
    # Sample 2 - energy in 2
    ComputeEnergies!(energies_s2, graph2.J, graph2.h, sample2, graph1.q)
    E_2_2 = mean(energies_s2)

    dkl = (E_1_2 - E_1_1) + (E_2_1 - E_2_2)
    return(dkl, E_1_1, E_1_2, E_2_1, E_2_2)
end

"""
    ComputeEnergies!(energies::Array{Float64,1}, J::Array{Float64,2}, h::Array{Float64,1}, sample::Array{Int64,2})

Computes energies of all configurations in `sample` with graph `J` and `h`.
"""
function ComputeEnergies!(energies::Array{Float64,1}, J::Array{Float64,2}, h::Array{Float64,1}, sample::Array{Int64,2}, q::Int64)

    (M,L) = size(sample)
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
    ReadGraphStruct(parameters::String, L::Int64, q::Int64; format="mat",gauge=true)
Reads potts parameters in file `parameters`. `format` can be either `'mat"` or `"MCMC"`. Output is of custom type `graph`. 
"""
function ReadGraphStruct(parameters::String, L::Int64, q::Int64; format="mat",gauge=true)
    if format == "MCMC"
        (J,h) = ReadGraphMCMC(parameters, L, q)
    elseif format == "mat"
        (J,h) = ReadGraphMat(parameters, q)
    else
        error("ModelTools.jl - ReadGraph: Incorrect format string.\n")
    end
    if gauge
        SwitchGauge!(J,h,L,q)
    end
    return graph(J,h,L,q)
end

"""
    ReadGraph(parameters::String, L::Int64, q::Int64; format="mat",gauge=true)
"""
function ReadGraph(parameters::String, L::Int64, q::Int64; format="mat",gauge=true)
    if format == "MCMC"
        (J,h) = ReadGraphMCMC(parameters, L, q)
    elseif format == "mat"
        (J,h) = ReadGraphMat(parameters, q)
    else
        error("ModelTools.jl - ReadGraph: Incorrect format string.\n")
    end
    if gauge
        SwitchGauge!(J,h,L,q)
    end
    return (J,h)
end

"""
    ReadGraphMat(parameters::String, q::Int64)

Reads potts parameters in dlm format (ie stored as a matrix, with h being the last line).  
Outputs J and h. 
"""
function ReadGraphMat(parameters::String, q::Int64)
    f = open(parameters)
    J::Array{Float64,2} = readdlm(f,Float64)
    close(f)
    h::Array{Float64,1} = J[end,:]
    J = J[1:(end-1),:]
    return (J,h)
end

"""
    SwitchGauge!(J::Array{Float64,2}, h::Array{Float64,1}, L::Int64, q::Int64)
"""
function SwitchGauge!(J::Array{Float64,2}, h::Array{Float64,1}, L::Int64, q::Int64)
    t = zeros(Float64,q,q)
    u = zeros(Float64,q,q)

    for i in 1:L
        h[(i-1)*q+(1:q)] -=  mean(h[(i-1)*q+(1:q)])
        for j in 1:L
            t = view(J,(i-1)*q+(1:q), (j-1)*q+(1:q))
            h[(i-1)*q+(1:q)] += mean(t,2)
            h[(i-1)*q+(1:q)] .-= mean(t)
        end
    end

    for i in 1:L
        for j in (i+1):L
            u = view(J,(i-1)*q+(1:q), (j-1)*q+(1:q))
            J[(i-1)*q+(1:q), (j-1)*q+(1:q)] .-= mean(u,1)
            J[(i-1)*q+(1:q), (j-1)*q+(1:q)] .-= mean(u,2)
            J[(i-1)*q+(1:q), (j-1)*q+(1:q)] .+= mean(u)
        end
    end 
    for i in 1:L*q
        for j in (i+1):L*q
            J[j,i] = J[i,j]
        end
    end
    return ()
end

"""
    PseudoLikelihood(Y::Array{Int64,2}, J::Array{Float64,2}, h::Array{Float64,1}, q::Int64)

Y should be the transpose of the actual alignment. 
"""
function PseudoLikelihood(Y::Array{Int64,2}, w::Array{Float64,1}, J::Array{Float64,2}, h::Array{Float64,1}, q::Int64)
    (L,M) = size(Y)
    PL = 0.
    Z = 0.
    p = 0.
    E = 0.
    for m in 1:M
        for i in 1:L
            # Local partition function
            Z = 0.
            for a in 1:q
                E = 0.
                for j in 1:L
                    E += J[(j-1)*q+Y[j,m], (i-1)*q+a] + h[(j-1)*q+Y[j,m]]
                end
                Z += exp(E)
            end
            # Local likelihood of data
            E = 0.
            for j in 1:L
                E += J[(j-1)*q+Y[j,m], (i-1)*q+Y[i,m]] + h[(j-1)*q+Y[j,m]]
            end
            PL += log(exp(E)/Z) * w[m]
        end
    end
    PL = PL/sum(w)
    return PL
end