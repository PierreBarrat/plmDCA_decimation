include("ProcessingTools.jl")

# fam = "PF00014_raw"
fam = ARGS[1]
cd(@sprintf("%s_results",fam))
N = parse(Int64,readstring(pipeline(`head -n 1 msa_numerical.txt`,`awk '{print NF}'`)))
q = 21

# List of samples
list_samples =  readstring(pipeline(`find . -regextype sed -regex '.*sample_\([0-9]\|[0-9][0-9]\).txt'`,`sort -t"_" -k2n`))
list_samples = map(String, split(list_samples, "\n")[1:end-1]) 

# Finding different decimation steps for samples
ld = map(x->parse(Int64,String(match(r"[0-9]*[0-9]",x).match)),list_samples)
K = size(ld)[1]

# List of inferred graphs which have a sample
list_graph_t = readstring(pipeline(`find . -regextype sed -regex '.*plmInf_\([0-9]\|[0-9][0-9]\)_mat.txt'`,`sort -t"_" -k3n`))
list_graph_t = map(String, split(list_graph_t, "\n")[1:end-1]) 
ld_graph = map(x->parse(Int64,String(match(r"[0-9]*[0-9]",x).match)),list_graph_t)
Kg = size(ld_graph)[1]
list_graph = Array{String,1}(0)
for i = 1:Kg
	if in(ld_graph[i],ld)
		push!(list_graph, list_graph_t[i])
	end
end


sample_nat = readdlm("msa_numerical.txt", Int64)
(f1_nat,f2_nat) = ReturnFreqs(sample_nat, weights="weights.txt",verbose=false)
out_kl = zeros(Float64,K,K)
out_fitqual = zeros(K,5)
for i = 1:K
	println("i = $i/$K")
	begin
		# Fitting quality
		sample_i = readdlm(list_samples[i],Int64,skipstart=1)+1
		(f1,f2) = ReturnFreqs(sample_i,verbose=false)
		out_fitqual[i,1] = ld[i]
		(out_fitqual[i,2],out_fitqual[i,3],out_fitqual[i,4],out_fitqual[i,5]) = FittingQuality(f2_nat, f1_nat, f2, f1, q)
		# Pairwise KL-div 
		graph_i = ReadGraphStruct(list_graph[i],N,q)
		for j = (i+1):K
			sample_j = readdlm(list_samples[j],Int64,skipstart=1)+1
			graph_j = ReadGraphStruct(list_graph[j],N,q)
			out_kl[i,j] = SymKLDiv(graph_i, graph_j, sample_i, sample_j)[1]
			out_kl[j,i] = out_kl[i,j]
		end
	end
end

run(`mkdir -p processing_results`)
save_ld = @sprintf("processing_results/%s_samplesteps.txt",fam)
save_kl = @sprintf("processing_results/%s_KLdiv.txt",fam)
save_fitqual = @sprintf("processing_results/%s_fittingquality.txt",fam)

# Saving list of decimation indexes
writedlm(save_ld,ld)
# Saving KLdiv
writedlm(save_kl, out_kl)
# Saving fitting quality
f = open(save_fitqual, "w")
write(f, "Decimation_index CC_pears CC_fronorm Mag_pears Mag_fronorm\n")
writedlm(f,out_fitqual)
close(f)

cd("../")
