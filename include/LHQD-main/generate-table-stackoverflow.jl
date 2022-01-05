using Statistics
using JLD,FileIO

include("common.jl")
include("local-hyper.jl")

method_map = Dict()
method_map["LH-2.0"] = "LH-2.0"
method_map["LH-1.4"] = "LH-1.4"
method_map["OneHop+flow"] = "OneHop+flow"
method_map["weighted_clique+acl"] = "WCE+ACL"
method_map["HGCRD"] = "HGCRD"
method_map["clique+acl"] = "UCE+ACL"
method_map["acl"] = "star+ACL"
method_map["LH-2.0+flow"] = "LH-2.0+flow"

method_ordering = ["LH-2.0","LH-1.4","LH-2.0+flow","acl","weighted_clique+acl","clique+acl","OneHop+flow","HGCRD"]
H,clusters = read_dataset("stackoverflow")
records = load("results/stackoverflow_records.jld")
wptr = open("stackoverflow.tex","w")
fields = ["Alg","runtime","PR","RC","F1","Cond"]
format = join(["l" for i = 1:(1+length(fields))],"")

write(wptr,"""\\noindent \\begin{tabularx}{0.75\\linewidth}{$format}\n""")
write(wptr,"""\\toprule\n""")
for (i,field) in enumerate(fields)
    if i != length(fields)
        write(wptr,"""$field &""")
    else
        write(wptr,"""$field \\\\\n""")
    end
end
write(wptr,"""\\midrule\n""")

stats = Dict()
stats["runtime"] = zeros(length(method_ordering),length(keys(clusters)))
for (i,method) in enumerate(method_ordering)
    for (j,label) in enumerate(collect(keys(clusters))[6:end])
        stats["runtime"][i,j] = mean([t[end] for t in records[string(label)][method]])
    end
end

stats["PR"] = zeros(length(method_ordering),length(keys(clusters)))
for (i,method) in enumerate(method_ordering)
    for (j,label) in enumerate(collect(keys(clusters))[6:end])
        stats["PR"][i,j] = median([t[end-3] for t in records[string(label)][method]])
    end
end

stats["RC"] = zeros(length(method_ordering),length(keys(clusters)))
for (i,method) in enumerate(method_ordering)
    for (j,label) in enumerate(collect(keys(clusters))[6:end])
        stats["RC"][i,j] = median([t[end-2] for t in records[string(label)][method]])
    end
end

stats["Cond"] = zeros(length(method_ordering),length(keys(clusters)))
for (i,method) in enumerate(method_ordering)
    for (j,label) in enumerate(collect(keys(clusters))[6:end])
        stats["Cond"][i,j] = median([t[2] for t in records[string(label)][method]])
    end
end

for (k,method) in enumerate(method_ordering)
    name = method_map[method]
    curr_row = """$name &"""
    for (i,field) in enumerate(fields)
        if field == "Alg"
            continue
        elseif field == "runtime"
            t = round(median(stats[field][k,:]),digits=3)
            curr_row *= """$t &"""
        elseif field == "PR"
            t = round(median(stats["PR"][k,:]),digits=2)
            curr_row *= """$t &"""
        elseif field == "RC"
            t = round(median(stats["RC"][k,:]),digits=2)
            curr_row *= """$t &"""
        elseif field == "Cond"
            t = round(median(stats["Cond"][k,:]),digits=2)
            curr_row *= """$t &"""
        elseif field == "F1"
            t1 = round(median(stats["PR"][k,:]),digits=2)
            t2 = round(median(stats["RC"][k,:]),digits=2)
            t = round(t1*t2*2/(t1+t2),digits=2)
            curr_row *= """$t &"""
        else
            _,pos = findmax(stats["F1"],dims=1)
            ranks = collections.Counter(vec([t[1] for t in pos]))
            @show ranks
            if haskey(ranks,k) t = ranks[k] else t=0 end
            curr_row *= """$t"""
        end
    end
    write(wptr,curr_row*"""\\\\ \n""")
end

write(wptr,"""\\bottomrule \n""")
write(wptr,"""\\end{tabularx} \n\n\n\n""")


close(wptr)