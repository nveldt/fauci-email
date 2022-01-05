using JLD,FileIO,Statistics

method_map = Dict()
method_map["LH-2.0"] = "LH-2.0"
method_map["LH-1.4"] = "LH-1.4"
method_map["OneHop+flow"] = "OneHop+flow"
method_map["acl"] = "star+ACL"
method_map["HGCRD"] = "HGCRD"
method_map["clique+acl"] = "UCE+ACL"
method_map["weighted_clique+acl"] = "WCE+ACL"
method_map["LH-2.0+flow"] = "LH-2.0+flow"

method_ordering = ["LH-2.0","LH-1.4","LH-2.0+flow","acl","weighted_clique+acl","clique+acl","OneHop+flow-005","HGCRD"]
labels = [12,18,17,25,15,24]

wptr = open("bulk_eval_table_amazon.tex","w")

write(wptr,"""\\noindent \\begin{tabularx}{\\linewidth}{@{\\,\\,}l@{\\,\\,}*{6}{@{}l@{}X@{}}@{}}\n""")
write(wptr,"""\\toprule\n""")
write(wptr,""" Alg & \n""")
write(wptr,"""\\mbox{\\rlap{12}} & &\n""")
write(wptr,"""\\mbox{\\rlap{18}} & &\n""")
write(wptr,"""\\mbox{\\rlap{17}} & &\n""")
write(wptr,"""\\mbox{\\rlap{25}} & &\n""")
write(wptr,"""\\mbox{\\rlap{15}} & &\n""")
write(wptr,"""\\mbox{\\rlap{24}} &  \\\\\n""")
write(wptr,""" & \\mbox{\\rlap{\\footnotesize F1 \\& Med.}}  & &\n""")
write(wptr,"""\\mbox{\\rlap{\\footnotesize F1 \\& Med.}} & &\n""")
write(wptr,"""\\mbox{\\rlap{\\footnotesize F1 \\& Med.}} & &\n""")
write(wptr,"""\\mbox{\\rlap{\\footnotesize F1 \\& Med.}} & &\n""")
write(wptr,"""\\mbox{\\rlap{\\footnotesize F1 \\& Med.}} & &\n""")
write(wptr,"""\\mbox{\\rlap{\\footnotesize F1 \\& Med.}} &    \\\\\n""")
write(wptr,"""\\midrule\n""")

function entry_function(median_val,img)
    # offset the median text entry
    width=20
    include = "\\includegraphics[width=$(width)pt,height=10pt]{$img}"
    img = "\\mbox{\\rlap{\\raisebox{-2pt}{$include}}}"
    offset = round(Int,median_val*15)
    text = "\\mbox{\\rlap{\\mbox{\\hspace{$(offset)pt}$median_val}}}"
    row = "& $img$text & "
    return row
end

records = load("results/amazon_records.jld")
for (i,method) in enumerate(method_ordering)
    curr_row = """"""
    name = method_map[method]
    curr_row *= """ $name """
    for label in labels
        @show label
        f1 = []
        for tmp in records[string(label)][method]
            push!(f1,tmp[end-1])
        end
        median_val = round(median(f1),digits=2)
        figmethod = replace(method,"."=>"")
        curr_row *= entry_function(median_val,"f1-$file-$label-$figmethod")
    end
    write(wptr,curr_row*"""\\\\ \n""")
end

write(wptr,"""\\bottomrule \n""")
write(wptr,"""\\end{tabularx} \n\n\n\n""")


close(wptr)


wptr = open("runtime_amazon.tex","w")
format = join(["l" for i = 1:(1+length(labels))],"")

write(wptr,"""\\noindent \\begin{tabularx}{0.95\\linewidth}{$format}\n""")
write(wptr,"""\\toprule\n""")
write(wptr,""" Alg & """)
for (i,label) in enumerate(labels)
    if i != length(labels)
        write(wptr,"""$label &""")
    else
        write(wptr,"""$label \\\\\n""")
    end
end
write(wptr,"""\\midrule\n""")

for method in method_ordering
    name = method_map[method]
    curr_row = """$name &"""
    for (i,label) in enumerate(labels)
        record = records[string(label)][method]
        runtime = round(median([t[end] for t in record]),digits=1)
        if i == length(labels)
            curr_row *= """$runtime"""
        else
            curr_row *= """$runtime &"""
        end
    end
    write(wptr,curr_row*"""\\\\ \n""")
end

write(wptr,"""\\bottomrule \n""")
write(wptr,"""\\end{tabularx} \n\n\n\n""")

close(wptr)