include("methods.jl")
include("methods_temporal.jl")
##
println("Threads: $(length(data["emails"]))")
println("Emails: $(sum(length.(data["emails"])))")
##
println("Duplicates: $(_emails_by_day_and_time(data).duplicates)")
