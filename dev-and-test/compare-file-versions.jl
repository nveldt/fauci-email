## Test email-2.json against email.json...
include("../methods.jl")
##
data1 = JSON.parsefile("fauci-email-graph.json")
data2 = JSON.parsefile("fauci-email-graph-2.json")
##
for tid in 1:length(data1["emails"])
  t1 = data1["emails"][tid]
  t2 = data2["emails"][tid]
  if length(t1) != length(t2)
    println("Threads differ in length")
  end
end
##
println("Testing difference in recipients")
mydiff = false
for tid in 1:length(data1["emails"])
  t1 = data1["emails"][tid]
  t2 = data2["emails"][tid]
  for eid in 1:length(t1)
    if t1[eid]["recipients"] != t2[eid]["recipients"]
      println("Difference in recipients")
      mydiff = true
      break
    end
  end
end
@assert(mydiff==false)
##
println("Testing difference in body")
mydiff = false
for tid in 1:length(data1["emails"])
  t1 = data1["emails"][tid]
  t2 = data2["emails"][tid]
  for eid in 1:length(t1)
    if t1[eid]["text"] != t2[eid]["body"]
      println("Difference in body")
      mydiff = true
      break
    end
  end
end
@assert(mydiff==false)
##
println("Testing difference in time")
mydiff = false
for tid in 1:length(data1["emails"])
  t1 = data1["emails"][tid]
  t2 = data2["emails"][tid]
  for eid in 1:length(t1)
    if t1[eid]["time"] != t2[eid]["time"]
      println("Difference in time")
      mydiff = true
      break
    end
  end
end
@assert(mydiff==false)
##
println("Testing difference in cc")
mydiff = false
for tid in 1:length(data1["emails"])
  t1 = data1["emails"][tid]
  t2 = data2["emails"][tid]
  for eid in 1:length(t1)
    if t1[eid]["cc"] != t2[eid]["cc"]
      println("Difference in cc")
      mydiff = true
      break
    end
  end
end
@assert(mydiff==false)
##
println("Testing difference in subject")
mydiff = false
for tid in 1:length(data1["emails"])
  t1 = data1["emails"][tid]
  t2 = data2["emails"][tid]
  for eid in 1:length(t1)
    if t1[eid]["subject"] != t2[eid]["subject"]
      println("Difference in subject")
      mydiff = true
      break
    end
  end
end
@assert(mydiff==false)
##
