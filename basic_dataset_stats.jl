using LinearAlgebra
using MatrixNetworks
using JSON
using DelimitedFiles
using Printf
data = JSON.parsefile("fauci-email-data.json")

## First email for display in figure
email1 = data["emails"][1][1]

# Save how many emails each person sends, receives, and is CC'ed on.
Names = data["names"]
n = length(Names)

T = zeros(Int64,n,3) # sent, received, cc'ed

emails = data["emails"]

for thread in emails
    for email in thread
        sender = email["sender"] 
        rec = email["recipients"] 
        cced = email["cc"] 
        T[sender+1,1] += 1
        for k in rec
            T[k+1,2] += 1
        end
        for k in cced
            T[k+1,3] +=1
        end
    end
end



## Tables

p1 = sortperm(T[:,1],rev = true)
p2 = sortperm(T[:,2],rev = true)
p3 = sortperm(T[:,3],rev = true)

P = [p1 p2 p3]
k  = 10

for ind = 1:3
for i = 1:k
    j = P[i,ind]
    println("$(Names[j]) $(T[j,ind])")
end
end


## Print tables
keyset = ["Sender", "Receiver", "CC"]

for i = 1:3
    println("%")
    println("% -- ", keyset[i])
    println("%")
    println("\\begin{tabular}{*{3}{p{16pt}@{}}p{112pt}}")
    println("\\toprule")
    println("\\multicolumn{4}{c}{$(keyset[i])} \\\\")
    println("\\midrule")
    for j = 1:k 
      id = P[j,i] # index of the jth top person in the ith category
      nm = Names[id]
      val = T[id,i] #P[j,i]  # number of emails they are in with this role
      for w = 1:3
        if w == i
          print("\\textcolor{LightGray}{")
        end
        # print("$(T[id,w])")
        rankin = findfirst(x->x==id,P[:,w])
        valw = T[id,w]
        if valw == 0
          print("--")
        else
          print("$rankin")
        end
        if w == i
          print("}")
        end
        print(" & ")
      end
      println("$(Names[id]), $val \\\\")
    end
    println("\\bottomrule")
    println("\\end{tabular}")
end

