using JSON
data = JSON.parsefile("../fauci-email-graph.json")


# include("methods.jl")
emails = data["emails"]
names = data["names"]

## Collect empty subject errors
EmptySubject = Vector{Tuple{Int64,Int64}}()
Subject1 = Vector{Tuple{Int64,Int64}}()
for i = 1:length(emails)
    thread = emails[i]
    for j = 1:length(thread)
        email = thread[j] 
        sub = email["subject"]
        if length(sub) < 1
            push!(EmptySubject,(i,j))
        end
        if length(sub) == 1
            # There are two one-character emails. These exactly match the pdf, there is no error
            push!(Subject1,(i,j))
        end
    end
end

## Save the original way the email was processed in a text file.
# Save another text file with abbreviated information about sender/receiver/cc/date/subject (subject will be empty)
# Save two blank .txt files, one for the fixed version of the sender/receiver/cc/date/subject info, and one for the fixed email body text.
# Then manually go make fixes and put correct processing into the new files.

for k = 1:length(EmptySubject)

    thread = EmptySubject[k][1]
    emailno = EmptySubject[k][2]
    email = emails[thread][emailno]
    mkdir("nosubject-emails/$thread-$emailno")
    io = open("nosubject-emails/$thread-$emailno/original-allprint.txt", "w")
    fromname = names[email["sender"]+1]
    tonames = join(map(x-> names[x+1], email["recipients"]), "; ")
    ccnames = join(map(x-> names[x+1], email["cc"]), "; ")
    println(io,"From: $fromname")
    println(io,"To: $tonames")
    println(io,"CC: $ccnames")
    println(io,"Date: ", email["time"])
    println(io,"Subject: ", email["subject"])
    println(io,"Body: \n ")
    println(io,email["body"])
    close(io)

    io = open("nosubject-emails/$thread-$emailno/original-short.txt", "w")
    println(io,email["sender"])
    print(io,"")
    for j = 1:length(email["recipients"])-1
        J = email["recipients"][j]
        print(io,"$J,")
    end
    print(io,"$(email["recipients"][end])")
    println(io)
    for j = 1:length(email["cc"])-1
        J = email["cc"][j]
        print(io,"$J,")
    end
    if length(email["cc"]) > 0
        print(io,"$(email["cc"][end])")
    end
    println(io)
    println(io,email["time"])
    println(io,email["subject"])
    close(io)

    io = open("nosubject-emails/$thread-$emailno/fixed-short.txt", "w")
    io = open("nosubject-emails/$thread-$emailno/fixed-body.txt", "w")

end

