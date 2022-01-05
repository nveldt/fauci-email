using JSON
data = JSON.parsefile("fauci-email-data-fixed.json")
emails = data["emails"]

dirs = readdir("nosubject-emails-done")
for k = 1:length(dirs)
    if ~isdir("nosubject-emails-done/$(dirs[k])")
        continue
    end
    threademail = parse.(Int64,split(dirs[k],"-")) 
    thread = threademail[1]
    emailno = threademail[2]
    email = emails[thread][emailno]
    io = open("nosubject-emails-done/$thread-$emailno/check-allprint.txt", "w")
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

    io = open("nosubject-emails-done/$thread-$emailno/check-short.txt", "w")
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
end


## Same checks for the longsubject errors
dirs = readdir("longsubjecterror-emails-done")
for k = 1:length(dirs)
    if ~isdir("longsubjecterror-emails-done/$(dirs[k])")
        continue
    end
    threademail = parse.(Int64,split(dirs[k],"-")) 
    thread = threademail[1]
    emailno = threademail[2]
    email = emails[thread][emailno]
    io = open("longsubjecterror-emails-done/$thread-$emailno/check-allprint.txt", "w")
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

    io = open("longsubjecterror-emails-done/$thread-$emailno/check-short.txt", "w")
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
end