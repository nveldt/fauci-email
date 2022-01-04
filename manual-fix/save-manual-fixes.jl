## Here we will load in the manual fixes for the no-subject emails with errors

using JSON
data = JSON.parsefile("fauci-email-initial.json")

emails = data["emails"]
names = data["names"]
clusters = data["clusters"]


## First of all, there were a few new names we had to add for the no-subject errors
push!(names, "hicks, lauri")
push!(names, "roohi, shahrokh")
push!(names, "thompson, robert")
push!(names, "poznasky, mark")
push!(names, "richardson, juliana")
push!(clusters, 3)
push!(clusters, 3)
push!(clusters, 7)
push!(clusters, 7)
push!(clusters, 2)

## Secondly, we added one more name because of the long-subject errors
push!(names, "goldstein, sandra")
push!(clusters,7)

## Get all the subdirectories in the nosubject-emails directory

# use the tag "done" for the already manually processed fixes
dirs = readdir("longsubjecterror-emails-done")


## Go through and make corrections to the problematic emails
for k = 1:length(dirs)

    if ~isdir("longsubjecterror-emails-done/$(dirs[k])")
        continue
    end
    threademail = parse.(Int64,split(dirs[k],"-")) 
    shortinfo = readlines("longsubjecterror-emails-done/$(dirs[k])/fixed-short.txt")
    if shortinfo[1] == "none" || shortinfo[1] == "None"
        continue
    end
    body = read("longsubjecterror-emails-done/$(dirs[k])/fixed-body.txt",String)
    emails[threademail[1]][threademail[2]]["body"] = body
    sender = parse(Int64,shortinfo[1])
    recipients = parse.(Int64,split(shortinfo[2],","))
    recipients = convert(Vector{Any},recipients)
    if length(shortinfo[3]) > 0
        cc = parse.(Int64,split(shortinfo[3],","))
    else 
        cc = Vector{Any}()
    end
    cc = convert(Vector{Any},cc)
    if length(shortinfo) == 5
        subject = shortinfo[5]
    else
        subject = ""
    end
    timestamp = shortinfo[4]
    emails[threademail[1]][threademail[2]]["recipients"] = recipients
    emails[threademail[1]][threademail[2]]["time"] = timestamp
    emails[threademail[1]][threademail[2]]["cc"] = cc
    emails[threademail[1]][threademail[2]]["sender"] = sender
    emails[threademail[1]][threademail[2]]["subject"] = subject
end


data["emails"] = emails
data["clusters"] = clusters
data["names"] = names
stringdata = JSON.json(data)

# write the file with the stringdata variable information
open("../fauci-email-data.json", "w") do f
        write(f, stringdata)
end