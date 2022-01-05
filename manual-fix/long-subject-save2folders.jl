using JSON

# include("../methods.jl")
data = JSON.parsefile("../fauci-email-graph.json")


## List subjects for all emails that don't begin in "RE", "Re", "Fw", or "FW", and have a least two characters
emails = data["emails"]
for i = 1:length(emails)
  thread = emails[i]
  for j = 1:length(thread)
    email = thread[j]
    sub = email["subject"]
    if length(sub) > 1
      try
        if sub[1:2] == "RE" || sub[1:2] == "Re" || sub[1:2] == "FW" || sub[1:2] == "Fw"
          continue
        else
          print("$i, $j \t ")
          println(email["subject"])
        end
      catch
        print("$i, $j \t ")
        println(email["subject"])
      end
  end
  end
end

## Print possibly problematic long subject lines
# This list was manually put together based on subject lines that looked possibly problematic
# Each pair gives the thread number and the email-within-thread number of the problematic email
pairs = [
1257 1;
1212 1;
1169 1;
1132 2;
1088 2;
1087 2; 
745 1;
712 6;
712 5;
712 2;
706 1;
705 1;
679 1;
589 1;
366 2;
292 1;
255 2;
221 3;
214 1;
96 2;
65 2;
57 2;
52 2;
49 1;
3 2]

# Confirmed manually that these entries in "pairs" are not errors
NotProblematic = [23,22,21,18,16,13,10,9,8,6,5,4,2,1]

# These require some type of manual fix
ToFix = setdiff(1:size(pairs,1),NotProblematic)

problempairs = pairs[ToFix,:]
checkpairs = pairs[NotProblematic,:]    # just to double check
## Save these in folders, for manual editing
for k = 1:size(problempairs,1)
    names = data["names"]
    thread = problempairs[k,1]
    emailno = problempairs[k,2]
    email = emails[thread][emailno]
    mkdir("longsubjecterror-emails/$thread-$emailno")
    io = open("longsubjecterror-emails/$thread-$emailno/original-allprint.txt", "w")
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

    io = open("longsubjecterror-emails/$thread-$emailno/original-short.txt", "w")
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

    io = open("longsubjecterror-emails/$thread-$emailno/fixed-short.txt", "w")
    io = open("longsubjecterror-emails/$thread-$emailno/fixed-body.txt", "w")

end


## Check these again too, just to be extra careful

for k = 1:size(checkpairs,1)

    thread = checkpairs[k,1]
    emailno = checkpairs[k,2]
    email = emails[thread][emailno]
    mkdir("longsubject2-emails/$thread-$emailno")
    io = open("longsubject2-emails/$thread-$emailno/original-allprint.txt", "w")
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

    io = open("longsubject2-emails/$thread-$emailno/original-short.txt", "w")
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

    io = open("longsubject2-emails/$thread-$emailno/fixed-short.txt", "w")
    io = open("longsubject2-emails/$thread-$emailno/fixed-body.txt", "w")

end
