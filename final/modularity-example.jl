##
include("../methods.jl")
##
G = _read_final_with_products("fauci-email-tofrom-5.json")
G = (G..., groups=G.products.simple.modularity)
##
drawgraph(G, size=(350,350))
showlabel!(G,"collins",7)
showlabel!(G,"conrad",7)
showlabel!(G,"giroir",7)
showlabel!(G,"redfield",7)
showlabel!(G,"farrar",7)
showlabel!(G,"niaid",7)
showlabel!(G,"routh",7)
showlabel!(G,"niaid news",7)
showlabel!(G,"kadlec",7)
showlabel!(G,"new england",7)
showlabel!(G,"giroir",7)
showlabel!(G,"birx",7)
