import copy
import re
from string import ascii_lowercase

# Like rstrip on a, but checks for a whole word b instead of by character
def rstrip_whole(a, b):
    if a[-len(b):] == b: return a[:-len(b)]
    return a

def format_name(name):
    s = name.lower().strip()

    s = s.replace("{", "(").replace("}", ")")
    
    for alt in ["[el", "[e]", "[e)", "(e]", "[e]\""]:
        s = s.replace(alt, "[e]")

    for alt in ["n ia id", "n iaid", "njaid", "niai0", "nip.id", "nlaid", "niai0", "niai d", 
                "niaio", "nia id", "n1a id", "niald", "njatd", "n--ia--id", "n_ia _i_d"]:
        s = s.replace(alt, "niaid")    
        
    s = s.replace("(nih/ n", "(nih/n").replace("niaid )", "niaid)")
    s = s.replace(" , ", ", ").replace(" .", ".")
        
    for alt in ["n!h", "nlh", "n ih"]: s = s.replace(alt, "nih")
    for alt in ["m .d.", "m.d ."]: s = s.replace(alt, "m.d.") 
        
    s = s.replace("[e]", "")
    
    for alt in ["f auci,", "f.aucl,", "fa uci,", "fa11ci,", "fauc i,", "faud,"]:
        s = s.replace(alt, "fauci,")
                
    s = s.replace('\'', '').replace(". nih", ".nih")
    
    for extra in ["eop/ovp", "eop/who", "eop/nsc", "edp/who"]: s = s.replace(extra, "")
        
    if s.find("on behalf of") == -1:
        left, right = s.find("<"), s.find(">")
        if left == -1 and right != -1: s = s.replace(">", "").strip()
        if left != -1 and right != -1: s = s[:left].strip()
        
        left, right = s.find("("), s.find(")")
        if left != -1 and right != -1: s = s[:left].strip()        
            
    for artifact in ['=', '~', '~', '•', '!', '\"', '\'', ':', '-', '_', '·']:
        s = s.replace(artifact, '')

    # remove various artifacts
    s = s.strip()
    extras = ['  cb', '  cbh6)', '  cbh6j', " ')", ' ( ni_h', '( ni_h', ' (b', ' (bh', ' (n', 
              ' (nih/', ' (nih/niaid', ' /os/aspr/10)', ' <f', ' cb)', ' cbh', ' cbh6l', ' h 6)', 
              ' nih/niaid)', '!nih/od)', '(n', '(nih/', ', m.d', ', m.d.', ', m.d., ph.d.', ',m.d',
              ', ph.d.', ',m.d.', ',m.d.,m.p.h.', ',m.d.,ph.d.', '- global', '.....,,.,=-e', 
              '<0d -', '<od-', '=-cc-c', '[jjcus)', '[jjcus]', '[jrdbe]', '[jrdbe)', '[mailto]', 
              '[mailto:,', 'nih/n it1t', '333333', '=-cc-c', '(nih/:-:-', '( ni_h.:..', 'cb)~', '(n ~',
              ' ,ot', '....,..,...,....e', '(bh6l', ' r;;', 'cb    ) c6', 'l. h/7n', ' >',
             ' md, mph', ',m.d.,m.p.h. cb)c', ' ccc', ' (os·', '(b,.... c', 'cb)c6', '7 7 7']
    for extra in extras: s = rstrip_whole(s, extra)
    s = s.strip().strip(" ()-,._<-")
    for extra in extras: s = rstrip_whole(s, extra)
    s = s.strip().strip(" ()-,._<-")
    
    # Remove double spaces
    while s.find("  ") != -1: s = s.replace("  ", " ")

    s = s.strip()
    return s

# regex that matches a name with any spaces between letters
def spaced_regex(name):
    return re.compile(".*" + "[-,;\s]*".join(list(name)) + ".*")

# if the name appears in the original, use the replacement name
def match_one(orig_name, name1):
    regex1 = spaced_regex(name1)
    return regex1.match(orig_name)

# if both names appear in the original, use the replacement name
def match_two(orig_name, name1, name2):
    regex1 = spaced_regex(name1)
    regex2 = spaced_regex(name2)
    return regex1.match(orig_name) and regex2.match(orig_name)

def canonicalize_name(name):
    orig = copy.deepcopy(name)
    s = format_name(name).strip()

    s = s.lower().replace("\n", " ")
    s = s.replace("at.ithin d oss", "auchincloss")
    s = s.replace("facui", "fauci")
    s = s.replace("pauc i", "fauci")
    s = s.replace("faud,", "fauci")
    s = s.replace("b~rasch", "barasch")
    s = s.replace("0dam", "odam")
    s = s.replace("ljst", "list")
    s = s.replace("ooepel", "doepel")
    s = s.replace("c -on_r_a d", "conrad")
    s = s.replace("0cgr", "ocgr")
    s = s.replace("e 6).. is in.. ger", "eisinger")
    s = s.replace("godai", "godal")
    s = s.replace("tane", "lane")
    s = s.replace("m ichael", "michael")
    s = s.replace("wholley", "whalley")
    s = s.replace("da yid", "david")
    s = s.replace("arantfno", "tarantino")
    s = s.replace("bemhardf", "bernard f")

    if s.find('@') != -1:
        s = s.replace(" ", "")

    for variant in ["fo iker", "foiker", "folkor"]:
        s = s.replace(variant, "folker")
    for variant in ["marl efred", "m. ar iefred", "marl efred", "martefred", "marl efred", "martefred"]:
        s = s.replace(variant, "mariefred")
    
    if s.find("on behalf of") != -1:
        s1, s2 = s.split("on behalf of")
        return " on behalf of ".join([canonicalize_name(s1), canonicalize_name(s2)])
        
    # Give up on these
    if s in ["$", "*@", "00", "1iij", "\jdpqhliclnquiries", "abcdefghijk",
             "cb", "cbh", "cbh6l", "ru.:.,pyo", "ww", "fm/§", "0woropa tcrw1c11"
             "olson. janet e.• ph.d.: cicek. m ine: olivier elemento : thibodeau. steve : gebo . kelly (nih/od) [gl"
             "b", "el", "eli", "dalt222", "ma", "nl", "peak",
             "nia d)... [.... el / ia 1 n..  d ) [. e]"]:
        return ''

    if len([ch for ch in list(s) if ch in ascii_lowercase]) < 4:
        return ''
    
    # Checks for certain people
    if match_one(s, "fauci"): return "fauci, anthony"

    if match_one(s, "auchinc"): return  "auchincloss, hugh"
    if match_one(s, "awwad"): return "awwad, david"
    if match_one(s, "barasch"): return "barasch, kimberly"
    if match_one(s, "erbelding"): return "erbelding, emily"
    if match_one(s, "evaklee"): return "eva k lee"
    if match_one(s, "folkers"): return "folkers, greg"
    if match_one(s, "niaidodam"): return "niaid odam"
    if match_one(s, "eisinger"): return "eisinger, robert"
    if match_one(s, "farrar"): return "farrar, jeremey"
    if match_one(s, "zhuchen"): return "chen, zhu"
    if match_one(s, "minhas"): return "minhas, raman"

    if match_two(s, "lane", "cliff"): return "lane, cliff"
    if match_two(s, "goldman", "lee"): return "goldman, lee"
    if match_two(s, "gilles", "sharon"): return "gilles, sharon"
    if match_two(s, "allen", "johnson"): return "johnson, richard allen"
    if match_two(s, "anderson", "liz"): return "anderson, liz"
    if match_two(s, "angela", "dunn"): return "dunn, angela"
    if match_two(s, "austin", "christopher"): return "austin, christopher"
    if match_two(s, "barillas", "car"): return "barillas, carolina"
    if match_two(s, "bei", "john"): return "beigel, john"
    if match_two(s, "belding", "emily"): return "erbelding, emily"
    if match_two(s, "bill", "court"): return "billet, courtney"
    if match_two(s, "birnbaum", "linda"): return "birnbaum, linda"
    if match_two(s, "birx", "deb"): return "birx, deborah"
    if match_two(s, "blonda", "quirico"): return "blonda quirico"
    if match_two(s, "bonds", "michelle"): return "bonds, michelle"
    if match_two(s, "bord", "kristin"): return "debord, kristin"
    if match_two(s, "bori", "uciana"): return "borio, luciana"
    if match_two(s, "burklow", "john"): return "burklow, john"
    if match_two(s, "callahan", "michael"): return "callahan, michael"
    if match_two(s, "caneva", "duane"): return "caneva, duane"
    if match_two(s, "carter", "mecher"): return "mecher, carter"
    if match_two(s, "cassetti", "cri"): return "cassetti, cristina"
    if match_two(s, "cetron", "marty"): return "cetron, marty"
    if match_two(s, "collin", "franc"): return "collins, francis"
    if match_two(s, "connie", "halkovich"): return "halkovich, connie"
    if match_two(s, "conra", "patr"): return "conrad, patricia"
    if match_two(s, "cox", "paul"): return "cox, paul"
    if match_two(s, "david", "gruber"): return "gruber, david"
    if match_two(s, "david", "marcozzi"): return "marcozzi, david"
    if match_two(s, "dean", "charity"): return "dean, charity"
    if match_two(s, "deatrick", "elizabeth"): return "deatrick, elizabeth"
    if match_two(s, "director", "general"): return "office of the director general"
    if match_two(s, "disbrow", "gary"): return "disbrow, gary"
    if match_two(s, "dodg", "daniel"): return "dodgen, daniel"
    if match_two(s, "donald", "eric"): return "mcdonald, eric"
    if match_two(s, "eastman", "alexander"): return "eastman, alexander"
    if match_two(s, "evans", "mariefred"): return "evans, mariefred"
    if match_two(s, "fall", "brahima"): return "fall, ibrahima soce"
    if match_two(s, "fraissy", "jean"): return "delfraissy, jean-francois"
    if match_two(s, "franken", "bob"): return "bob franken"
    if match_two(s, "gifford", "patrice"): return "allen-gifford, patrice"
    if match_two(s, "giroir", "brett"): return "giroir, brett"
    if match_two(s, "gottesman", "michael"): return "gottesman, michael"
    if match_two(s, "hahn", "stephen"): return "hahn, stephen"
    if match_two(s, "hamel", "joseph"): return "hamel, joseph"
    if match_two(s, "hanfling", "dan"): return "hanfling, dan"
    if match_two(s, "harrison", "brian"): return "harrison, brian"
    if match_two(s, "harvey", "melissa"): return "harvey, melissa"
    if match_two(s, "haskins", "melinda"): return "haskins, melinda"
    if match_two(s, "hassell", "david"): return "hassell, david"
    if match_two(s, "henning", "alexa"): return "henning, alexa"
    if match_two(s, "hepburn", "matthew"): return "hepburn, matthew"
    if match_two(s, "jeremy", "farrar"): return "jeremy farrar"
    if match_two(s, "john", "lauerman"): return "lauerman, john"
    if match_two(s, "johnson", "martin"): return "johnson, martin"
    if match_two(s, "johnson", "robert"): return "johnson, robert"
    if match_two(s, "kadlec", "robert"): return "kadlec, robert"
    if match_two(s, "kaushik", "sangeeta"): return "kaushik, sangeeta"
    if match_two(s, "lapook", "jon"): return "lapook, jon"
    if match_two(s, "lauer", "michael"): return "lauer, michael"
    if match_two(s, "laurie", "doepel"): return "doepel, laurie"
    if match_two(s, "lawler", "james"): return "lawler, james"
    if match_two(s, "lerner", "andrea"): return "lerner, andrea"
    if match_two(s, "lesley", "cahillroy"): return "lesley cahill roy"
    if match_two(s, "lisa", "koonin"): return "lisa koonin"
    if match_two(s, "maria", "kerkhove"): return "van kerkhove, maria"
    if match_two(s, "marston", "hilary"): return "marston, hilary"
    if match_two(s, "martin", "gregory"): return "martin, gregory"
    if match_two(s, "mascola", "john"): return "mascola, john"
    if match_two(s, "mcgowan", "john"): return "mcgowan, john"
    if match_two(s, "mcneil", "donald"): return "mcneil, donald"
    if match_two(s, "media", "inquiries"): return "niaid media inquiries"
    if match_two(s, "public", "inquiries"): return "niaid public inquiries"
    if match_two(s, "niaid", "ocgrleg"): return "niaid ocgr leg"
    if match_two(s, "myles", "renate"): return "myles, renate"
    if match_two(s, "nathaniel", "upert"): return "hupert, nathaniel"
    if match_two(s, "phillips", "sally"): return "phillips, sally"
    if match_two(s, "redfield", "robert"): return "redfield, robert"
    if match_two(s, "richard", "hatchett"): return "hatchett, richard"
    if match_two(s, "rioux", "amelie"): return "rioux, amelie"
    if match_two(s, "rout", "jennifer"): return "routh, jennifer"
    if match_two(s, "stewart", "simonso"): return "simonson, stewart"
    if match_two(s, "tabak", "lawrence"): return "tabak, lawrence"
    if match_two(s, "tarantino", "david"): return "tarantino, david"
    if match_two(s, "tia", "miller@turner.com"): return "tia.miller@turner.com"
    if match_two(s, "tracey", " mcnamara"): return "mcnamara, tracey"
    if match_two(s, "wade", "david"): return "wade, david"
    if match_two(s, "walters", "william"): return "walters, william"
    if match_two(s, "wargo", "michael"): return "wargo, michael"
    if match_two(s, "wilkinson", "thomas"): return "wilkinson, thomas"
    if match_two(s, "wolfe", "herbert"): return "wolfe, herbert"
    if match_two(s, "wood", "gretchen"): return "wood, gretchen"
    if match_two(s, "yeskey", "kevin"): return "yeskey, kevin"  
    if match_two(s, "zau", "victor"): return "dzau, victor"
    if match_two(s, "wolinetz", "carrie"): return "wolinetz, carrie"
    if match_two(s, "tromberg", "bruce"): return "tromberg, bruce"
    if match_two(s, "stover", "kathy"): return "stover, kathy"
    if match_two(s, "smith", "steven"): return "smith, steven"
    if match_two(s, "pate", "muhamed"): return "pate, muhammad ali"
    if match_two(s, "katie", "miller"): return "miller, katie"
    if match_two(s, "messonnier", "nancy"): return "messonnier, nancy"
    if match_two(s, "jason", "gale"): return "gale, jason"
    if match_two(s, "lorsch", "jon"): return "lorsch, jon"
    if match_two(s, "murphy", "ryan"): return "murphy, ryan"
    if match_two(s, "fine", "amanda"): return "fine, amanda"
    if match_two(s, "shu", "bryan"): return "shuy, bryan"
    if match_two(s, "redd", "john"): return "redd, john"
    if match_two(s, "arwady", "allison"): return "arwady, allison"
    if match_two(s, "hynds", "joanna"): return "hynds, joanna"
    if match_two(s, "ryan", "michael"): return "michael, ryan"
    if match_two(s, "walsh", "elizabeth"): return "walsh, elizabeth"
    if match_two(s, "reynolds", "wayne"): return "reynolds, wayne"
    if match_two(s, "schwartlander", "ber"): return "schwartlander, bernhard"
    if match_two(s, "schuchat", "anne"): return "schuchat, anne"
    if match_two(s, "morrison", "stephen"): return "morrison, stephen"
    if match_two(s, "morens", "david"): return "morens, david"
    if match_two(s, "crawford", "chase"): return "crawford, chase"
    if match_two(s, "del rio", "carlos"): return "del rio, carlos"
    if match_two(s, "shorbaj", "farah"): return "al-shorbaji, farah"
    if match_two(s, "haynes", "barton"): return "haynes, barton"
    if match_two(s, "howard", "bauchner"): return "bauchner, howard"
    if match_two(s, "grigsby", "garrett"): return "grigsby, garrett"
    if match_two(s, "williams", "judee"): return "williams, judee"
    if match_two(s, "hallett", "adrienne"): return "hallett, adrienne"
    if match_two(s, "mahjour", "jaouad"): return "mahjour, jaouad"
    if match_two(s, "dieffenbach", "carl"): return "dieffenbach, carl"
    if match_two(s, "scolnick", "edward"): return "scolnick, edward"
    if match_two(s, "edward", "holmes"): return "holmes, edward"
    if match_two(s, "cecconi", "hunimed"): return "cecconi, hunimed"
    if match_two(s, "b", "doherty"): return "doherty, briana"
    if match_two(s, "nabel", "gary"): return "nabel, gary"
    if match_two(s, "maryjane", "walker"): return "walker, maryjane"
    if match_two(s, "mcginnis", "michael"): return "mcginnis, michael"
    if match_two(s, "mollet", "melissa"): return "mollet, melissa"
    if match_two(s, "coomes", "stephanie"): return "coomes, stephanie"
    if match_two(s, "eron", "joseph"): return "eron, joseph"
    if match_two(s, "howard", "bauchner"): return "bauchner, howard"
    if match_two(s, "jordan", "mary"): return "jordan, mary"
    if match_two(s, "pettigrew", "roderic"): return "pettigrew, roderic"
    if match_two(s, "vijayraghavan", "krishna"): return "krishnaswamy, vijayraghavan"
    if match_two(s, "ocgr", "correspond"): return "niaid ocgr correspondence"
    if match_two(s, "ilkinson", "thomas"): return "wilkinson, thomas"
    if match_two(s, "kristian", "andersen"): return "andersen, kristian"
    if match_two(s, "president", "resolve"): return "president@resolvetosavelives.org"
    if match_two(s, "john", "brooks"): return "brooks, john"
    if match_two(s, "joan", "hussey"): return "hussey, joan"
    if match_two(s, "jake", "liang"): return "liang, jake"
    if match_two(s, "walensky", "rochelle"): return "walensky, rochelle"
    if match_two(s, "pratt", "michael"): return "pratt, michael"
    if match_two(s, "carter", "mec"): return "mecher, carter"
    if match_two(s, "vazquez", "nancy"): return "vazquez-maldonado, nancy"
    if match_two(s, "omalley", "devin"): return "o'malley, devin"
    if match_two(s, "frieden", "thomas"): return "frieden, thomas"
    if match_two(s, "frieden", "tom"): return "frieden, thomas"    

    if match_two(s, "lbaden", "nejm"): s = "lbaden@nejm.org"
    
    artifacts = ['(bh6j', ' cb >', '.  cb    ) c6', ' cb>', ', ph.d.', ' (o s', '6)', '(b as',
                 'cdr usn whmo/whmu', 'cdr usnwhmo/whmu', 'ph.d. (bh6j', ', ph.d.',
                 '.. .,;', '(b,.... c', '(bt(6h', '<eli', 'l. h/7n', 'cb)c',
                 '7 7 7',
                 '>',
                 ' 8', ' 6', ' 4', ';']
    for artifact in artifacts:
        s = s.replace(artifact, '')
    s = s.strip()

    num_spaces = len([ch for ch in s if ch == ' '])
    if num_spaces == 0 and s.find(",") != -1 and s.find("@") == -1 and s.find("niaid") == -1:
        s = " ".join(s.split(","))
    if num_spaces == 1 and s.find(",") == -1 and s.find("@") == -1 and s.find("niaid") == -1:
        s1, s2 = s.split(" ")
        s = ", ".join([s2, s1])

    if s in ["conra d", "conrad"]: return "conrad, patricia"
    if s == "routh": return "routh, jennifer"
    if s == "kadlec": return "kadlec, robert"
    if s == "stover": return "stover, kathy"
    if s == "awwad": return "awwad, david"
    if s == "lepore": return "lepore, loretta"
    if s == "eisinger": return "eisinger, robert"
    if s == "gottesman": return "gottesman, michael"
    if s == "grigsby": return "grigsby, garrett"
    if s == "handley": return "handley, gray"
    if s == "haskins": return "haskins, melinda"
    if s == "wood": return "wood, gretchen"
    if s == "marston": return "marston, hilary"
    if s == "crawford": return "crawford, chase"
    if s == "lane": return "lane, cliff"
    if s == "dieffen bach": return "dieffenbach, carl"
    if s == "cindy": return "burnett, cindy"
    if s == "heather": return "janik, heather"
    if s == "herb": return "johnson, herbert e"
    if s == "gmail, david": return "pryce, david"

    # Edited Autogen
    if s == "joan hussey,": s = "hussey, joan"
    if s == "jake liang": s = "liang, jake"
    if s == "pratt, michae l": s = "pratt, michael"
    if s == "walensky rochelle": s = "walensky, rochelle"
    if s == "c, sho": s = "shoc"
    if s == "k  . lawrence": s = "kerr, lawrence"
    if s == "dela cruz, charlescb": s = "dela cruz, charles"
    if s == "shawn onea il": s = "oneail, shawn"
    if s == "lovetrue, ldonae": s = "lovetrue, idonae"
    if s == "gaarv, diane": s = "gaary, diane"
    if s == "l.hoffman, stephen" or s == "stephen l. hoffman": s = "stephen l hoffman"
    if s == "\vilkinso n, thomas": s = "wilkinson, thomas"
    if s == "abuta leb, yasmeen": s = "abutaleb, yasmeen"
    if s == "acey, tr": s = "tracey"
    if s == "adams, jero me m": s = "adams, jerome"
    if s == "adr ian h ill": s = "hill, adrian"
    if s == "akinso, woleo la": s = "akinso, woleola"
    if s == "alb ert i saverio": s = "saverio, alberti"
    if s == "alexander patterso n": s = "patterson, alexander"
    if s == "alison galvan i": s = "galvani, alison"
    if s == "amanda.seaiy@cnn.com": s = "amanda.sealy@cnn.com"
    if s == "andrew cvon eschenbach": s = "andrew c von eschenbach"
    if s == "andrew turl ey": s = "turley, andrew"
    if s == "antoniak, cynt hia": s = "antoniak, cynthia"
    if s == "aofu, g": s = "fu, gao"
    if s == "phylli s arthur": s = "arthur, phyllis"
    if s == "asha m. ge org e": s = "asha m. george"
    if s == "ies, aspadeput": s = "aspadeputies"
    if s == "augustine m. k. choi": s = "choi, augustine"
    if s == "austi n, james": s = "austin, james"
    if s == "austin, cb ristopher": s = "austin, christopher"
    if s == "azar, alex m": s = "azar, alex"
    if s == "ly nn banks": s = "banks, lynn"
    if s == "barton hay nes": s = "haynes, barton"
    if s == "becks karen": s = "becks, karen"
    if s == "bowman,  lauren": s = "bowman, lauren k"
    if s == "daniel bednar ik": s = "bednarik, daniel"
    if s == "seigel, john": s = "beigel, john"
    if s == "jonatha n bennett": s = "bennett, jonathan"
    if s == "benson, conn ie": s = "benson, connie"
    if s == "jeremyberg": s = "berg, jeremy"
    if s == "bertuzzi stefano": s = "bertuzzi, stefano"
    if s == "mart in blaser": s = "blaser, martin"
    if s == "bonner, maria k": s = "bonner, maria"
    if s == "kevin bowe n": s = "bowen, kevin"
    if s == "t im boyd": s = "boyd, tim"
    if s == "boyse, nata lie": s = "boyse, natalie"
    if s == "suzannebradley": s = "bradley, suzanne"
    if s == "brennan, patr ick": s = "brennan, patrick"
    if s == "gro brund t land": s = "brundtland, gro"
    if s == "cindy burnett": s = "burnett, cindy"
    if s == "bushar, nicho las": s = "bushar, nicholas"
    if s == "thomas cahil l": s = "cahill, thomas"
    if s == "fabrizio cantin i": s = "cantini, fabrizio"
    if s == "char les holmes": s = "holmes, charles"
    if s == "stephen chia rello": s = "chiarello, stephen"
    if s == "chris so rg": s = "sorg, chris"
    if s == "lati ka chugh": s = "chugh, latika"
    if s == "cohen, eli zabeth": s = "cohen, elizabeth"
    if s == "jon coh en": s = "cohen, jon"
    if s == "i, autotell": s = "coletti, paul"
    if s == "conover, craigl": s = "conover, craig"
    if s == "corey md, larry": s = "corey, larry"
    if s == "guidelines, covid19treatment": s = "covid19 treatment guidelines"
    if s == "daley,george q": s = "daley, george q"
    if s == "daniel gagn on": s = "gagnon, daniel"
    if s == "david r. liu": s = "liu, david"
    if s == "davidwillman": s = "willman, david"
    if s == "oeatrick, elizabeth": s = "deatrick, elizabeth"
    if s == "golett i delia": s = "delia, goletti"
    if s == "dicasimir ro, gemma": s = "dicasimirro, gemma"
    if s == "dmid wo rd nerds": s = "dmid word nerds"
    if s == "dr. josh backen": s = "dr. josh backon"
    if s == "drosten, christ ian": s = "drosten, christian"
    if s == "drury,patrick anthony": s = "drury, patrick anthony"
    if s == "duch in, jeff": s = "duchin, jeff"
    if s == "ee, s cott": s = "lee, scott"
    if s == "eidex, rachelbarw ick": s = "eidex, rachel barwick"
    if s == "elgabalawy, nadia": s = "eigabalawy, nadia"
    if s == "elisafdieh": s = "eli j. safdieh"
    if s == "eliperencevich": s = "perencevich, eli"
    if s == "eliseo perez stable": s = "perez-stable, eliseo"
    if s == "perezstable, eliseo": s = "perez-stable, eliseo"
    if s == "emanuel, ezekielj": s = "emanuel, ezekiel j"
    if s == "emil io emini": s = "emini, emilio"
    if s == "engels, thoma s": s = "engels, thomas"
    if s == "erbeldlng, emily": s = "erbelding, emily"
    if s == "eric otte sen": s = "ottesen, eric"
    if s == "ev n a s ,  mariefred": s = "evans, mariefred"
    if s == "figlio la, mike": s = "figliola, mike"
    if s == "fordbarnes, arwent hia": s = "ford-barnes, arwenthia"
    if s == "fordbarnes, arwenthia": s = "ford-barnes, arwenthia"    
    if s == "fore henrietta o": s = "henrietta, fore"
    if s == "francis, chisari.": s = "francis v. chisari"
    if s == "real.francisco": s = "francisco, real."
    if s == "the odo re friedman n": s = "friedmann, theodore"
    if s == "g lass, roger": s = "glass, roger"
    if s == "gall in, john": s = "gallin, john"
    if s == "georg e gao": s = "gao, george"
    if s == "gay le smith": s = "smith, gayle"
    if s == "m artin gelbau m": s = "gelbaum, martin"
    if s == "gelfand, jeffrey a.,m": s = "gelfand, jeffrey a"
    if s ==  "genesis regalado cbh": s = "regalado, genesis"
    if s == "go ld, jeffrey p": s = "gold, jeffrey p"
    if s == "goldne r, shannah": s = "goldner, shannah"
    if s == "gr aham, barney": s = "graham, barney"
    if s == "graaff,peter jan": s = "graaff, peter jan"
    if s == "grein thomas": s = "grein, thomas"
    if s == "grogan, josephj": s = "grogan, joseph j"
    if s == "matthias gromeie r": s = "gromeier, matthias"
    if s == "h o xie  ,james": s = "hoxie, james"
    if s == "ha rper, jill": s = "harper, jill"
    if s == "hand ley, gray": s = "handley, gray"
    if s == "hunter handsfie ld": s = "handsfield, hunter"
    if s == "hannon , emm a": s = "hannon, emma"
    if s == "harold c. slavkin": s = "slavkin, harold"
    if s == "harr is, kara": s = "harris, kara"
    if s == "hasenk rug, kim": s = "hasenkrug, kim"
    if s == "ho lland, steven": s = "holland, steven"
    if s == "j1, edward holmes": s = "holmes, edward"
    if s == "yndjoanna, h": s = "hynds, joanna"
    if s == "ilona kiekbusch ,ot": s = "kiekbusch, ilona"
    if s == "jameson,james l": s = "jameson, james l"
    if s == "janet tob ias": s = "tobias, janet"
    if s == "jennife r nieman": s = "weisman, jennifer"
    if s ==  "jernigan, daniel b": s = "jernigan, daniel"
    if s == "johnson, .f.. alfre d": s = "johnson, alfred"
    if s == "kayjohnson": s = "johnson, kay"
    if s == "rob ert jone s": s = "jones, robert"
    if s == "kabir sophia": s = "kabir, sophia"
    if s == "kadiec, robert": s = "kadlec, robert"
    if s == "kai kupferschm idt": s = "kupferschmidt, kai"
    if s == "kalil, andre c": s = "kalil, andre"
    if s == "karen adamstay lo r": s = "karen adams taylor"
    if s == "kaush ck, sangee ta": s = "kaushik, sangeeta"
    if s == "ker r, lawrence": s = "kerr, lawrence"
    if s == "robert knob ler": s = "knobler, robert"
    if s == "mar ion koopmans": s = "koopmans, marion"
    if s == "krellensteih, james": s = "krellenstein, james"
    if s == "l ipkin, ian w": s = "lipkin, ian w"
    if s == "l usso, paolo": s = "lusso, paolo"
    if s == "laura landermanga rber": s = "landermangarber, laura"
    if s ==  "lavelle, jud ith": s = "lavelle, judith"
    if s == "lee, kun lin": s = "lee, kun-lin"
    if s == "leitman, laura": s = "leifman, laura"
    if s == "lepore loretta": s = "lepore, loretta"
    if s == "love, kelly a": s = "love, kelly"
    if s == "love, kellya": s = "love, kelly"    
    if s == "m ichae l oldstone": s = "oldstone, michael"
    if s == "m itch ell, michelle": s = "mitchell, michelle"
    if s == "sabrina ma lhi": s = "malhi, sabrina"
    if s == "mango, pau l": s = "mango, paul"
    if s == "mar k zuckerberg": s = "zuckerberg, mark"
    if s == "mar ks, peter": s = "marks, peter"
    if s == "margolis, leon id": s = "margolis, leonid"
    if s == "marovic h, mary": s = "marovich, mary"
    if s == "stephanie marsha ll": s = "marshall, stephanie"
    if s == "mccance  katz, elinore": s = "mccance katz, elinore"
    if s == "mcgowan , robert": s = "mcgowan, robert"
    if s == "mcmanus, ayan na": s = "mcmanus, ayanna"
    if s == "traceymcnamara": s = "mcnamara, tracey"
    if s == "teri m cpeak": s = "mcpeak, teri"
    if s == "me llors, john w": s = "mellors, john w"
    if s == "mecl1er, carter": s = "mecher, carter"
    if s == "melissa miller ... m l": s = "miller, melissa"
    if s == "mi nhas, raman": s = "minhas, raman"
    if s == "mitchell , alexander": s = "mitchell, alexander"
    if s == "seyed mo ghadas": s = "moghadas, seyed"
    if s == "morgan oliver": s = "morgan, oliver"
    if s == "mougha lian, jen": s = "moughalian, jen"
    if s == "mu n, jenny": s = "mun, jenny"
    if s == "trevor m undel": s = "mundel, trevor"
    if s == "new england journal of medi cine": s = "new england journal of medicine"
    if s == "pizzoli, nicola": s = "nico la pizzoli"
    if s == "nih directo rs executive committee": s = "nih directors executive committee"
    if s == "oma lley, devin m": s = "omalley, devin m"
    if s == "opl inger, anne": s = "oplinger, anne"
    if s == "parikj1, purvi": s = "parikh, purvi"
    if s == "patr ick, vanessa": s = "patrick, vanessa"
    if s == "pence, laur a": s = "pence, laura"
    if s == "suzannepeskin": s = "peskin, suzanne"
    if s == "phi l sklar": s = "sklar, phil"
    if s == "poole marcia": s = "poole, marcia"
    if s == "porter, macaulay v. cb": s = "porter, macaulay v"
    if s == "prof trevorm jones cbefmedsci": s = "prof trevor m jones cbefmedsci"
    if s == "thomas qu inn": s = "quinn, thomas"
    if s == "sierpinski, radostaw": s = "rados law sierpinski"
    if s == "ramnik xav ier": s = "xavier, ramnik"
    if s == "rebek a yasmin  cepl": s = "rebeka yasmin cepi"
    if s == "ricksawaya": s = "sawaya, rick"
    if s == "ro stami, nahid": s = "rostami, nahid"
    if s == "rotro sen, daniel": s = "rotrosen, daniel"
    if s == "saukkonen, jussi j": s = "saukkonen, jussi"
    if s == "stuart schrei ber": s = "schreiber, stuart"
    if s == "seigfre id, kim": s = "seigfreid, kim"
    if s == "sharp less, norman": s = "sharpless, norman"
    if s == "shaya ce...cile": s = "shaya, cecile"
    if s == "short, marct": s = "short, marc t"
    if s == "shuy, cait rin": s = "shuy, caitrin"
    if s == "siege l, marc": s = "siegel, marc"
    if s == "skinner, james b": s = "skinner, james"
    if s == "spitaln iak, lawa": s = "spitalniak, laura"
    if s == "steck er, judy": s = "stecker, judy"
    if s == "stee le, danielle": s = "steele, danielle"
    if s == "steinberg, danie lle": s = "steinberg, danielle"
    if s == "strauss, nico le": s = "strauss, nicole"
    if s == "swami nathan, soumya": s = "swaminathan, soumya"
    if s == "t uler, matias": s = "tuler, matias"
    if s == "william temple t on": s = "templeton, william"
    if s == "teresa m iller de vega": s = "teresa miller de vega"
    if s == "thom as r. frieden": s = "frieden, thomas"
    if s == "tru eman, laura": s = "trueman, laura"
    if s == "upton, frede": s = "upton, fred"
    if s == "vanhoof, johan": s = "van hoof, johan"
    if s == "vasquez, aureli o": s = "vasquez, aurelio"
    if s == "w ha lley, david": s = "wholley, david"
    if s == "wa lensky, loren d": s = "walensky, loren d"
    if s == "walens ky, rochelle": s = "walensky, rochelle"
    if s == "walke r, robert": s = "walker, robert"
    if s == "weahkee, michael cbh6": s = "weahkee, michael"
    if s == "wol, anki": s = "wolf, anki"
    if s == "wolfe, herber": s = "wolfe, herbert"
    if s == "y ewdell, jon": s = "yewdell, jon"

    if s == "buono,lucia<lbuono@imf.orgonbehalfofgopinath,gita":
        s = "buono, lucia on behalf of gopinath, gita"
    if s == "corey md, lany": s = "corey, larry"
    if s == "da ni el lucey  (bh  cbh": s = "lucey, daniel"
    if s == "dodgen,..,r dan ie i": s = "dodgen, daniel"
    if s == "dr. eva k": s = "eva k lee"
    if s == "e in.. is ger": s = "eisinger, robert"
    if s == "glim cher, laurie,m.d. .as": s = "glimcher, laurie"
    if s == "hall, biii": s = "hall, bill"
    if s == "koopmans, m.p.g.": s = "koopmans, marion"
    if s == "mmwr media list": s = "mmwrmedia@listserv.cdc.gov"
    if s == "ocpostoffice, niaid": s = "ocpostoffice@niaid.nih.gov"
    if s == "pixton<mailto, heather": s = "pixton, heather"
    if s == "faud, anthony": s = "fauci, anthony"
    if s in ["c ,la ne", "c lif f h/ niaid", "la ,c1i ... ff",
             "clane@niaid.nih.gov"]:
        s = "lane, cliff"

    # Cleanup...
    for artifact in ['lt usn jsj4', 'c  ..........',
                     '<mailto', '[mailto', '[ma ilto',
                     '.....ccs', ', nlm/ncbi', '(nih/fic',
                     ' <1',  ',m.d. .as', '... ..,as< (b',
                     '  (h.. h /a.. s s . l',
                     '(bh6l s', 'cbh6j', 'cb  h ]  cb',
                     'cb)<6j', 'cbh6', 'cbc6', ' h6', '7',
                     'm.d. m.p.h.',
                     'nih/od',
                     '1,',
                     ',,']:
        s = s.replace(artifact, '')
    s = s.strip(' .,;')

    # ambiguous
    if s == "allen": return ''
    if s == "aspa": return ''
    if s == "cecpep6ekob": return ''
    if s == "charles": return ''
    if s == "evans": return ''
    if s =="francis": return ''
    if s =="ifflff": return ''
    if s =="jennifer": return ''
    if s =="joseph": return ''
    if s =="kathy": return ''
    if s =="kevin": return ''
    if s =="lawrence": return ''
    if s =="laura": return ''
    if s =="megan": return ''    
    if s =="richard": return ''
    if s =="robert": return ''
    if s =="ryan": return ''
    if s =="stuart": return ''
    if s == "olson. janet e. ph.d. cicek. m ine olivier elemento thibodeau. steve gebo. kelly":
        return ["olson, janet", "cicek, mine", "elemento, olivier", "thibodeau, steve", "gebo, kelly"]
    if s == "vazquez": return ''
    

    # bad splits
    if s == "kabir, sophia fares,christine youssef":
        return ["kabir, sophia", "fares, christine youssef"]
    if s == "amanda.sealy@cnn.comneel.khairzada@turner.com":
        return ["amanda.sealy@cnn.com", "neel.khairzada@turner.com"]
    if s == "thomas r. frieden lynn banks": return ["frieden, thomas", "banks, lynn"]
    if s == "david rubensteindonaho e, john": return ["rubenstein, david", "donahoe, john"]
    if s == "michael rosbash ramnik xavier": return ["xavier, ramnik", "rosbash, michael"]

    # emails
    if s == "amanda.sealy@cnn.com": return "sealy, amanda"
    if s == "amb@cbsnews.com": return "birnbaum, amy"
    if s == "hfore@unicef.org": return "henrietta, fore"
    if s == "l...i@sisrchey11.ntacod.e": s = "kircheis dr. ralf"
    if s =="lbaden@nejm.org": return "baden, lindsey"
    if s =="mtones8@myseneca.ca": s = "mtorres8@myseneca.ca"
    if s =="naugenstein@wtop.com": return "augenstein, neal"
    if s =="perencevich@uiowa.edu": return "perencevich, eli"
    if s =="pnicholas@theatlantic.com": return "nicholas, peter"
    if s == "schoofs@usc.edu": return "schoofs, mark"
    if s == "secretary@hhs.gov": return "azar, alex"
    if s == "vdzau@nas.edu": return "dzau, victor"
    if s == "ocpostoffice@niaid.nih.gov": s = "niaid ocpostoffice"

    if s == "amisimms": s = "simms, ami"
    if s == "anthony": s = "fauci, anthony"
    if s == "arthur": s = "bobrove, arthur"  # redacted last name
    if s == "aspadeputies": s = "aspa deputies"
    if s == "blatner": s = "gretta, blatner"
    if s == "caneva": s = "caneva, duane"
    if s == "carrie": s = "wolinetz, carrie"
    if s == "catherine": s = "bird, catherine"
    if s == "chris.elias" or s == "chris elias anthony": s = "elias, chris"
    if s == "courtney": s = "billet, courtney"
    if s == "daniellaleger": s = "leger, daniella"
    if s == "drury": s = "drury, patrick anthony"
    if s == "eastman": s = "eastman, alexander"
    if s == "emoryford": s = "ford, emory"
    if s == "fabien": s = "sordet, fabien" # redacted last name
    if s == "fares" or s == "fares,christine youssef": s = "fares, christine youssef"
    if s == "hatchett": s = "hatchett, richard"
    if s == "jacquelyn": s = "madry-taylor, jacquelyn"
    if s == "jayshaylor": s = "shaylor, jay"
    if s == "jinwanzhu": s = "wan-zhu, jin"
    if s == "joelmeyer": s = "meyer, joel"
    if s == "kaushik": s = "kaushik, sangeeta"
    if s == "kenglen": s = "glen, ken"
    if s == "miriamrodriguez": s = "rodriguez, miriam"
    if s == "moreno": s = "moreno, rafael"
    if s == "rafael": s = "moreno, rafael"
    if s == "richsilverman": s = "silverman, rich"
    if s == "robynsnyder": s = "snyder, robyn"
    if s == "ryanmeyer": s = "meyer, ryan"
    if s == "schreiber": s = "schreiber, stuart"
    if s == "segal": s = "segal, allen"
    if s == "shirley.gathers": s = "gathers, shirley"
    if s == "short": s = "short, marc t"
    if s == "tabak": s = "tabak, lawrence"
    if s == "tracey": s = "mcnamara, tracey"
    if s == "weisman": s = "weisman, jennifer"
    if s == "wilkinson": s = "wilkinson, thomas"
    if s == "zhijian": s = "chen, zhijian"
        
    # From dgleich
    if s == "mecher, cai1er": s = "mecher, carter"
    if s == "hiattf@washpost.com": s = "hiatt, fred"
    if s == "bdoherty@mrns.org": s = "doherty, briana"
    if s == "sheila.kaplan@nytimes.com": s = "kaplan, sheila"
    if s == "julia.belluz@voxmedia.com": s = "belluz, julia"
    if s == "jackscientek": s = "scientek, jack"
    if s == "muhamed, pate": s = "pate, muhammad ali"
    if s == "lapook, cbs": s = "lapook, jon"

    # Ones we can't figure out still...
    if s.find(".... el") != -1: return ''
    if s in [")c6jr, cb",
             "a.. s p 0 r. / 1",
             "e in.. is ger",
             "i... rc.. ra r... . .. f< k ..r",
             "ia 1 n..  d ) [. e]",
             "ipjjlli1",
             "mou., ...  n,jen",
             "ru.. pyo",
             ]:
        return ''

    if s.find("@") == -1:
        s = s.replace('.', '')

    # Hand-written
    if s == "lipkin, ian w" or s == "lipkin, lan": s = "lipkin, ian"
    if s == "marti n, gregorv j": s = "martin, gregory"
    if s == "tromberg, bruce on behalf of tromberg": s = "tromberg, bruce"
    if s == "comad, patri cia": s = "conrad, patricia"
    if s == "pel, l,, aurie": s = "doepel, laurie"
    if s == "palea, joe": s = "palca, joe"
    if s == "stephen l. hoffman": s = "stephen l hoffman"
    if s == "be rkowitz, avra hm j eopfwho" or s == "berkowitz, avrahm j": s = "berkowitz, avrahm"
    if s == "frieden, tom": s = "frieden, thomas"
    if s == "sy, elhadi": return "sy, elhadj"
    if s == "tore godal godal, tore": return "godal, tore"
    if s == "troye, nsc": return "troye, olivia"
    if s == "m ichelle whitten global": return "whitten, michelle"
    if s == "alexander morden md": return "morden, alexander"    
    if s == "lawrence 0 brown": return "brown, lawrence"
    if s == "william h sherman, md": return "sherman, william"
    if s == "mccoll um, jeffrey t": return "mccollum, jeffrey"
    if s == "ogan gurel": return "gurel, ogan"
    if s == "jeongsun seo md phd": return "seo, jeongsun"
    if s == "forde, michae l": return "forde, michael"
    if s == "imbria le, samuel": return "imbriale, samuel"
    if s == "valant ine, hannah": return "valantine, hannah"
    if s == "judith  wasserhei t": return "wasserheit, judith"
    if s == "br idbord, ken": return "bridbord, ken"
    if s == "elaine j abrams, m d": return "abrams, elaine"
    if s == "gibbo ns, gary": return "gibbons, gary"
    if s == "koroshetz, walt er": return "koroshetz, walter"
    if s == "goodcohn, mered ith": return "good-cohn, meredith"
    if s == "lankfo rd, hannah a": return "lankford, hannah"
    if s == "gordon,bruce allan": return "gordon, bruce"
    if s == "skip vi rgin": return "virgin, skip"
    if s == "char les dinarello": return "dinarello, charles"
    if s == "be th abramson": return "abramson, beth"
    #if s == "martin,  robert": return "martin, robert"
    if s == "martin,  robert": return ["kadlec, robert", "martin, gregory"]
    if s == "kontoyiannis,dimitrios p": return "kontoyiannis, dimitrios"
    if s == "stre ngthmcgaughey,  tracie": return "strength-mcgaughey, tracie"
    if s == "dr josh backon": return "backon, josh"
    if s == "m cguffee, tyler ann a": return "mcguffee, tyler ann"
    if s == "alerts, google": return "google alerts"
    if s == "aylward , raymond bruce j": return "bruce, raymond"
    if s == "gro harlem brundtland": return "brundtland, gro"
    if s == "ma rcia bache": return "bache, marcia"
    if s == "hanson, e lizabet h": return "hanson, elizabeth"
    if s == "nicho las agresti": return "agresti, nicholas"
    if s == "deepak bhatt manager communicat ions": return "bhatt, deepak"
    if s == "pablo l pena antonmarch i": return "antomarchi, pablo"
    if s == "dust i rainey": return "rainey, dusti"
    if s in ["dr art kamm", "kamm,a rt"]: return "kamm, art"
    if s == "d isbrow, gmy": return "disbrow, gary"
    if s == "avr aham halbrei ch": return "halbreich, avraham"
    if s == "balatbat, celynn e": return "balatbat, celynne" 
    if s == "nigam, m inali": return "nigam, minali"
    if s == "dan watt endorf jennifer": return "wattendorf, dan"
    if s == "wat t s, ary ee": return "watts, mary lee"
    if s == "nlald ocgr nswb": return "niaid ocgr nswb"
    if s == "labc, od": return "od labc"
    if s == "r, cieslakpaul": return "cieslak, paul"
    if s == "secretariat, hivr4p": return "hivr4p secretariat"
    if s == "d, patricia": return "conrad, patricia"
    if s == "t, bille": return "billet, courtney"
    if s == "teresa, c": return "teresa miller de vega"
    if s == "w, gary": return "disbrow, gary"
    if s == "tcrw1c1 0woropa": return "tsoli, theodora"
    if s == "kumar shah, md": return "shah, kumar"
    if s == "singer,peter alexander": return "singer, peter alexander"
    if s == "prof trevor m jones cbefmedsci": return "jones, trevor"
    if s == "oleary, brendan stenzei, timothy": return ["oleary, brendan", "stenzel, timothy"]
    if s == "sharon hillier phd": return "hillier, sharon"
    if s == "alex wolf, esq": return "wolf, alex"
    if s == "baldw in, brittany l": return "baldwin, brittany"
    if s == "multilateral, oga": return "oga multilateral"
    if s == "contacts, odolpaleg": return "odolpaleg contacts"
    if s == "ra n court": return "rancourt, anne"

    if s == "baier, bret nih000513": return "baier, bret"
    if s == "kabir, sophia nih00135": return "kabir, sophia"
    if s == "kanarek, morgan nih001052": return "kanarek, morgan"
    if s == "lewis m dmsin": return "lewis m drusin"
    if s == "muhammad ali pate": return "pate, muhammad ali"
    if s == "nih002 108 stecker, judy": return "stecker, judy"
    if s == "oplinger, anne nih0005 14": return "oplinger, anne"
    if s == "leger, daniellanih001111": return "leger, daniella"
    if s == "kabir, sophia nih001252 fares,christine youssef":
        return ["kabir, sophia", "fares, christine youssef"]
    if s == "rob erts, jacqueline": return "roberts, jacqueline"
    if s == "brasch, kir,1bt?rly": return "barasch, kimberly"
    if s == "w right, janet": return "wright, janet"
    if s == "secretariat, gpmb": return "gpmb secretariat"
    if s == "ncalio@airlines.org": return "calio, nicholas"
    if s == "kane, eiieen": return "kane, elleen"
    if s == "divs, op": return "op divs"
    if s == "gmail, wrb": return "brody, willian"
    if s == "nih000040 toomas palu": return "palu, toomas"
    if s == "nih001658, kevin": return "bowen, kevin"
    if s == "hawk ja m ar": return "hawkins, jamar"
    if s == "sec, asprexec": return "aspr exec sec"
    if s == "herrman, jack": return "herrmann, jack"
    if s == "pekoe, ken": return "pekoc, ken"
    if s == "icddirl@list.nih.gov": return "icddir-l@list.nih.gov"
    if s == "lee, kunlin": return "lee, kun-lin"
    if s == "simmon sbutler, kirk": return "simmons-butler, kirk"
    if s == "nihstaff@list.nih.gov": return "nih-staff@list.nih.gov"
    if s == "warne r, agnes": return "warner, agnes"
    if s == "peti llo, jay": return "petillo, jay"
    if s == "whalley, david": return "wholley, david"
    
    if s in ["sy, as", "m, jon", "george", "cb, david", "s sarah",
             "mich ae i j", "i, president", "fh mu hal to",
             "con t act"]:
        return ''

    # Lots of variations (could possible be fixed)
    if s == "niaid": return ''
    
    return s.strip()

def parse_name(name):
     s = canonicalize_name(name)
     return s
