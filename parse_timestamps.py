from datetime import datetime

def parse_timestamp(ts):
    # Super special cases
    ts = ts.replace("-----      (b)( 6)   (NIH/N IA ID) [E)Sat, 28 Mar 2020 19 :43:52 +0000",
                    "Sat, 28 Mar 2020 19:43:52 +0000")
    ts = ts.replace("Saturday, February 29, 2020 8_:5_8_P----------\n_M_____           ~~\n(b)(",
                    "Saturday, February 29, 2020 8:58 PM")
    ts = ts.replace("Saturday, March 14, 2020 9:04 A\"'M;..;;.....\n____       ~~",
                    "Saturday, March 14, 2020 9:04 AM")
    ts = ts.replace("Wednesday, March 4, 2020 12:l..9-'-P.;;...\nM.;;......\n___       ~~",
                    "Wednesday, March 4, 2020 12:l9 PM")
    ts = ts.replace("Thursday, February               PM    ----------20, 2020 4 :16\n(b)(",
                    "Thursday, February 20, 2020 4:16PM")
    ts = ts.replace("Fri, 2 8 Fe b 2020 1 1:21 :22 -0500", "Fri, 28 Feb 2020 11:21:22 -0500")
    ts = ts.replace("Monda y, Febrnary l 7, 2020 3: 15 PM", "Monday, February 17, 2020 3:15 PM")
    ts = ts.replace("Tuesday, February 25 , 2020 9:0 l PM", "Tuesday, February 25, 2020 9:01 PM")
    ts = ts.replace("Sat, 28 Mar 20201 7:17:12 +0000", "Sat, 28 Mar 2020 17:17:12 +0000")
    ts = ts.replace("Monday, February 10, 2020 11: I I AM", "Monday, February 10, 2020 11:11 AM")
    ts = ts.replace("Wednesday, February 26, 2020 I :43 PM", "Wednesday, February 26, 2020 1:43 PM")
    ts = ts.replace("Friday, April l 0, 2020 6:0 l AM", "Friday, April 10, 2020 6:01 AM")
    ts = ts.replace("Sunday, April I 9, 202 0 5: 18 PM", "Sunday, April 19, 2020 5:18 PM")
    ts = ts.replace("Tuesday, February l l , 2020 1:19 PM", "Tuesday, February 11, 2020 1:19 PM")
    ts = ts.replace("Monday, February I 0, 2020 2: 17 PM", "Monday, February 10, 2020 2:17 PM")
    ts = ts.replace("Tuesday, February 11, 2020 8:0 I AM", "Tuesday, February 11, 2020 8:01 AM")
    ts = ts.replace("Thursday ,\nDecember     19 ,\n20 1 9 4 : 24 PM", "Thursday, December 19, 2019 4:24 PM")
    ts = ts.replace("Tuesday , Apri l 2 l, 2020 5 :06 PM", "Tuesday, April 21, 2020 5:06 PM")
    ts = ts.replace("18Febmary202015:                26", "18 February 2020 15:26")
    ts = ts.replace("k:,-o-n-da_y_, -A-p-ri-11- 3-,-2-0-20\n- 7-:S_S_A_M\n__",
                    "Monday, April 13, 2020 7:55 AM")    

    ts = ts.replace("\n", " ")
    ts = ts.strip('=-~-~ ')
    # Cases where we couldn't really find a timestamp anyway
    if len(ts) == 0: return ''
    
    # Thursday, February 27, 2020 I 0:42 AM
    ts = ts.replace(" I 0:", " 10:").replace(" I ", " 1 ").replace(" I,", " 1,")
    ts = ts.replace(" l,", " 1,").replace("l9", "19").replace("J 1", "11").replace("l J", "11")
    ts = ts.replace("January 3 1", "January 31")

    ts = ts.replace("-----      (b)( 6)   (NIH/N IA ID) [E)", "")
    ts = ts.replace("-------------", "").replace("--------", "").replace("-------", "")
    ts = ts.replace("---", "").replace("--", " ")

    while ts.find("  ") != -1: ts = ts.replace("  ", " ")
    for s in ["Atvl", "Arvl", "Alvl", "At\'vl", "Afvl", "Al", "AJ"]:
        ts = ts.replace(s, "AM")
    ts = ts.replace("_", "").replace("- (b)(6)>", "").replace("(b) (6)· ======~-: :-::-:::::===~~ :::::", "")
    ts = ts.strip()

    # Sunday. March 29. 2020 8:53 PM
    # Sunday , March 1, 2020 3:00 PM
    ts = ts.replace(".", ",").replace(" ,", ",")
    
    # Wednesday,March 4, 2020 8:56 PM
    ts = ts.replace("y,M", "y, M")
    
    # Wed, 4 Mar 2020 02 :55:36 +0000
    # Tuesday, March 3, 2020 9 :48 PM
    # Saturday, February 29, 2020 7: 14 PM
    # Mon, 9 Mar 2020 13:49:31-0400
    for i in range(0, 10):
        for j in range(0, 10):
            ts = ts.replace(repr(i) + " :" + repr(j), repr(i) + ":" + repr(j))
            ts = ts.replace(repr(i) + ": " + repr(j), repr(i) + ":" + repr(j))
            ts = ts.replace(":" + repr(i) + " " + repr(j), ":" + repr(i) + repr(j))
            ts = ts.replace(repr(i) + "-" + repr(j), repr(i) + " -" + repr(j))
            ts = ts.replace(repr(i) + "+" + repr(j), repr(i) + " +" + repr(j))
            ts = ts.replace(repr(i) + "," + repr(j), repr(i) + ", " + repr(j))
            
    # Thur sday, Febru ary 27, 2020 7: 16 PM
    # Thu rsday, February 20, 2020 5:29 PM
    # Tues day, March 3, 2020 12:07 AM
    # Mo n, 4 May 2020 21:05:39 +0000
    ts = ts.replace("Mo n", "Mon") .replace("M on", "Mon")
    ts = ts.replace("da y", "day").replace("s day", "sday").replace("sd ay", "sday")
    ts = ts.replace("Tu es", "Tues").replace("ue sd", "uesd")
    ts = ts.replace("Thu rsday", "Thursday").replace("Thur sday", "Thursday")
    ts = ts.replace("W ed", "Wed").replace("Wed nesday", "Wednesday").replace("Wedne sday", "Wednesday")
    ts = ts.replace("Wedne:sday", "Wednesday")
    ts = ts.replace("F1;day", "Friday")
    ts = ts.replace("Sanrrday", "Saturday").replace("Saturd ay", "Saturday").replace("ur da", "urda")

    # Fri, 2 8 Fe b 2020 1 1:21 :22 -0500
    # Thursday, Febru ary 27, 2020 2:51 PM
    # Saturday , Febrnary 29, 2020 2:25 PM
    # Monday, Ma rch 2, 2020 5:29 AM
    # Monday , Februar 03, 2020 12:21 PM
    ts = ts.replace("Fe b", "Feb").replace("b r", "br").replace("Februar 0", "February 0")
    for s in ["f ebr1.1a1y", "Febrnary", "Febru ary", "Februaiy", "Februar y",
             "Fcbniary", "Februa ry", "f ebr1,1a1y", "febrnary", "f ebr1, 1a1y"]:
        ts = ts.replace(s, "February")
    ts = ts.replace("Ma r", "Mar").replace("M ar", "Mar")
    for s in ["Marc h", "Mar ch", "Marcb"]: ts = ts.replace(s, "March")
    ts = ts.replace("Ap r", "Apr")

    # Thu, 27 Feb2020 04:44:07 +0000
    ts = ts.replace("Feb2", "Feb 2").replace("March2", "March 2").replace("April13", "April 13")
    ts = ts.replace("n,2", "n, 2").replace("y,J", "y, J")

    # Tuesday, February4, 2020 10:24 PM
    for i in range(0, 10): ts = ts.replace("y" + repr(i), "y " + repr(i))

    # Wednesday,February 19, 202011 :04 PM
    ts = ts.replace("y,F", "y, F").replace("y,A", "y, A").replace("y,2", "y, 2").replace("Thu,2", "Thu, 2")

    # Wednesday,February 19, 202011 :04 PM
    for i in range(0, 10): ts = ts.replace("2020" + repr(i), "2020 " + repr(i))

    # Friday, M arch 27, 202 0 1:29 PM
    # Thu, 19 Mar 2020 0 1:53:32 +0000
    ts = ts.replace("2020 0 ", "2020 0").replace("202 0", "2020").replace(" 20 20 ", " 2020 ").replace("2020,", "2020")
    ts = ts.replace(" 20 17 ", " 2017 ")

    # Tue, 10 Mar 2020 02:11:06 +000 0
    for i in range(0, 10):
        ts = ts.replace("+0" + repr(i) + "0 0", "+0" + repr(i) + "00")
        ts = ts.replace("-0" + repr(i) + "0 0", "-0" + repr(i) + "00")
    ts = ts.replace("+GOOO", "+0000")
    ts = ts.replace("-0 400", "-0400").replace("4- 0400", "4 -0400").replace("+00 00", "+0000")

    ts = ts.replace("P  M", "PM").replace("P M", "PM")
    if ts[-1] == "A": ts = ts + "M"

    for i in range(0, 10):
        ts = ts.replace(repr(i) + "PM", repr(i) + " PM").replace(repr(i) + "AM", repr(i) + " AM")

    pm_ind = ts.find(" PM")
    if pm_ind != -1: ts = ts[:pm_ind + 3]
    am_ind = ts.find(" AM")
    if am_ind != -1: ts = ts[:am_ind + 3]
    zone_ind = ts.find("+0000")
    if zone_ind != -1: ts = ts[:zone_ind + 5]

    # Fri, 6 Mar 2020 03:49:45 +0000
    try:    return datetime.strptime(ts, "%a, %d %b %Y %H:%M:%S %z")
    except: pass
    # Thursday, March 5, 2020 9:53 AM
    try:    return datetime.strptime(ts, "%A, %B %d, %Y %H:%M %p")
    except: pass
    # Tuesday, March 3, 2020 11:16:51 AM
    try: return datetime.strptime(ts, "%A, %B %d, %Y %H:%M:%S %p")
    except: pass
    # Wed, 4 Mar 2020 02 :55:36 +0000
    try:    return datetime.strptime(ts, "%a, %d %b %Y %H:%M:%S %z")
    except: pass
    # Tuesday, 25 February 2020 10:18
    try:    return datetime.strptime(ts, "%A, %d %B %Y %H:%M")
    except: pass
    # Monday, February 17, 2020
    try:    return datetime.strptime(ts, "%A, %B %d, %Y")
    except: pass
    # 09 February 2020 16:15
    try:    return datetime.strptime(ts, "%d %B %Y %H:%M")
    except: pass
    # Wednesday, February 5, 2020 9:17
    try:    return datetime.strptime(ts, "%A, %B %d, %Y %H:%M")
    except: pass
    # Sunday, 26 January 2020 8:59 PM
    try:    return datetime.strptime(ts, "%A, %d %B %Y %H:%M %p")
    except: pass
    # Monday, February 17, 2020 3:15 PM
    try:    return datetime.strptime(ts, "%A, %B %d, %Y %H:%M %p")
    except: pass
    
    # ... and give up on the last ~10 cases
    return ''
