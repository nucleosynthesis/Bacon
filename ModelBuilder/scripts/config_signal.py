import ROOT as r
directory = "category_0"
signals = {	 
	   # "th1f_WH0":  ["WH",r.kGreen	,0]       
	   #,"th1f_ZH0":  ["ZH",r.kGreen+1	,0] 
           "th1f_GGH0": ["ggH",r.kAzure+10	,0] 
           ,"th1f_VBFH0":["VBF",r.kRed	,0] 
	   }

key_order = ["Dibosons","top","W#rightarrow #mu#nu","Z#rightarrow #nu#nu"]
#key_order = ["W#rightarrow #mu#nu","Z#rightarrow #nu#nu"]

backgrounds = { 
		"top":			  [["th1f_TT","th1f_T"],		r.kBlue-2,0]
		,"Dibosons":		  [["th1f_ZZ","th1f_WW","th1f_WZ"],	r.kPink-4,    0]
		#,"Z#rightarrow ll":	  [["th1f_RDY"],			r.kYellow,  0]
		,"W#rightarrow #mu#nu":	  [["th1f_corrected_W"], 				r.kOrange-2,   0]
		#,"W#rightarrow #mu#nu":	  [["th1f_WHT","th1f_W"], 				r.kOrange-2,   0]
		#,"Z#rightarrow #nu#nu":	  [["th1f_DY"],				r.kBlue-9,  0]
		#,"Z#rightarrow #nu#nu":	  [["th1f_corrected_DY"],				r.kBlue-9,  0]
		,"Z#rightarrow #nu#nu":	  [["photon_zvv_model"],				r.kBlue-9,  0]
		#,"Z#rightarrow ll":	  [["th1f_DYLL"],				r.kBlue+5,  0]

	      }

dataname  = "th1f_data_obs"
