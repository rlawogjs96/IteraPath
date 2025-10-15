rule [
   ruleID "BGC1000000_m3_to_m4__r2"
   left [
      node [ id 16 label "*" ]
      node [ id 1 label "C" ]
      node [ id 14 label "*" ]
      node [ id 15 label "*" ]
   ]
   context [
   ]
   right [
      edge [ source 16 target 1 label "-" ]
      edge [ source 16 target 14 label "=" ]
      edge [ source 16 target 15 label "-" ]
      node [ id 16 label "C" ]
      node [ id 1 label "C" ]
      node [ id 14 label "O" ]
      node [ id 15 label "C" ]
   ]
]