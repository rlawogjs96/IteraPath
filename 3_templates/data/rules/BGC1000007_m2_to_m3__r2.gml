rule [
   ruleID "BGC1000007_m2_to_m3__r2"
   left [
      edge [ source 8 target 10 label "-" ]
      node [ id 8 label "C" ]
      node [ id 10 label "S" ]
      node [ id 11 label "*" ]
      node [ id 12 label "*" ]
      node [ id 13 label "*" ]
      node [ id 14 label "*" ]
   ]
   context [
   ]
   right [
      edge [ source 8 target 11 label "-" ]
      edge [ source 11 target 12 label "-" ]
      edge [ source 12 target 14 label "-" ]
      edge [ source 12 target 13 label "=" ]
      node [ id 8 label "C" ]
      node [ id 10 label "*" ]
      node [ id 11 label "C" ]
      node [ id 12 label "C" ]
      node [ id 13 label "O" ]
      node [ id 14 label "S" ]
   ]
]