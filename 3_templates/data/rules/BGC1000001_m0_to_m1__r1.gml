rule [
   ruleID "BGC1000001_m0_to_m1__r1"
   left [
      node [ id 6 label "*" ]
      node [ id 7 label "*" ]
      node [ id 1 label "C" ]
      node [ id 5 label "*" ]
   ]
   context [
   ]
   right [
      edge [ source 6 target 7 label "-" ]
      edge [ source 6 target 5 label "-" ]
      edge [ source 6 target 1 label "-" ]
      node [ id 6 label "C" ]
      node [ id 7 label "O" ]
      node [ id 1 label "C" ]
      node [ id 5 label "C" ]
   ]
]