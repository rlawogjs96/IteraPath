rule [
   ruleID "BGC1000004_m1_to_m2__r3"
   left [
      node [ id 8 label "*" ]
      node [ id 1 label "C" ]
      node [ id 10 label "*" ]
      node [ id 9 label "*" ]
   ]
   context [
   ]
   right [
      edge [ source 8 target 9 label "-" ]
      edge [ source 1 target 9 label "-" ]
      edge [ source 10 target 9 label "=" ]
      node [ id 8 label "C" ]
      node [ id 1 label "C" ]
      node [ id 10 label "O" ]
      node [ id 9 label "C" ]
   ]
]