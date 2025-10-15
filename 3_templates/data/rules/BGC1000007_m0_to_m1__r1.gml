rule [
   ruleID "BGC1000007_m0_to_m1__r1"
   left [
      node [ id 5 label "*" ]
      node [ id 6 label "*" ]
      node [ id 7 label "*" ]
      node [ id 1 label "C" ]
   ]
   context [
   ]
   right [
      edge [ source 5 target 6 label "-" ]
      edge [ source 6 target 7 label "=" ]
      edge [ source 6 target 1 label "-" ]
      node [ id 5 label "C" ]
      node [ id 6 label "C" ]
      node [ id 7 label "O" ]
      node [ id 1 label "C" ]
   ]
]