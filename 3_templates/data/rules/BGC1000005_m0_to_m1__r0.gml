rule [
   ruleID "BGC1000005_m0_to_m1"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "=" ]
      edge [ source 2 target 7 label "-" ]
   ]
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "C" ]
      node [ id 3 label "O" ]
      node [ id 4 label "*" ]
      node [ id 5 label "*" ]
      node [ id 6 label "*" ]
      node [ id 7 label "S" ]
   ]
   right [
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "=" ]
      edge [ source 2 target 4 label "-" ]
      edge [ source 4 target 5 label "-" ]
      edge [ source 5 target 6 label "=" ]
      edge [ source 5 target 7 label "-" ]
   ]
]