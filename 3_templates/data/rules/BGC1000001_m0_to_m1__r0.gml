rule [
   ruleID "BGC1000001_m0_to_m1"
   left [
      edge [ source 1 target 5 label "-" ]
      edge [ source 5 target 3 label "=" ]
      edge [ source 5 target 7 label "-" ]
   ]
   context [
      node [ id 1 label "C" ]
      node [ id 5 label "C" ]
      node [ id 3 label "O" ]
      node [ id 7 label "S" ]
      node [ id 2 label "*" ]
      node [ id 4 label "*" ]
      node [ id 6 label "*" ]
      node [ id 8 label "*" ]
   ]
   right [
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "-" ]
      edge [ source 2 target 4 label "-" ]
      edge [ source 3 target 8 label "-" ]
      edge [ source 4 target 5 label "-" ]
      edge [ source 5 target 6 label "=" ]
      edge [ source 5 target 7 label "-" ]
   ]
]