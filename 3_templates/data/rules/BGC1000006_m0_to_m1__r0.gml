rule [
   ruleID "BGC1000006_m0_to_m1__r0"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "=" ]
      edge [ source 2 target 4 label "-" ]
      node [ id 5 label "*" ]
      node [ id 6 label "*" ]
      node [ id 7 label "*" ]
      node [ id 1 label "C" ]
   ]
   context [
      node [ id 2 label "C" ]
      node [ id 3 label "O" ]
      node [ id 4 label "S" ]
   ]
   right [
      edge [ source 5 target 6 label "-" ]
      edge [ source 6 target 7 label "=" ]
      edge [ source 6 target 1 label "-" ]
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "=" ]
      edge [ source 2 target 4 label "-" ]
      node [ id 5 label "C" ]
      node [ id 6 label "C" ]
      node [ id 7 label "O" ]
      node [ id 1 label "C" ]
   ]
]