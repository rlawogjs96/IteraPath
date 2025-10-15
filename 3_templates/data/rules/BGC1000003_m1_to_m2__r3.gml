rule [
   ruleID "BGC1000003_m1_to_m2__r3"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "=" ]
      node [ id 8 label "*" ]
      node [ id 9 label "*" ]
      node [ id 1 label "C" ]
      node [ id 2 label "C" ]
      node [ id 3 label "O" ]
      node [ id 10 label "*" ]
   ]
   context [
   ]
   right [
      edge [ source 8 target 9 label "-" ]
      edge [ source 9 target 1 label "-" ]
      edge [ source 9 target 10 label "=" ]
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "-" ]
      node [ id 8 label "C" ]
      node [ id 9 label "C" ]
      node [ id 1 label "C" ]
      node [ id 2 label "C" ]
      node [ id 3 label "O" ]
      node [ id 10 label "O" ]
   ]
]