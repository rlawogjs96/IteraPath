rule [
   ruleID "BGC1000001_m1_to_m2__r2"
   left [
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "-" ]
      node [ id 9 label "*" ]
      node [ id 1 label "C" ]
      node [ id 2 label "C" ]
      node [ id 3 label "O" ]
      node [ id 8 label "*" ]
      node [ id 10 label "*" ]
   ]
   context [
   ]
   right [
      edge [ source 9 target 1 label "-" ]
      edge [ source 9 target 10 label "-" ]
      edge [ source 9 target 8 label "-" ]
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "=" ]
      node [ id 9 label "C" ]
      node [ id 1 label "C" ]
      node [ id 2 label "C" ]
      node [ id 3 label "O" ]
      node [ id 8 label "C" ]
      node [ id 10 label "O" ]
   ]
]