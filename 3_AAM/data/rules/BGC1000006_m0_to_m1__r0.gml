rule [
   ruleID "BGC1000006_m0_to_m1__r0"
   left [
      edge [ source 2 target 4 label "-" ]
      edge [ source 6 target 7 label "-" ]
      edge [ source 6 target 8 label "-" ]
   ]
   context [
      node [ id 2 label "C" ]
      node [ id 4 label "S" ]
      node [ id 6 label "C" ]
      node [ id 7 label "O" ]
      node [ id 8 label "C" ]
      edge [ source 2 target 8 label "-" ]
      edge [ source 2 target 4 label "-" ]
      edge [ source 6 target 7 label "-" ]
      edge [ source 6 target 8 label "-" ]
   ]
   right [
      edge [ source 2 target 8 label "-" ]
      edge [ source 6 target 7 label "=" ]
   ]
]