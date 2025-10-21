rule [
   ruleID "BGC1000003_m1_to_m2__r0"
   left [
      edge [ source 5 target 6 label "=" ]
      edge [ source 5 target 7 label "-" ]
      edge [ source 9 target 10 label "-" ]
      edge [ source 9 target 11 label "-" ]
   ]
   context [
      node [ id 5 label "C" ]
      node [ id 6 label "O" ]
      node [ id 7 label "S" ]
      node [ id 9 label "C" ]
      node [ id 10 label "O" ]
      node [ id 11 label "C" ]
      edge [ source 5 target 7 label "-" ]
      edge [ source 5 target 11 label "-" ]
      edge [ source 5 target 6 label "=" ]
      edge [ source 9 target 10 label "-" ]
      edge [ source 9 target 11 label "-" ]
   ]
   right [
      edge [ source 5 target 6 label "-" ]
      edge [ source 5 target 11 label "-" ]
      edge [ source 9 target 10 label "=" ]
   ]
]