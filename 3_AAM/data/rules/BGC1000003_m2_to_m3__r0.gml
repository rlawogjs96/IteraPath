rule [
   ruleID "BGC1000003_m2_to_m3__r0"
   left [
      edge [ source 8 target 10 label "-" ]
      edge [ source 12 target 13 label "-" ]
      edge [ source 12 target 14 label "-" ]
   ]
   context [
      node [ id 8 label "C" ]
      node [ id 10 label "S" ]
      node [ id 12 label "C" ]
      node [ id 13 label "O" ]
      node [ id 14 label "C" ]
      edge [ source 8 target 14 label "-" ]
      edge [ source 8 target 10 label "-" ]
      edge [ source 12 target 14 label "-" ]
      edge [ source 12 target 13 label "-" ]
   ]
   right [
      edge [ source 8 target 14 label "-" ]
      edge [ source 12 target 13 label "=" ]
   ]
]