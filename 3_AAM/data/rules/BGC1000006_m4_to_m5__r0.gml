rule [
   ruleID "BGC1000006_m4_to_m5__r0"
   left [
      edge [ source 14 target 16 label "-" ]
      edge [ source 18 target 19 label "-" ]
      edge [ source 18 target 20 label "-" ]
   ]
   context [
      node [ id 14 label "C" ]
      node [ id 16 label "S" ]
      node [ id 18 label "C" ]
      node [ id 19 label "O" ]
      node [ id 20 label "C" ]
      edge [ source 14 target 20 label "-" ]
      edge [ source 14 target 16 label "-" ]
      edge [ source 18 target 19 label "-" ]
      edge [ source 18 target 20 label "-" ]
   ]
   right [
      edge [ source 14 target 20 label "-" ]
      edge [ source 18 target 19 label "=" ]
   ]
]