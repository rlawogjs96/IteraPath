rule [
   ruleID "BGC1000006_m3_to_m4__r0"
   left [
      edge [ source 11 target 12 label "=" ]
      edge [ source 11 target 13 label "-" ]
      edge [ source 15 target 16 label "-" ]
      edge [ source 15 target 17 label "-" ]
   ]
   context [
      node [ id 11 label "C" ]
      node [ id 12 label "O" ]
      node [ id 13 label "S" ]
      node [ id 15 label "C" ]
      node [ id 16 label "O" ]
      node [ id 17 label "C" ]
      edge [ source 11 target 13 label "-" ]
      edge [ source 11 target 12 label "=" ]
      edge [ source 11 target 17 label "-" ]
      edge [ source 15 target 16 label "-" ]
      edge [ source 15 target 17 label "-" ]
   ]
   right [
      edge [ source 11 target 12 label "-" ]
      edge [ source 11 target 17 label "-" ]
      edge [ source 15 target 16 label "=" ]
   ]
]