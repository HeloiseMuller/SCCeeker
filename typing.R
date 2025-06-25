#mec complexes:
classA <- c("mecA", "mecR1", "mecI")
classB <- c("mecA", "dmecR1", "IS1272")
#The class C mec gene complex contains mecA and truncated mecR1 by the insertion of IS431 upstream of mecA and HVR and IS431 downstream of mecA.
#Just write mecA alone as candidate class C. If class C, need to look whether mec-class-C1 or mec-class-C2 is present
classC <- c("mecA") #?
classD <- "" #not a single SCCmec with that class anyway
classE <- c("mecC", "mecR1", "mecI")

dt_class <- data.table(class = c("A", "A_rev", "B", "B_rev", "mecA_only", "E", "E_rev"),
                       comb = c(paste(classA, collapse = ":"), paste(rev(classA), collapse = ":"),
                                paste(classB, collapse = ":"), paste(rev(classB), collapse = ":"),
                                "mecA",
                                paste(classE, collapse = ":"), paste(rev(classE), collapse = ":"))
)

#ccr complexes
dt_complex <- data.table(complex = c(as.character(1:9), "1_rev", "2_rev", "3_rev", "4_rev", "7_rev", "8_rev"),
                         comb = c("ccrB1:ccrA1", "ccrB2:ccrA2", "ccrB3:ccrA3", "ccrB4:ccrA4", "ccrC1", "?", "ccrB6:ccrA1", "ccrB3:ccrA1", "ccrC2",
                                  "ccrA1:ccrB1", "ccrA2:ccrB2", "ccrA3:ccrB3", "ccrA4:ccrB4",             "ccrA1:ccrBB6", "ccrA1:ccrB3"    )
)

#Types
typing <- data.table(type = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV"),
                     comb = c("B:1", "A:2", "A:3", "B:2", "C2:5", "B:4", "5:C1", "A:4", "1:C2", "C1:7", "E:8", "9:C2", "9:A", "5:A", "A:7")
)

#Also take into account cassette assemblies from right to left
typing_rev <- data.table(type = typing$type,
                         comb =  c("1_rev:B_rev", "2_rev:A_rev", "3_rev:A_rev", "2_rev:B_rev", "5:C2", "4_rev:B_rev", "C1:5", "4_rev:A_rev", "C2:1_rev", "7_rev:C1", "8_rev:E_rev", "C2:9", "A_rev:9", "A_rev:5", "7_rev:A_rev")
)
typing <- rbind(typing, typing_rev)
