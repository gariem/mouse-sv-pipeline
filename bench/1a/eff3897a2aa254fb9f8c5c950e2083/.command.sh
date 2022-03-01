#!/bin/bash -ue
cat AKR_J.H2.mm39.bed AKR_J.H1.mm39.bed | awk '!x[$0]++' > AKR_J.DEL.validated.bed
