#!/bin/bash -ue
cat AKR_J.H6.mm39.bed AKR_J.H7.mm39.bed | awk '!x[$0]++' > AKR_J.INS.validated.bed
