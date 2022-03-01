#!/bin/bash -ue
cat DBA_2J.H6.mm39.bed DBA_2J.H7.mm39.bed | awk '!x[$0]++' > DBA_2J.INS.validated.bed
