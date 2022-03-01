#!/bin/bash -ue
cat DBA_2J.H1.mm39.bed DBA_2J.H2.mm39.bed | awk '!x[$0]++' > DBA_2J.DEL.validated.bed
