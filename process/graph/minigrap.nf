#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def strain_index(strain) {
    $ASSEMDIR/129S1_SvImJ_chromosomes_MT.fasta $ASSEMDIR/AKR_J_chromosomes_MT.fasta $ASSEMDIR/A_J_chromosomes_MT.fasta $ASSEMDIR/BALB_cJ_chromosomes_MT.fasta $ASSEMDIR/C3H_HeJ_chromosomes_MT.fasta $ASSEMDIR/C57BL_6NJ_chromosomes_MT.fasta $ASSEMDIR/CAST_EiJ_chromosomes_MT.fasta $ASSEMDIR/CBA_J_chromosomes_MT.fasta $ASSEMDIR/DBA_2J_chromosomes_MT.fasta $ASSEMDIR/FVB_NJ_chromosomes_MT.fasta $ASSEMDIR/JF1_MsJ_chromosomes_MT.fasta $ASSEMDIR/NOD_ShiLtJ_chromosomes_MT.fasta $ASSEMDIR/PWK_PhJ_chromosomes_MT.fasta $ASSEMDIR/SPRET_EiJ_chromosomes_MT.fasta $ASSEMDIR/WSB_EiJ_chromosomes_MT.fasta
}