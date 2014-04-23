#!/bin/tcsh

./runPDF.sh  1 0; rm -rf wwsel0j_wwmg;	    mkdir wwsel0j_wwmg;      mv pdf_*txt  wwsel0j_wwmg/;
./runPDF.sh  2 0; rm -rf wwsel0j_wwmcnlo;   mkdir wwsel0j_wwmcnlo;   mv pdf_*txt  wwsel0j_wwmcnlo/;
./runPDF.sh  3 0; rm -rf wwsel0j_wwpowheg;  mkdir wwsel0j_wwpowheg;  mv pdf_*txt  wwsel0j_wwpowheg/;
./runPDF.sh  4 0; rm -rf wwsel0j_wwggww;    mkdir wwsel0j_wwggww;    mv pdf_*txt  wwsel0j_wwggww/;
./runPDF.sh  5 0; rm -rf wwsel0j_hww;       mkdir wwsel0j_hww;       mv pdf_*txt  wwsel0j_hww/;

./runPDF.sh  1 1; rm -rf wwsel1j_wwmg;	    mkdir wwsel1j_wwmg;      mv pdf_*txt  wwsel1j_wwmg/;
./runPDF.sh  2 1; rm -rf wwsel1j_wwmcnlo;   mkdir wwsel1j_wwmcnlo;   mv pdf_*txt  wwsel1j_wwmcnlo/;
./runPDF.sh  3 1; rm -rf wwsel1j_wwpowheg;  mkdir wwsel1j_wwpowheg;  mv pdf_*txt  wwsel1j_wwpowheg/;
./runPDF.sh  4 1; rm -rf wwsel1j_wwggww;    mkdir wwsel1j_wwggww;    mv pdf_*txt  wwsel1j_wwggww/;
./runPDF.sh  5 1; rm -rf wwsel1j_hww;       mkdir wwsel1j_hww;       mv pdf_*txt  wwsel1j_hww/;

./runPDF.sh  6 0; rm -rf hwwsel0j_wwmg;	    mkdir hwwsel0j_wwmg;     mv pdf_*txt  hwwsel0j_wwmg/;
./runPDF.sh  7 0; rm -rf hwwsel0j_wwmcnlo;  mkdir hwwsel0j_wwmcnlo;  mv pdf_*txt  hwwsel0j_wwmcnlo/;
./runPDF.sh  8 0; rm -rf hwwsel0j_wwpowheg; mkdir hwwsel0j_wwpowheg; mv pdf_*txt  hwwsel0j_wwpowheg/;
./runPDF.sh  9 0; rm -rf hwwsel0j_wwggww;   mkdir hwwsel0j_wwggww;   mv pdf_*txt  hwwsel0j_wwggww/;
./runPDF.sh 10 0; rm -rf hwwsel0j_hww;      mkdir hwwsel0j_hww;      mv pdf_*txt  hwwsel0j_hww/;

./runPDF.sh  6 1; rm -rf hwwsel1j_wwmg;	    mkdir hwwsel1j_wwmg;     mv pdf_*txt  hwwsel1j_wwmg/;
./runPDF.sh  7 1; rm -rf hwwsel1j_wwmcnlo;  mkdir hwwsel1j_wwmcnlo;  mv pdf_*txt  hwwsel1j_wwmcnlo/;
./runPDF.sh  8 1; rm -rf hwwsel1j_wwpowheg; mkdir hwwsel1j_wwpowheg; mv pdf_*txt  hwwsel1j_wwpowheg/;
./runPDF.sh  9 1; rm -rf hwwsel1j_wwggww;   mkdir hwwsel1j_wwggww;   mv pdf_*txt  hwwsel1j_wwggww/;
./runPDF.sh 10 1; rm -rf hwwsel1j_hww;      mkdir hwwsel1j_hww;      mv pdf_*txt  hwwsel1j_hww/;

./runPDF.sh 11 0; rm -rf zzsel0j_zz2l;     mkdir zzsel0j_zz2l;      mv pdf_*txt  zzsel0j_zz2l/;
./runPDF.sh 12 0; rm -rf zzsel0j_wz3l;     mkdir zzsel0j_wz3l;      mv pdf_*txt  zzsel0j_wz3l/;
./runPDF.sh 13 0; rm -rf zzsel0j_zh;	   mkdir zzsel0j_zh;	    mv pdf_*txt  zzsel0j_zh/;

./runPDF.sh 11 1; rm -rf zzsel1j_zz2l;     mkdir zzsel1j_zz2l;      mv pdf_*txt  zzsel1j_zz2l/;
./runPDF.sh 12 1; rm -rf zzsel1j_wz3l;     mkdir zzsel1j_wz3l;      mv pdf_*txt  zzsel1j_wz3l/;
./runPDF.sh 13 1; rm -rf zzsel1j_zh;	   mkdir zzsel1j_zh;	    mv pdf_*txt  zzsel1j_zh/;

./runPDF.sh 14 0; rm -rf wwsssel_wwss;     mkdir wwsssel_wwss;      mv pdf_*txt  wwsssel_wwss/;
./runPDF.sh 15 0; rm -rf wwsssel_wz3l;     mkdir wwsssel_wz3l;      mv pdf_*txt  wwsssel_wz3l/;
