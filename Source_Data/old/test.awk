BEGIN { FS=OFS=", " }
NR==FNR {
    if ( NR == 1 ) {
        split($0,uniq)
    }
    for (inFldNr in uniq) {
        if ( seen[inFldNr,$inFldNr]++ ) {
            delete seen[inFldNr,$inFldNr]
            delete uniq[inFldNr]
        }
    }
    next
}
FNR==1 {
    for (inFldNr=1; inFldNr<=NF; inFldNr++) {
        if (inFldNr in uniq) {
            out2inFldNr[++numOutFlds] = inFldNr
        }
    }
}
{
    for (outFldNr=1; outFldNr<=numOutFlds; outFldNr++) {
        inFldNr = out2inFldNr[outFldNr]
        printf "%s%s", $inFldNr, (outFldNr<numOutFlds ? OFS : ORS)
    }
}
