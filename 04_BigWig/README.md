for f in *Rpa190HTP_wt_none*sam; do cat $f | awk '$3=="chrXII"{print}' | awk '($4>=460847) && ($4<=460849) && $6!~"S" && $6!~"N" && $6!~"D" {print}' | cut -f6 | sed 's/M//' | histogram.pl | sort -k1,1n > ${f%.sam}"_RDN37-1_peak.list" & done

for f in *Rpa190HTP_wt_none*sam; do cat $f | awk '$3=="chrXII"{print}' | awk '($4>=460712) && ($4<=467569) && $6!~"S" && $6!~"N" && $6!~"D" && $6!~"I"{print}' | cut -f6 | sed 's/M//' | histogram.pl | sort -k1,1n > ${f%.sam}"_RDN37-1_all.list" & done

head -100 *list
