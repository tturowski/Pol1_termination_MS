#!/usr/bin/awk -f
BEGIN {
OFS="\t"
}
NR%4 == 2 {
	lengths[length($0)]++
}
END {
	for (l in lengths) {
		print l, lengths[l]
}
}
