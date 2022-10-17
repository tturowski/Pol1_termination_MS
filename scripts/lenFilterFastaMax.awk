#!/usr/bin/awk -f

# use awk -v var=x file_name
# where x is max length
# scaffold downloaded from https://stackoverflow.com/questions/62066578/filter-out-fasta-files-by-specified-sequence-length-in-bash

/^>/{
  if(sign_val && strLen<=var){
    print sign_val ORS line
  }
  strLen=line=""
  sign_val=$0
  next
}
{
  strLen+=length($0)
  line=(line?line ORS:"")$0
}
END{
  if(sign_val && strLen<=var){
    print sign_val ORS line
  }
}