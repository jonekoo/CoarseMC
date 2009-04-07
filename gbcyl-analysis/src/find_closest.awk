#!/usr/bin/awk -f
#
# Finds the closest number from the first column of the given input compared to
# parameter value.
#
# To find the value closest to 0.1 in the first column of datafile run:
#
# find_closest.awk -v value=0.1 datafile
#
{
    current=$1;
    diffsqr=(value-current)*(value-current);
    if(NR==1) {temp=diffsqr;}
    if(diffsqr < temp) { 
      temp=diffsqr; 
      closest=current; 
      nr_closest=NR; 
    }
}
END{
    print closest " " nr_closest;
}
