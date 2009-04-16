#!/usr/bin/awk -f
#
# Prints every periodth field starting from offset.
#
# Usage: every-xth-field.awk -v period=periodvalue -v offset=offsetvalue < file
#
{
    for(i=offset; i<=NF; i=i+1) 
    {
	if(0==(i-offset)%period)
        {
	    printf("%s ", $i);
        };
    }; 
    printf("\n")
}

