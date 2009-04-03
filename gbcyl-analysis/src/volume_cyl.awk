#!/usr/bin/awk -f
{
  if($0 ~ /Lz/) {
    vol = $4*$2*$2*4.0*atan2(1.0, 1.0);
    print vol;
  } 
}


     