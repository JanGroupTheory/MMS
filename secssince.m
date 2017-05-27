function secs = secssince(datenumt)
if length(datenumt)>1
    datenumt = datenumt-datenumt(1);
end
secs = datenumt*3600*24;