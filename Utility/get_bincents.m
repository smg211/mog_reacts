function bincents = get_bincents(binedges)
if all(diff(diff(binedges)) == 0)
  bincents = binedges(2:end)-0.5*abs(diff(binedges(1:2)));
else
  binedges_diff = diff(binedges);
  bincents = binedges(1:end-1)+0.5*binedges_diff;
end