function [phiV, vbad] = volt_penalty(V, band)
below = max(0, band(1) - V); above = max(0, V - band(2));
phiV = sum(below + above);
vbad = any(below>0 | above>0);
end
