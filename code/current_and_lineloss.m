function [Iabs,Iang,Ploss,Qloss] = current_and_lineloss(Vm_i,ang_i_deg,Vm_j,ang_j_deg,Z,Sbase)
Vi = Vm_i * exp(1i*deg2rad(ang_i_deg));
Vj = Vm_j * exp(1i*deg2rad(ang_j_deg));
I  = (Vi - Vj) / Z;
Iabs = abs(I); Iang = rad2deg(angle(I));
S_loss_pu = (Iabs^2) * Z;      % I^2 * Z
Ploss = Sbase * real(S_loss_pu);
Qloss = Sbase * imag(S_loss_pu);
end
