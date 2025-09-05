function Pmw = pvPlantMW_datasheet(G_Wm2, TA_C, P)
s = max(0, G_Wm2/1000);                 % kW/m^2 (0..~1)
FF = (P.VMPP*P.IMPP)/(P.VOC*P.ISC);
Tcy = TA_C + s*((P.NOT-20)/0.8);
Vy  = P.VOC - P.Kv*Tcy;
Iy  = s * ( P.ISC + P.Ki*(Tcy-25) );
Po_W = P.N_modules * FF * Vy .* Iy;     % Watts
Pmw  = max(0, Po_W) / 1e6;              % MW
end
