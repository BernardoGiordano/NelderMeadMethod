
function Bz=funzione(R2, I2, Z2, z)
mu0 = 4*pi*1e-7; % vacuum permeability [H/m]
 Bz = (mu0/2)*(((I2*(R2.^2))./((((Z2-z).^2)+(R2.^2)).^(3/2))) + ((I2*(R2.^2))./((((Z2+z).^2)+(R2.^2)).^(3/2))));
end